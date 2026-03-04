package sim.simulation;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.IntRange;
import sim.simulation.SARSimulator.Status;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Gillespie（総レート）方式による SAR シミュレータ。
 * <p>
 * 各ステップで全アクティブイベントの総レートを計算し、
 * 次のイベント時刻と種類を確率的に決定する。
 * 各辺は感染源の感染期間中に最大1回しか伝播しないため、
 * イベント駆動型 {@link SARSimulator} と統計的に同値なモデルとなる。
 */
public final class SARGillespieSimulator {

    private final DirectedGraph g;
    private final double lambdaDirected;
    private final double lambdaNondirected;
    private final double mu;
    private final double tMax;
    private final int[] thresholdList;
    private final Random random;

    // ── ノード状態 ──
    private final Status[] status;
    private final int[] infectedCount;
    private final double[] tInfect;
    private final double[] tRecover;

    // ── アクティブ A ノード（swap-removal 用 indexed list）──
    private final int[] aList;
    private int aSize;
    private final int[] aPos; // aPos[node] = index in aList, -1 if inactive

    // ── アクティブ有向辺（swap-removal 用 indexed list）──
    private final int[] activeDirEdges; // graph edge index
    private int activeDirCount;
    private final int[] dirPos; // dirPos[edgeIdx] = index in activeDirEdges, -1 if inactive

    // ── アクティブ無向辺（swap-removal 用 indexed list）──
    private final int[] activeUndirEdges;
    private int activeUndirCount;
    private final int[] undirPos;

    // ── 記録用 ──
    private final ArrayList<Double> times = new ArrayList<>();
    private final ArrayList<Integer> S = new ArrayList<>();
    private final ArrayList<Integer> A = new ArrayList<>();
    private final ArrayList<Integer> R = new ArrayList<>();
    private final ArrayList<Integer> Phi = new ArrayList<>();
    private int Scount, Acount, Rcount, PhiCount;
    private double initialAdoptedTime, finalAdoptedTime;

    public SARGillespieSimulator(DirectedGraph g, double lambdaDirected, double lambdaNondirected,
            double mu, double tMax, int[] thresholdList, long seed) {
        if (g == null) throw new IllegalArgumentException("Directed graph must be non-null");
        if (lambdaDirected < 0 || lambdaNondirected < 0)
            throw new IllegalArgumentException("lambdaDirected and lambdaNondirected must be non-negative");
        if (tMax <= 0) throw new IllegalArgumentException("tMax must be positive");
        if (thresholdList == null || thresholdList.length != g.n)
            throw new IllegalArgumentException("thresholdList must be an array of length n");

        this.g = g;
        this.lambdaDirected = lambdaDirected;
        this.lambdaNondirected = lambdaNondirected;
        this.mu = mu;
        this.tMax = tMax;
        this.thresholdList = thresholdList;
        this.random = new Random(seed);

        int n = g.n;
        int m = g.m;
        this.status = new Status[n];
        this.infectedCount = new int[n];
        this.tInfect = new double[n];
        this.tRecover = new double[n];

        this.aList = new int[n];
        this.aPos = new int[n];

        this.activeDirEdges = new int[m];
        this.dirPos = new int[m];

        this.activeUndirEdges = new int[m];
        this.undirPos = new int[m];
    }

    /**
     * シミュレーションを実行する。
     *
     * @param initialInfecteds 初期感染者リスト
     * @return シミュレーション結果
     */
    public SARResult run(int[] initialInfecteds) {
        final int n = g.n;
        initialAdoptedTime = 0;
        finalAdoptedTime = 0;
        Scount = n;
        Acount = 0;
        Rcount = 0;
        PhiCount = 0;
        aSize = 0;
        activeDirCount = 0;
        activeUndirCount = 0;

        Arrays.fill(status, Status.S);
        Arrays.fill(infectedCount, 0);
        Arrays.fill(tInfect, Double.NaN);
        Arrays.fill(tRecover, Double.NaN);
        Arrays.fill(aPos, -1);
        Arrays.fill(dirPos, -1);
        Arrays.fill(undirPos, -1);

        for (int u = 0; u < n; u++) {
            if (thresholdList[u] == 1) PhiCount++;
        }

        // 初期感染者の投入
        boolean[] seen = new boolean[n];
        for (int u : initialInfecteds) {
            if (u < 0 || u >= n) throw new IllegalArgumentException("Invalid initial infected: " + u);
            if (seen[u]) continue;
            seen[u] = true;

            Scount--;
            Acount++;
            status[u] = Status.A;
            tInfect[u] = 0.0;
            addANode(u);
            addOutEdges(u);
        }
        record(0.0);

        // Gillespie メインループ
        double t = 0.0;
        while (true) {
            double rateRecovery = aSize * mu;
            double rateDir = activeDirCount * lambdaDirected;
            double rateUndir = activeUndirCount * lambdaNondirected;
            double totalRate = rateRecovery + rateDir + rateUndir;

            if (totalRate <= 0.0) break;

            double dt = -Math.log(1.0 - random.nextDouble()) / totalRate;
            t += dt;
            if (t >= tMax) break;

            // イベント選択
            double r = random.nextDouble() * totalRate;
            if (r < rateRecovery) {
                // 回復イベント
                int idx = random.nextInt(aSize);
                int u = aList[idx];
                processRecover(u, t);
            } else if (r < rateRecovery + rateDir) {
                // 有向辺での伝播
                int idx = random.nextInt(activeDirCount);
                int edgeIdx = activeDirEdges[idx];
                int target = g.getOutNeighbor(edgeIdx);
                removeDirEdge(edgeIdx);
                processTransmit(target, t);
            } else {
                // 無向辺での伝播
                int idx = random.nextInt(activeUndirCount);
                int edgeIdx = activeUndirEdges[idx];
                int target = g.getOutNeighbor(edgeIdx);
                removeUndirEdge(edgeIdx);
                processTransmit(target, t);
            }
        }

        return new SARResult(n, times, S, A, R, Phi, initialAdoptedTime, finalAdoptedTime, tInfect, tRecover);
    }

    // ── イベント処理 ──

    private void processTransmit(int target, double t) {
        infectedCount[target]++;
        if (infectedCount[target] == thresholdList[target]) {
            PhiCount--;
            if (status[target] == Status.S) {
                Scount--;
                Acount++;
                record(t);

                if (initialAdoptedTime == 0) initialAdoptedTime = t;
                finalAdoptedTime = t;

                status[target] = Status.A;
                tInfect[target] = t;
                addANode(target);
                addOutEdges(target);
            }
        } else if (infectedCount[target] == thresholdList[target] - 1) {
            PhiCount++;
        }
    }

    private void processRecover(int u, double t) {
        Acount--;
        Rcount++;
        record(t);

        status[u] = Status.R;
        tRecover[u] = t;

        removeANode(u);
        removeOutEdges(u);
    }

    // ── Indexed list 操作（O(1) swap-removal）──

    private void addANode(int u) {
        aPos[u] = aSize;
        aList[aSize] = u;
        aSize++;
    }

    private void removeANode(int u) {
        int pos = aPos[u];
        aSize--;
        int last = aList[aSize];
        aList[pos] = last;
        aPos[last] = pos;
        aPos[u] = -1;
    }

    private void addOutEdges(int u) {
        IntRange range = g.outNeighborRange(u);
        for (int e = range.start; e < range.end; e++) {
            if (g.isOutUndirected(e)) {
                undirPos[e] = activeUndirCount;
                activeUndirEdges[activeUndirCount] = e;
                activeUndirCount++;
            } else {
                dirPos[e] = activeDirCount;
                activeDirEdges[activeDirCount] = e;
                activeDirCount++;
            }
        }
    }

    private void removeOutEdges(int u) {
        IntRange range = g.outNeighborRange(u);
        for (int e = range.start; e < range.end; e++) {
            if (g.isOutUndirected(e)) {
                if (undirPos[e] >= 0) removeUndirEdge(e);
            } else {
                if (dirPos[e] >= 0) removeDirEdge(e);
            }
        }
    }

    private void removeDirEdge(int edgeIdx) {
        int pos = dirPos[edgeIdx];
        activeDirCount--;
        int lastEdge = activeDirEdges[activeDirCount];
        activeDirEdges[pos] = lastEdge;
        dirPos[lastEdge] = pos;
        dirPos[edgeIdx] = -1;
    }

    private void removeUndirEdge(int edgeIdx) {
        int pos = undirPos[edgeIdx];
        activeUndirCount--;
        int lastEdge = activeUndirEdges[activeUndirCount];
        activeUndirEdges[pos] = lastEdge;
        undirPos[lastEdge] = pos;
        undirPos[edgeIdx] = -1;
    }

    private void record(double t) {
        times.add(t);
        S.add(Scount);
        A.add(Acount);
        R.add(Rcount);
        Phi.add(PhiCount);
    }

    /**
     * SAR シミュレーション（Gillespie 方式）を実行する。
     *
     * @param g                   グラフ
     * @param lambdaDirected      有向辺の感染率
     * @param lambdaNondirected   無向辺の感染率
     * @param mu                  回復率
     * @param tMax                シミュレーション終了時刻
     * @param thresholdList        各ノードの閾値リスト
     * @param initialInfecteds    初期感染者リスト
     * @param seed                乱数シード
     * @return シミュレーション結果
     */
    public static SARResult simulate(DirectedGraph g, double lambdaDirected, double lambdaNondirected, double mu,
            double tMax, int[] thresholdList, int[] initialInfecteds, long seed) {
        return new SARGillespieSimulator(g, lambdaDirected, lambdaNondirected, mu, tMax, thresholdList, seed)
                .run(initialInfecteds);
    }
}
