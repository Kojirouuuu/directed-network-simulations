package sim.simulation;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.IntRange;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Random;

/**
 * SAR（Susceptible-Adopted-Recovered）シミュレーションを実行するクラス。
 */
public final class SARSimulator {
    /** ノードの状態 */
    public enum Status {
        S, // Susceptible（感受性）
        A, // Adopted（採用）
        R  // Recovered（回復）
    }

    /** イベントの種類 */
    public enum EventType {
        TRANSMIT, // 感染伝播
        RECOVER   // 回復
    }

    private static final class Event {
        final double time; // イベント発生時刻
        final int node; // 対象ノード
        final EventType type; // イベントの種類
        final long seq; // シーケンス番号（同時刻のイベントの順序付け用）

        Event(double time, int node, EventType type, long seq) {
            this.time = time;
            this.node = node;
            this.type = type;
            this.seq = seq;
        }
    }

    private final DirectedGraph g; // グラフ
    private final double lambdaDirected; // 有向辺の感染率
    private final double lambdaNondirected; // 無向辺の感染率
    private final double mu; // 回復率
    private final double tMax; // シミュレーション終了時刻
    private final int[] thresholdList; // 各ノードの閾値リスト
    private final Random random; // 乱数ジェネレータ

    private final Status[] status; // 各ノードの状態
    private final double[] predInfTime; // 各ノードの予測感染時刻
    private final double[] recTime; // 各ノードの回復時刻
    private final double[] tInfect; // 各ノードの感染成立時刻
    private final double[] tRecover; // 各ノードの回復成立時刻

    private final ArrayList<Double> times = new ArrayList<>(); // 記録時刻のリスト
    private final int[] infectedCount; // 各ノードの感染カウント
    private final ArrayList<Integer> S = new ArrayList<>(); // 各時刻の Susceptible 数
    private final ArrayList<Integer> A = new ArrayList<>(); // 各時刻の Adopted 数
    private final ArrayList<Integer> R = new ArrayList<>(); // 各時刻の Recovered 数
    private final ArrayList<Integer> Phi = new ArrayList<>(); // 各時刻のあと一回の伝播で採用されるノード数

    private int Scount; // 現在の Susceptible 数
    private int Acount; // 現在の Adopted 数
    private int Rcount; // 現在の Recovered 数
    private int PhiCount; // 現在のあと一回の伝播で採用されるノード数
    private double initialAdoptedTime; // 初期採用時刻
    private double finalAdoptedTime; // 最終採用時刻
    
    /**
     * SARSimulator を構築する。
     *
     * @param g                   グラフ
     * @param lambdaDirected      有向辺の感染率
     * @param lambdaNondirected   無向辺の感染率
     * @param mu                  回復率
     * @param tMax                シミュレーション終了時刻
     * @param thresholdList        各ノードの閾値リスト
     * @param seed                乱数シード
     */
    public SARSimulator(DirectedGraph g, double lambdaDirected, double lambdaNondirected, double mu, double tMax,
            int[] thresholdList, long seed) {
        if (g == null) throw new IllegalArgumentException("Directed graph must be non-null");
        if (lambdaDirected < 0 || lambdaNondirected < 0) throw new IllegalArgumentException("lambdaDirected and lambdaNondirected must be non-negative");
        if (tMax <= 0) throw new IllegalArgumentException("tMax must be positive");
        if (thresholdList == null || thresholdList.length != g.n) throw new IllegalArgumentException("thresholdList must be an array of length n");

        this.g = g;
        this.lambdaDirected = lambdaDirected;
        this.lambdaNondirected = lambdaNondirected;
        this.mu = mu;
        this.tMax = tMax;
        this.thresholdList = thresholdList;
        this.random = new Random(seed);

        int n = g.n;
        this.status = new Status[n];
        this.predInfTime = new double[n];
        this.recTime = new double[n];
        this.tInfect = new double[n];
        this.tRecover = new double[n];

        Arrays.fill(tInfect, Double.NaN);
        Arrays.fill(tRecover, Double.NaN);
        this.infectedCount = new int[n];
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

        Arrays.fill(infectedCount, 0);

        for (int u = 0; u < n; u++) {
            status[u] = Status.S;
            predInfTime[u] = Double.POSITIVE_INFINITY;
            recTime[u] = Double.POSITIVE_INFINITY;
        }

        final Comparator<Event> cmp = (a, b) -> {
            int c = Double.compare(a.time, b.time);
            if (c != 0) return c;
            c = a.type.compareTo(b.type);
            if (c != 0) return c;
            return Long.compare(a.seq, b.seq);
        };

        final PriorityQueue<Event> Q = new PriorityQueue<>(cmp);
        final SeqGen seqGen = new SeqGen() {
            private long seq = 0L;

            @Override
            public long next() {
                return seq++;
            }
        };

        // 初期感染者の投入
        boolean[] seen = new boolean[n];
        for (int u : initialInfecteds) {
            if (u < 0 || u >= n) {
                throw new IllegalArgumentException("Invalid initial infected: " + u);
            }
            if (seen[u]) {
                continue;
            }
            seen[u] = true;
            // 初期感染者は threshold に達しているとみなす
            infectedCount[u] = thresholdList[u];
            // 初期感染者を直接感染者（A）に設定
            Scount--;
            Acount++;
            status[u] = Status.A;
            tInfect[u] = 0.0;

            // 回復イベントをスケジュール
            double tRec = exp(mu);
            recTime[u] = tRec;
            if (tRec < tMax) {
                Q.add(new Event(tRec, u, EventType.RECOVER, seqGen.next()));
            }

            // 隣接ノードへの感染伝播を開始
            IntRange outNeighborRange = g.outNeighborRange(u);
            for (int e = outNeighborRange.start; e < outNeighborRange.end; e++) {
                int v = g.getOutNeighbor(e);
                if (g.isOutUndirected(e)) {
                    findTransmit(Q, 0.0, u, v, seqGen, lambdaNondirected);
                } else {
                    findTransmit(Q, 0.0, u, v, seqGen, lambdaDirected);
                }
            }
        }
        record(0.0);

        while (!Q.isEmpty()) {
            Event ev = Q.poll();
            final int u = ev.node;
            final double t = ev.time;

            if (t >= tMax) break;

            if (ev.type == EventType.TRANSMIT) {
                processTransmit(u, t, Q, () -> seqGen.next());
            } else {
                // EventType.RECOVER
                if (status[u] == Status.A && t == recTime[u]) {
                    processRecover(u, t);
                }
            }
        }

        return new SARResult(n, times, S, A, R, Phi, initialAdoptedTime, finalAdoptedTime, tInfect, tRecover);
    }

    /**
     * 感染伝播イベントを処理する。
     *
     * @param u      対象ノード
     * @param t      現在時刻
     * @param Q      イベントキュー
     * @param seqGen シーケンス番号ジェネレータ
     */
    private void processTransmit(int u, double t, PriorityQueue<Event> Q, SeqGen seqGen) {
        infectedCount[u]++;
        if (infectedCount[u] >= thresholdList[u]) {
            if (status[u] == Status.S) {
                Scount--;
                Acount++;
                PhiCount--;
                record(t);

                if (initialAdoptedTime == 0) {
                    initialAdoptedTime = t;
                }
                finalAdoptedTime = t;

                status[u] = Status.A;
                tInfect[u] = t;

                double tRec = t + exp(mu);
                recTime[u] = tRec;
                if (tRec < tMax) {
                    Q.add(new Event(tRec, u, EventType.RECOVER, seqGen.next()));
                }

                IntRange outNeighborRange = g.outNeighborRange(u);
                for (int e = outNeighborRange.start; e < outNeighborRange.end; e++) {
                    int v = g.getOutNeighbor(e);
                    if (g.isOutUndirected(e)) {
                        findTransmit(Q, t, u, v, seqGen, lambdaNondirected);
                    } else {
                        findTransmit(Q, t, u, v, seqGen, lambdaDirected);
                    }
                }
            }
        } else if (infectedCount[u] == thresholdList[u] - 1) {
            PhiCount++;
        }
    }

    /**
     * 感染伝播イベントをスケジュールする。
     *
     * @param Q      イベントキュー
     * @param t      現在時刻
     * @param source 感染源ノード
     * @param target 対象ノード
     * @param seqGen シーケンス番号ジェネレータ
     * @param lambda 感染率
     */
    private void findTransmit(PriorityQueue<Event> Q, double t, int source, int target, SeqGen seqGen, double lambda) {
        if (status[target] != Status.S) {
            return;
        }

        double tInf = t + exp(lambda);

        // predInfTime[target] は最初の感染、次点以降の感染でも threshold=1 であれば Status=S で無視イベントとなる
        // double bound = Math.min(recTime[source], Math.min(predInfTime[target], tMax));
        double bound = Math.min(recTime[source], tMax);
        if (tInf < bound) {
            // threshold=1 の場合は、最初の感染のみイベント追加
            if (thresholdList[target] == 1) {
                if (tInf < predInfTime[target]) {
                    predInfTime[target] = tInf;
                    Q.add(new Event(tInf, target, EventType.TRANSMIT, seqGen.next()));
                }
            } else {
                predInfTime[target] = tInf;
                Q.add(new Event(tInf, target, EventType.TRANSMIT, seqGen.next()));
            }
        }
    }

    /**
     * 回復イベントを処理する。
     *
     * @param u 対象ノード
     * @param t 現在時刻
     */
    private void processRecover(int u, double t) {
        Acount--;
        Rcount++;
        record(t);

        status[u] = Status.R;
        tRecover[u] = t;
    }

    /** シーケンス番号を生成するインターフェース */
    private interface SeqGen {
        long next();
    }

    /**
     * 指数分布に従う乱数を生成する。
     *
     * @param rate レート
     * @return 指数分布に従う乱数値
     */
    private double exp(double rate) {
        if (rate <= 0.0) {
            return Double.POSITIVE_INFINITY;
        }

        double u = 1.0 - random.nextDouble();
        return -Math.log(u) / rate;
    }

    /**
     * 現在の状態を記録する。
     *
     * @param t 現在時刻
     */
    private void record(double t) {
        times.add(t);
        S.add(Scount);
        A.add(Acount);
        R.add(Rcount);
        Phi.add(PhiCount);
    }
    
    /**
     * SAR シミュレーションを実行する。
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
        return new SARSimulator(g, lambdaDirected, lambdaNondirected, mu, tMax, thresholdList, seed)
                .run(initialInfecteds);
    }
}
