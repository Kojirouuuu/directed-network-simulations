package sim.simulation;

import sim.network.DirectedGraph;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.Random;

/**
 * ボンドパーコレーションシミュレータ。
 * 各辺を確率 p で残し、生き残った辺から成るサブグラフの
 * GWCC（最大弱連結成分）と GSCC（最大強連結成分）を計測する。
 */
public final class BondPercolationSimulator {
    private final DirectedGraph g;

    public BondPercolationSimulator(DirectedGraph g) {
        this.g = g;
    }

    /**
     * パーコレーションを実行して結果を返す。
     *
     * @param p   辺の生存確率
     * @param rng 乱数生成器
     * @return パーコレーション結果
     */
    @SuppressWarnings("unchecked")
    public BondPercolationResult run(double p, Random rng) {
        int n = g.n;

        // Step 1: 辺サバイバルと adjOut / adjIn の構築
        List<Integer>[] adjOut = new ArrayList[n];
        List<Integer>[] adjIn = new ArrayList[n];
        for (int i = 0; i < n; i++) {
            adjOut[i] = new ArrayList<>();
            adjIn[i] = new ArrayList<>();
        }

        for (int u = 0; u < n; u++) {
            DirectedGraph.IntRange range = g.outNeighborRange(u);
            for (int e = range.start; e < range.end; e++) {
                int v = g.getOutNeighbor(e);
                if (u == v) continue; // 自己ループ除外

                if (g.isOutUndirected(e)) {
                    // 無向辺は u < v のガードでコイントスを1回だけ
                    if (u < v) {
                        if (rng.nextDouble() < p) {
                            adjOut[u].add(v);
                            adjOut[v].add(u);
                            adjIn[v].add(u);
                            adjIn[u].add(v);
                        }
                    }
                } else {
                    // 有向辺は独立にコイントス
                    if (rng.nextDouble() < p) {
                        adjOut[u].add(v);
                        adjIn[v].add(u);
                    }
                }
            }
        }

        // Step 2: GWCC — BFS（adjOut ∪ adjIn）
        int[] wccComp = new int[n];
        Arrays.fill(wccComp, -1);
        int wccCount = 0;
        int[] wccSizes = new int[n];
        ArrayDeque<Integer> queue = new ArrayDeque<>();

        for (int start = 0; start < n; start++) {
            if (wccComp[start] != -1) continue;
            int compId = wccCount++;
            wccComp[start] = compId;
            int size = 1;
            queue.clear();
            queue.add(start);
            while (!queue.isEmpty()) {
                int u = queue.poll();
                for (int v : adjOut[u]) {
                    if (wccComp[v] == -1) {
                        wccComp[v] = compId;
                        size++;
                        queue.add(v);
                    }
                }
                for (int v : adjIn[u]) {
                    if (wccComp[v] == -1) {
                        wccComp[v] = compId;
                        size++;
                        queue.add(v);
                    }
                }
            }
            wccSizes[compId] = size;
        }

        int gwccSize = 0;
        int gwccSecondSize = 0;
        long wccSizeSum = 0;
        for (int i = 0; i < wccCount; i++) {
            int s = wccSizes[i];
            wccSizeSum += s;
            if (s > gwccSize) {
                gwccSecondSize = gwccSize;
                gwccSize = s;
            } else if (s > gwccSecondSize) {
                gwccSecondSize = s;
            }
        }
        double gwccMeanSize = wccCount == 0 ? 0.0 : (double) wccSizeSum / wccCount;

        // Step 3: GSCC — Kosaraju's アルゴリズム（反復 DFS）

        // パス1: adjOut で後退順（finish order）を記録
        int[] state = new int[n]; // 0=未訪問, 1=処理中, 2=完了
        ArrayDeque<Integer> finishOrder = new ArrayDeque<>();
        Deque<int[]> dfsStack = new ArrayDeque<>();

        for (int start = 0; start < n; start++) {
            if (state[start] != 0) continue;
            dfsStack.clear();
            dfsStack.push(new int[]{start, 0});
            state[start] = 1;
            while (!dfsStack.isEmpty()) {
                int[] top = dfsStack.peek();
                int u = top[0];
                int idx = top[1];
                List<Integer> neighbors = adjOut[u];
                if (idx < neighbors.size()) {
                    top[1]++;
                    int v = neighbors.get(idx);
                    if (state[v] == 0) {
                        state[v] = 1;
                        dfsStack.push(new int[]{v, 0});
                    }
                } else {
                    dfsStack.pop();
                    state[u] = 2;
                    finishOrder.addFirst(u);
                }
            }
        }

        // パス2: adjIn（転置グラフ）で finishOrder の降順に BFS
        int[] sccComp = new int[n];
        Arrays.fill(sccComp, -1);
        int sccCount = 0;
        int[] sccSizes = new int[n];

        for (int start : finishOrder) {
            if (sccComp[start] != -1) continue;
            int compId = sccCount++;
            sccComp[start] = compId;
            int size = 1;
            queue.clear();
            queue.add(start);
            while (!queue.isEmpty()) {
                int u = queue.poll();
                for (int v : adjIn[u]) {
                    if (sccComp[v] == -1) {
                        sccComp[v] = compId;
                        size++;
                        queue.add(v);
                    }
                }
            }
            sccSizes[compId] = size;
        }

        int gsccSize = 0;
        int gsccSecondSize = 0;
        long sccSizeSum = 0;
        for (int i = 0; i < sccCount; i++) {
            int s = sccSizes[i];
            sccSizeSum += s;
            if (s > gsccSize) {
                gsccSecondSize = gsccSize;
                gsccSize = s;
            } else if (s > gsccSecondSize) {
                gsccSecondSize = s;
            }
        }
        double gsccMeanSize = sccCount == 0 ? 0.0 : (double) sccSizeSum / sccCount;

        return new BondPercolationResult(n,
                gwccSize, gwccSecondSize, gwccMeanSize,
                gsccSize, gsccSecondSize, gsccMeanSize);
    }

    /**
     * 静的ファクトリメソッド。
     *
     * @param g    グラフ
     * @param p    辺の生存確率
     * @param seed 乱数シード
     * @return パーコレーション結果
     */
    public static BondPercolationResult simulate(DirectedGraph g, double p, long seed) {
        return new BondPercolationSimulator(g).run(p, new Random(seed));
    }
}
