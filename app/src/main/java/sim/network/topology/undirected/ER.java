package sim.network.topology.undirected;

import sim.network.DirectedGraph;

import java.util.*;

/**
 * ERモデル（Erdős–Rényi型ランダムグラフ）を生成する。
 * 全辺を無向として扱い、有効辺（有向のみの辺）は持たない。
 */
public class ER {

    /**
     * 確率 p で各ペアに辺を張る ER グラフを生成する。
     *
     * @param N ノード数
     * @param p エッジ生成確率（0.0〜1.0）
     * @param seed 乱数シード
     * @return 全辺無向の DirectedGraph（有効辺なし）
     */
    public static DirectedGraph generateERFromP(int N, double p, long seed) {
        if (N <= 0) {
            throw new IllegalArgumentException("ノード数Nは正の整数である必要があります");
        }
        if (p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("確率pは0.0〜1.0の範囲で指定してください");
        }

        Random random = new Random(seed);

        List<Integer> src = new ArrayList<>();
        List<Integer> dst = new ArrayList<>();

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (random.nextDouble() < p) {
                    src.add(i);
                    dst.add(j);
                }
            }
        }

        int m = src.size();
        int[] s = new int[m];
        int[] d = new int[m];
        boolean[] undirected = new boolean[m];
        Arrays.fill(undirected, true);

        for (int k = 0; k < m; k++) {
            s[k] = src.get(k);
            d[k] = dst.get(k);
        }
        return DirectedGraph.fromEdgeListWithUndirectedFlag("ER", N, s, d, undirected);
    }

    /**
     * シード省略版
     */
    public static DirectedGraph generateERFromP(int N, double p) {
        return generateERFromP(N, p, System.currentTimeMillis());
    }

    /**
     * 平均次数 kAve になるように、無作為に辺を m = floor(N * kAve / 2) 本選ぶ ER グラフを生成する。
     * 重複辺・自己ループは作らない。
     *
     * @param N ノード数
     * @param kAve 目標とする平均次数（辺数 m は N*kAve/2 に近づく）
     * @param seed 乱数シード
     * @return 全辺無向の DirectedGraph（有効辺なし）
     */
    public static DirectedGraph generateERFromKAve(int N, double kAve, long seed) {
        if (N <= 0) {
            throw new IllegalArgumentException("ノード数Nは正の整数である必要があります");
        }
        if (kAve < 0.0) {
            throw new IllegalArgumentException("平均次数kAveは0以上で指定してください");
        }

        Random random = new Random(seed);
        int maxEdges = (int) Math.floor(N * kAve);
        int m = maxEdges / 2;

        if (m <= 0) {
            int[] emptyS = new int[0];
            int[] emptyD = new int[0];
            boolean[] emptyU = new boolean[0];
            return DirectedGraph.fromEdgeListWithUndirectedFlag("ER", N, emptyS, emptyD, emptyU);
        }

        // 可能な無向辺の最大数
        long maxPossible = (long) N * (N - 1) / 2;
        if (m > maxPossible) {
            m = (int) maxPossible;
        }

        Set<String> selectedEdges = new HashSet<>();
        int[] s = new int[m];
        int[] d = new int[m];
        int edgeCount = 0;

        while (edgeCount < m) {
            int u = random.nextInt(N);
            int v = random.nextInt(N);
            if (u == v) {
                continue;
            }

            String edgeKey = (u < v) ? u + "-" + v : v + "-" + u;
            if (!selectedEdges.contains(edgeKey)) {
                selectedEdges.add(edgeKey);
                s[edgeCount] = u;
                d[edgeCount] = v;
                edgeCount++;
            }
        }

        boolean[] undirected = new boolean[m];
        Arrays.fill(undirected, true);
        return DirectedGraph.fromEdgeListWithUndirectedFlag("ER", N, s, d, undirected);
    }

    /**
     * シード省略版
     */
    public static DirectedGraph generateERFromKAve(int N, double kAve) {
        return generateERFromKAve(N, kAve, System.currentTimeMillis());
    }
}
