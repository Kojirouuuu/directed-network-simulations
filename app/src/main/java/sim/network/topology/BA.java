package sim.network.topology;

import sim.network.DirectedGraph;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Barabási–Albert 型のスケールフリーグラフ生成。
 * <ul>
 *   <li>Directed: 新規ノードは既存ノードを <strong>入次数</strong>に比例する確率で選び、そのノードへ有向辺を張る。</li>
 *   <li>Undirected: 新規ノードは既存ノードを <strong>無向次数</strong>に比例する確率で選び、辺を張る。</li>
 * </ul>
 */
public class BA {
    private BA() {}

    /**
     * @param name グラフ名
     * @param N 頂点数
     * @param m0 初期完全グラフの頂点数
     * @param m 各新規ノードが接続する辺（弧）の数
     * @param isDirected true で有向（入次数に比例）、false で無向（次数に比例）
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(String name, int N, int m0, int m, boolean isDirected, long seed) {
        if (name == null) throw new IllegalArgumentException("name must be non-null");
        if (N <= 0) throw new IllegalArgumentException("Number of nodes N must be positive");
        if (m0 <= 0 || m0 > N) throw new IllegalArgumentException("m0 must be in (0, N]");
        if (m < 0 || m > m0) throw new IllegalArgumentException("m must be in [0, m0]");

        // 有向: 弧数 = 初期 m0(m0-1) + 新規 m(N-m0). 無向: 辺数 = 初期 m0(m0-1)/2 + 新規 m(N-m0)（逆向きはAPIが追加）
        int numEdges = isDirected ? m0 * (m0 - 1) + m * (N - m0) : m0 * (m0 - 1) / 2 + m * (N - m0);
        final int[] src = new int[numEdges];
        final int[] dst = new int[numEdges];
        final boolean[] isUndirected = new boolean[numEdges];
        if (!isDirected) {
            for (int i = 0; i < numEdges; i++) isUndirected[i] = true;
        }

        // 優先的選択用リスト長: 有向＝弧数（終点のみ）、無向＝2*辺数（両端点）
        int prefCapacity = numEdges * 2;
        int[] pref = new int[prefCapacity];
        int prefLen = 0;

        Random rng = new Random(seed);
        int curEdge = 0;

        // 初期完全グラフ
        for (int i = 0; i < m0; i++) {
            for (int j = i + 1; j < m0; j++) {
                src[curEdge] = i;
                dst[curEdge] = j;
                curEdge++;
                pref[prefLen++] = j;  // 入次数: 終点
                pref[prefLen++] = i;

                if (isDirected) {
                    src[curEdge] = j;
                    dst[curEdge] = i;
                    curEdge++;
                }
            }
        }

        // 新規ノードの付加
        for (int u = m0; u < N; u++) {
            Set<Integer> targets = new HashSet<>();
            while (targets.size() < m) {
                int v = pref[rng.nextInt(prefLen)];
                targets.add(v);
            }
            for (int v : targets) {
                src[curEdge] = u;
                dst[curEdge] = v;
                curEdge++;
                pref[prefLen++] = u;
                pref[prefLen++] = v;
            }
        }

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, N, src, dst, isUndirected);
    }
}
