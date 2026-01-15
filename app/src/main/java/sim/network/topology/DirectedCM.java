package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;

import java.util.Random;

/**
 * 論文 2.2.1 の Directed Networks のネットワーク（Directed + Nondirected の CM）を生成するユーティリティ。
 *
 * 生成モデルの要点:
 * - 各ノードに (k_i, k_o, k_n) を割り当てる。
 * - directed: out-stub と in-stub をランダムにマッチングして u->v を作る。
 * - nondirected: nondirected stub をランダムにペアリングして u--v を作る（DirectedGraph では両向き2アークに展開）。
 *
 * 注意:
 * - Configuration Model と同様に多重辺・自己ループを許容する（論文の前提に整合）。
 * - 「単純グラフ（多重辺なし・自己ループなし）」が必要なら、生成後にリジェクト/リワイヤリングが必要。
 */

public class DirectedCM {
    private DirectedCM() {}

    /**
     * 一般の次数列 (k_i, k_o, k_n) から Directed + Nondirected CM を生成する。
     *
     * @param name グラフ名
     * @param ki   各ノードの in-directed 次数
     * @param ko   各ノードの out-directed 次数
     * @param kn   各ノードの nondirected 次数
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generateFromDegreeSequence(String name, int[] ki, int[] ko, int[] kn, long seed) {
        if (name == null) {
            throw new IllegalArgumentException("name must be non-null");
        }
        if (ki == null || ko == null || kn == null) {
            throw new IllegalArgumentException("Degree arrays must be non-null");
        }
        if (ki.length != ko.length || ki.length != kn.length) {
            throw new IllegalArgumentException("Degree arrays must have the same length");
        }

        final int n = ki.length;

        long sumKi = 0, sumKo = 0, sumKn = 0;
        for (int u = 0; u < n; u++) {
            int a = ki[u], b = ko[u], c = kn[u];
            if (a < 0 || b < 0 || c < 0) throw new IllegalArgumentException("degree values must be non-negative");
            sumKi += a;
            sumKo += b;
            sumKn += c;
        }

        // directed の stub 数は一致していないとマッチングできない
        if (sumKi != sumKo) {
            throw new IllegalArgumentException("Sum of in-degrees must equal sum of out-degrees");
        }

        // nondirected はペアリングするために偶数が必要
        if ((sumKn & 1L) != 0L) {
            throw new IllegalArgumentException("Sum of nondirected degrees must be even");
        }
        if (sumKo > Integer.MAX_VALUE || sumKn > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Sum of degrees must be less than or equal to Integer.MAX_VALUE");
        }

        final int mDir = (int) sumKo; // directed アーク数（base edge としては isUndirected=false のエントリ数）
        final int mUnd = (int) (sumKn / 2); // nondirected 辺数（base edge として isUndirected=true のエントリ数）
        final int mBase = mDir + mUnd;

        int[] outStubs = new int[mDir];
        int[] inStubs = new int[mDir];

        int pOut = 0, pIn = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < ko[u]; t++) {
                outStubs[pOut++] = u;
            }
            for (int t = 0; t < ki[u]; t++) {
                inStubs[pIn++] = u;
            }
        }
        if (pOut != mDir || pIn != mDir) {
            throw new IllegalStateException("Stub count mismatch: pOut=" + pOut + ", pIn=" + pIn + ", mDir=" + mDir);
        }

        ArrayUtils.shuffle(outStubs, seed);
        ArrayUtils.shuffle(inStubs, seed);

        // nondirected: kn-stubs を作ってシャッフルし、(0,1),(2,3),... をペアリング
        int[] undStubs = new int[mUnd * 2];
        int pUnd = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < kn[u]; t++) {
                undStubs[pUnd++] = u;
            }
        }
        if (pUnd != undStubs.length) {
            throw new IllegalStateException("Undirected stub count mismatch: pUnd=" + pUnd + ", expected=" + undStubs.length);
        }

        ArrayUtils.shuffle(undStubs, seed);

        // DirectedGraph に渡す base edge list を作る
        int[] srcs = new int[mBase];
        int[] dsts = new int[mBase];
        boolean[] isUndirected = new boolean[mBase];

        int e = 0;

        // directed base edges（展開なし）
        for (int i = 0; i < mDir; i++) {
            srcs[e] = outStubs[i];
            dsts[e] = inStubs[i];
            isUndirected[e] = false;
            e++;
        }

        // nondirected base edges（DirectedGraph が逆向きを追加）
        for (int j = 0; j < mUnd; j++) {
            int u = undStubs[j * 2];
            int v = undStubs[j * 2 + 1];
            srcs[e] = u;
            dsts[e] = v;
            isUndirected[e] = true;
            e++;
        }

        if (e != mBase) {
            throw new IllegalStateException("Edge count mismatch: e=" + e + ", mBase=" + mBase);
        }

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, n, srcs, dsts, isUndirected);
    }

    /**
     * 指定されたパラメータから Directed + Nondirected CM を生成する。
     *
     * @param name グラフ名
     * @param n    頂点数
     * @param kHat 平均次数
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(String name, int n, int kHat, long seed) {
        if (name == null) throw new IllegalArgumentException("name must be non-null");
        if (n < 0) throw new IllegalArgumentException("n must be non-negative");
        if (kHat < 0) throw new IllegalArgumentException("kHat must be non-negative");

        Random random = new Random(seed);

        long target = (long) n * kHat;
        // kn の総和は（最終的に）n*kHat になるので、偶数でないと無向ペアリングできない
        if ((target & 1L) != 0L) {
            throw new IllegalArgumentException("n * kHat must be even");
        }

        int[] ki = new int[n];
        int[] ko = new int[n];
        int[] kn = new int[n];

        // sum(ko) を target=n*kHat に合わせる（分布の歪みは O(sqrt(n)) 程度で、n が大きいほど無視できる）
        long sumKo = 0;
        for (int u = 0; u < n; u++) {
            ki[u] = kHat;
            ko[u] = random.nextInt(2 * kHat + 1); // [0, 2kHat] の一様乱数
            sumKo += ko[u];
        }

        long delta = sumKo - target;

        // delta>0: ko を減らす必要がある / delta<0: ko を増やす必要がある
        while (delta != 0) {
            int u = random.nextInt(n);
            if (delta > 0) {
                if (ko[u] > 0) {
                    ko[u]--;
                    delta--;
                }
            } else {
                if (ko[u] < 2 * kHat) {
                    ko[u]++;
                    delta++;
                }
            }
        }

        // kn = 2*kHat - ko（常に 0..2kHat）
        for (int u = 0; u < n; u++) {
            kn[u] = 2 * kHat - ko[u];
        }

        return generateFromDegreeSequence(name, ki, ko, kn, seed ^ 0x9E3779B97F4A7C15L);
    }
}
