package sim.network.topology.undirected;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;

import java.util.Arrays;
import java.util.Random;

/**
 * 無向コンフィギュレーションモデル。
 * 与えられた次数列（またはべき則分布からサンプリングした次数列）に従い、
 * 無向辺のみを持つネットワークを生成する。
 * 出力は DirectedGraph で、全辺は無向として扱われる（isUndirected=true）。
 */
public final class CM {

    private static final long SEED_MIX_UNDIRECTED = 0x517CC1B727220A95L;

    private CM() {
    }

    // =========================================================
    // Public API
    // =========================================================

    /**
     * 無向次数列から無向コンフィギュレーションモデルを生成する。
     * このクラスは無向専用のため、ki と ko はすべて 0 である必要がある。
     *
     * @param name グラフ名
     * @param ki 各ノードの入次数配列（無向専用のためすべて 0）
     * @param ko 各ノードの出次数配列（無向専用のためすべて 0）
     * @param kn 各ノードの無向次数配列
     * @param seed 乱数シード
     * @return 生成された DirectedGraph（全辺無向）
     */
    public static DirectedGraph generateFromDegreeSequence(
            String name, int[] ki, int[] ko, int[] kn, long seed) {

        validateDegreeSequence(name, ki, ko, kn);

        final int n = ki.length;
        final int mDir = (int) sum(ko);
        final int mUnd = (int) (sum(kn) / 2);
        final int mBase = mDir + mUnd;

        // 無向専用: 有向辺は使わない（mDir == 0 である必要あり）
        int[] outStubs = new int[mDir];
        int[] inStubs = new int[mDir];

        // --- 無向スタブ構築・シャッフル ---
        int[] undirectedStubs = buildUndirectedStubs(n, kn, mUnd);
        ArrayUtils.shuffle(undirectedStubs, seed ^ SEED_MIX_UNDIRECTED);

        // --- 辺配列を生成 ---
        int[] srcs = new int[mBase];
        int[] dsts = new int[mBase];
        boolean[] isUndirected = new boolean[mBase];
        fillBaseEdgeArrays(srcs, dsts, isUndirected, outStubs, inStubs, undirectedStubs, mDir, mUnd);

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, n, srcs, dsts, isUndirected);
    }

    /**
     * べき則次数分布 P(k) ∝ k^{-gamma} に従う無向コンフィギュレーションモデルを生成する。
     *
     * <p>
     * アルゴリズム概要:
     * <ol>
     * <li>べき則分布から各ノードの次数 ku をサンプリングする。</li>
     * <li>Σku が偶数になるよう 1 ノードの次数を ±1 調整する。</li>
     * <li>無向スタブをランダムにシャッフルし、隣接するペアを辺として張る。</li>
     * </ol>
     *
     * @param name グラフ名
     * @param n 頂点数
     * @param kuMin 最小次数 (≥ 1)
     * @param kuMax 最大次数 (≥ kuMin)
     * @param gamma べき則指数 (> 1.0)
     * @param seed 乱数シード
     * @return 生成された DirectedGraph（全辺無向）
     */
    public static DirectedGraph generate(
            String name, int n, int kuMin, int kuMax,
            double gamma, long seed) {

        validateGenerateParams(name, n, kuMin, kuMax, gamma);

        // ① べき則からそれぞれ独立にサンプリング
        int[] ku = samplePowerLawDistribution(n, kuMin, kuMax, gamma, seed);

        // ② Σkuが偶数になるように調整する（コンフィギュレーションモデルの必要条件）
        ensureEvenNondirectedSum(ku, seed);

        // ③ 無向辺のみ（ki, ko は 0、次数は ku のみ）
        int[] ki = new int[n];
        int[] ko = new int[n];

        return generateFromDegreeSequence(name, ki, ko, ku, seed ^ SEED_MIX_UNDIRECTED);
    }

    // =========================================================
    // Validation
    // =========================================================

    private static void validateDegreeSequence(String name, int[] ki, int[] ko, int[] kn) {
        if (name == null)
            throw new IllegalArgumentException("name must be non-null");
        if (ki == null || ko == null || kn == null)
            throw new IllegalArgumentException("Degree arrays must be non-null");
        if (ki.length != ko.length || ki.length != kn.length)
            throw new IllegalArgumentException("Degree arrays must have the same length");

        int n = ki.length;
        long sumKi = 0, sumKo = 0, sumKn = 0;
        for (int u = 0; u < n; u++) {
            if (ki[u] < 0 || ko[u] < 0 || kn[u] < 0)
                throw new IllegalArgumentException("Degree values must be non-negative");
            sumKi += ki[u];
            sumKo += ko[u];
            sumKn += kn[u];
        }
        if (sumKi != sumKo)
            throw new IllegalArgumentException(
                    "Sum of in-degrees must equal sum of out-degrees: sumKi=" + sumKi + ", sumKo=" + sumKo);
        if (sumKi != 0 || sumKo != 0)
            throw new IllegalArgumentException(
                    "Undirected CM only: ki and ko must be all zeros");
        if ((sumKn & 1L) != 0L)
            throw new IllegalArgumentException("Sum of nondirected degrees must be even");
        if (sumKo > Integer.MAX_VALUE || sumKn > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Sum of degrees must be <= Integer.MAX_VALUE");
    }

    private static void validateGenerateParams(
            String name, int n, int kuMin, int kuMax, double gamma) {
        if (name == null)
            throw new IllegalArgumentException("name must be non-null");
        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");
        if (kuMin < 1)
            throw new IllegalArgumentException("kuMin must be >= 1");
        if (kuMax < kuMin)
            throw new IllegalArgumentException("kuMax must be >= kuMin");
        if (gamma <= 1.0)
            throw new IllegalArgumentException("gamma must be > 1.0");
    }

    // =========================================================
    // Degree sequence helpers
    // =========================================================

    private static long sum(int[] arr) {
        long s = 0;
        for (int x : arr)
            s += x;
        return s;
    }

    /**
     * 無向次数の合計を偶数にする（ペアリングのため）。ku を破壊的に更新する。
     */
    private static void ensureEvenNondirectedSum(int[] ku, long seed) {
        long sumKu = sum(ku);
        if ((sumKu & 1L) == 0L) return;
        Random rng = new Random(seed);
        int n = ku.length;
        for (int attempt = 0; attempt < n * 2; attempt++) {
            int u = rng.nextInt(n);
            if (ku[u] > 0) {
                ku[u]--;
                return;
            }
        }
        ku[rng.nextInt(n)]++;
    }

    /**
     * べき則分布 P(k) ∝ k^{-gamma} に従う次数列をサンプリングする。
     */
    private static int[] samplePowerLawDistribution(
            int n, int kuMin, int kuMax, double gamma, long seed) {
        Random rng = new Random(seed);
        double[] cdf = buildPowerLawCdf(kuMin, kuMax, gamma);
        int[] k = new int[n];
        for (int u = 0; u < n; u++) {
            double r = rng.nextDouble();
            int idx = Arrays.binarySearch(cdf, r);
            // binarySearch が負の場合: -(insertion point) - 1 → insertion point = -idx - 1
            k[u] = (idx >= 0 ? idx : -idx - 1) + kuMin;
        }
        return k;
    }

    private static double[] buildPowerLawCdf(int kuMin, int kuMax, double gamma) {
        int len = kuMax - kuMin + 1;
        double[] cdf = new double[len];
        double s = 0.0;
        for (int k = kuMin; k <= kuMax; k++) {
            s += Math.pow(k, -gamma);
            cdf[k - kuMin] = s;
        }
        for (int i = 0; i < len; i++) {
            cdf[i] /= s;
        }
        cdf[len - 1] = 1.0; // 浮動小数点誤差対策
        return cdf;
    }

    // =========================================================
    // Stub and edge construction
    // =========================================================

    /**
     * 無向スタブ配列を構築する。
     * undStubs[2j] と undStubs[2j+1] が j 番目の辺の両端点スタブ。
     */
    private static int[] buildUndirectedStubs(int n, int[] kn, int mUnd) {
        int[] undStubs = new int[mUnd * 2];
        int p = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < kn[u]; t++) {
                undStubs[p++] = u;
            }
        }
        if (p != mUnd * 2) {
            throw new IllegalStateException(
                    "Undirected stub count mismatch: p=" + p + ", 2*mUnd=" + (mUnd * 2));
        }
        return undStubs;
    }

    /**
     * 有向・無向スタブから辺配列 (srcs, dsts, isUndirected) を埋める。
     */
    private static void fillBaseEdgeArrays(
            int[] srcs, int[] dsts, boolean[] isUndirected,
            int[] outStubs, int[] inStubs, int[] undStubs,
            int mDir, int mUnd) {
        int e = 0;
        for (int i = 0; i < mDir; i++) {
            srcs[e] = outStubs[i];
            dsts[e] = inStubs[i];
            isUndirected[e] = false;
            e++;
        }
        for (int j = 0; j < mUnd; j++) {
            srcs[e] = undStubs[j * 2];
            dsts[e] = undStubs[j * 2 + 1];
            isUndirected[e] = true;
            e++;
        }
        if (e != srcs.length) {
            throw new IllegalStateException("Edge count mismatch: e=" + e + ", mBase=" + srcs.length);
        }
    }
}
