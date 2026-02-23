package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;

import java.util.Arrays;
import java.util.Random;

/**
 * べき則次数分布を持つ有向ネットワークを生成するクラス。
 */
public final class SameInOut {

    private static final long SEED_MIX_DIRECTED = 0x9E3779B97F4A7C15L;
    private static final long SEED_MIX_UNDIRECTED = 0x517CC1B727220A95L;

    private SameInOut() {
    }

    // =========================================================
    // Public API
    // =========================================================

    /**
     * 一般の次数列 (ki, ko, kn) から有向コンフィギュレーションモデルを生成する。
     *
     * @param name グラフ名
     * @param ki 各ノードの入次数配列
     * @param ko 各ノードの出次数配列
     * @param kn 各ノードの無向次数配列
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generateFromDegreeSequence(
            String name, int[] ki, int[] ko, int[] kn, long seed) {

        validateDegreeSequence(name, ki, ko, kn);

        final int n = ki.length;
        final int mDir = (int) sum(ko); // Σko = Σki が保証済み
        final int mUnd = (int) (sum(kn) / 2);
        final int mBase = mDir + mUnd;

        // --- 有向スタブ構築・シャッフル ---
        int[][] directedStubs = buildDirectedStubs(n, ki, ko, mDir);
        int[] outStubs = directedStubs[0];
        int[] inStubs = directedStubs[1];
        ArrayUtils.shuffle(outStubs, seed);
        ArrayUtils.shuffle(inStubs, seed ^ SEED_MIX_DIRECTED);

        // --- 無向スタブ構築・シャッフル ---
        int[] undStubs = buildUndirectedStubs(n, kn, mUnd);
        ArrayUtils.shuffle(undStubs, seed ^ SEED_MIX_UNDIRECTED);

        // --- 辺配列を生成 ---
        int[] srcs = new int[mBase];
        int[] dsts = new int[mBase];
        boolean[] isUndirected = new boolean[mBase];
        fillBaseEdgeArrays(srcs, dsts, isUndirected, outStubs, inStubs, undStubs, mDir, mUnd);

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, n, srcs, dsts, isUndirected);
    }

    /**
     * べき則次数分布を持つ有向ネットワークを生成する（無向辺なし）。
     *
     * @param name グラフ名
     * @param n 頂点数
     * @param kMin 最小次数 (≥ 1)
     * @param kMax 最大次数 (≥ kMin)
     * @param gamma べき則指数 (> 1.0)
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(
            String name, int n, int kMin, int kMax,
            double gamma, long seed) {

        validateGenerateParams(name, n, kMin, kMax, gamma);

        // ① べき則からそれぞれ独立にサンプリング
        int[] ko = samplePowerLawDistribution(n, kMin, kMax, gamma, seed);
        int stubLength = (int) ko.length;
        int[] ki = new int[stubLength];
        for (int i = 0; i < stubLength; i++) {
            ki[i] = ko[i];
        }

        // ④ 有向辺のみ（kn はすべて 0）
        int[] kn = new int[n];

        return generateFromDegreeSequence(name, ki, ko, kn, seed ^ SEED_MIX_DIRECTED);
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
        if ((sumKn & 1L) != 0L)
            throw new IllegalArgumentException("Sum of nondirected degrees must be even");
        if (sumKo > Integer.MAX_VALUE || sumKn > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Sum of degrees must be <= Integer.MAX_VALUE");
    }

    private static void validateGenerateParams(
            String name, int n, int kMin, int kMax, double gamma) {
        if (name == null)
            throw new IllegalArgumentException("name must be non-null");
        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");
        if (kMin < 1)
            throw new IllegalArgumentException("kMin must be >= 1");
        if (kMax < kMin)
            throw new IllegalArgumentException("kMax must be >= kMin");
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
     * べき則分布 P(k) ∝ k^{-gamma} に従う次数列をサンプリングする。
     */
    private static int[] samplePowerLawDistribution(
            int n, int kMin, int kMax, double gamma, long seed) {
        Random rng = new Random(seed);
        double[] cdf = buildPowerLawCdf(kMin, kMax, gamma);
        int[] k = new int[n];
        for (int u = 0; u < n; u++) {
            double r = rng.nextDouble();
            int idx = Arrays.binarySearch(cdf, r);
            // binarySearch が負の場合: -(insertion point) - 1 → insertion point = -idx - 1
            k[u] = (idx >= 0 ? idx : -idx - 1) + kMin;
        }
        return k;
    }

    private static double[] buildPowerLawCdf(int kMin, int kMax, double gamma) {
        int len = kMax - kMin + 1;
        double[] cdf = new double[len];
        double s = 0.0;
        for (int k = kMin; k <= kMax; k++) {
            s += Math.pow(k, -gamma);
            cdf[k - kMin] = s;
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
     * 有向スタブ配列を構築する。
     *
     * @return { outStubs, inStubs }
     */
    private static int[][] buildDirectedStubs(int n, int[] ki, int[] ko, int mDir) {
        int[] outStubs = new int[mDir];
        int[] inStubs = new int[mDir];
        int pOut = 0, pIn = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < ko[u]; t++)
                outStubs[pOut++] = u;
            for (int t = 0; t < ki[u]; t++)
                inStubs[pIn++] = u;
        }
        if (pOut != mDir || pIn != mDir) {
            throw new IllegalStateException(
                    "Directed stub count mismatch: pOut=" + pOut + ", pIn=" + pIn + ", mDir=" + mDir);
        }
        return new int[][] { outStubs, inStubs };
    }

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
