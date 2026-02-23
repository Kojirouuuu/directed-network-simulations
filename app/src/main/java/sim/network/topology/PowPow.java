package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;

import java.util.Arrays;
import java.util.Random;

/**
 * べき則次数分布を持つ有向ネットワークを生成するクラス。
 * swapNum により、同一ノードの入次数・出次数の相関を制御できる。
 *
 * swapNum > 0 : 正相関（swapNum 回スワップ、高出次数ノードが高入次数も持ちやすくなる）
 * swapNum = 0 : 無相関（スワップなし）
 * swapNum < 0 : 負相関（|swapNum| 回スワップ、高出次数ノードが低入次数を持ちやすくなる）
 */
public final class PowPow {

    private static final long SEED_MIX_DIRECTED = 0x9E3779B97F4A7C15L;
    private static final long SEED_MIX_UNDIRECTED = 0x517CC1B727220A95L;
    private static final long SEED_MIX_BALANCE = 0xDEADBEEFCAFEBABEL;
    private static final long SEED_MIX_SWAP = 0x1234567890ABCDEFL;

    private PowPow() {
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
     * <p>
     * アルゴリズム概要:
     * <ol>
     * <li>べき則分布から ko（出次数列）と ki（入次数列）を独立にサンプリングする。</li>
     * <li>Σko = Σki になるよう ki を微調整する（balanceSums）。</li>
     * <li>correlateBySwap で 2 頂点を繰り返し選び、swapNum の符号に応じて
     * ki を並び替えることで ko と ki の相関を注入する（swapNum 回または |swapNum| 回スワップ）。</li>
     * <li>コンフィギュレーションモデルで辺を張る。</li>
     * </ol>
     *
     * @param name グラフ名
     * @param n 頂点数
     * @param kdMin 最小次数 (≥ 1)
     * @param kdMax 最大次数 (≥ kdMin)
     * @param gamma べき則指数 (> 1.0)
     * @param swapNum スワップ回数。正: 正相関（swapNum 回）、負: 負相関（|swapNum| 回）、0: 無相関
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(
            String name, int n, int kdMin, int kdMax,
            double gamma, int swapNum, long seed) {

        validateGenerateParams(name, n, kdMin, kdMax, gamma, swapNum);

        // ① べき則からそれぞれ独立にサンプリング
        int[] ko = samplePowerLawDistribution(n, kdMin, kdMax, gamma, seed);
        int[] ki = samplePowerLawDistribution(n, kdMin, kdMax, gamma, seed ^ SEED_MIX_DIRECTED);

        // ② Σko = Σki を保証する（コンフィギュレーションモデルの必要条件）
        balanceSums(ko, ki, kdMin, kdMax, seed ^ SEED_MIX_BALANCE);

        // ③ スワップで同一ノードの ko と ki に相関を注入する
        if (swapNum != 0) {
            correlateBySwap(ko, ki, swapNum, seed ^ SEED_MIX_SWAP);
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
            String name, int n, int kdMin, int kdMax, double gamma, int swapNum) {
        if (name == null)
            throw new IllegalArgumentException("name must be non-null");
        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");
        if (kdMin < 1)
            throw new IllegalArgumentException("kdMin must be >= 1");
        if (kdMax < kdMin)
            throw new IllegalArgumentException("kdMax must be >= kdMin");
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
            int n, int kdMin, int kdMax, double gamma, long seed) {
        Random rng = new Random(seed);
        double[] cdf = buildPowerLawCdf(kdMin, kdMax, gamma);
        int[] k = new int[n];
        for (int u = 0; u < n; u++) {
            double r = rng.nextDouble();
            int idx = Arrays.binarySearch(cdf, r);
            // binarySearch が負の場合: -(insertion point) - 1 → insertion point = -idx - 1
            k[u] = (idx >= 0 ? idx : -idx - 1) + kdMin;
        }
        return k;
    }

    private static double[] buildPowerLawCdf(int kdMin, int kdMax, double gamma) {
        int len = kdMax - kdMin + 1;
        double[] cdf = new double[len];
        double s = 0.0;
        for (int k = kdMin; k <= kdMax; k++) {
            s += Math.pow(k, -gamma);
            cdf[k - kdMin] = s;
        }
        for (int i = 0; i < len; i++) {
            cdf[i] /= s;
        }
        cdf[len - 1] = 1.0; // 浮動小数点誤差対策
        return cdf;
    }

    /**
     * Σko = Σki となるよう ki を ±1 ずつ調整する。
     * 次数分布の形状をできるだけ保つため、1 回の調整は 1 ずつとする。
     */
    private static void balanceSums(int[] ko, int[] ki, int kdMin, int kdMax, long seed) {
        Random rng = new Random(seed);
        int n = ki.length;
        long sumKo = sum(ko);
        long sumKi = sum(ki);

        // ki が不足 → ランダムなノードに +1
        while (sumKi < sumKo) {
            int idx = rng.nextInt(n);
            if (ki[idx] < kdMax) {
                ki[idx]++;
                sumKi++;
            }
        }
        // ki が過剰 → 次数 > 0 のランダムなノードから -1
        while (sumKi > sumKo) {
            int idx = rng.nextInt(n);
            if (ki[idx] > kdMin) {
                ki[idx]--;
                sumKi--;
            }
        }
    }

    /**
     * 先生のアルゴリズム：2 頂点を無作為に選び ki をスワップすることで
     * 同一ノードの ko（出次数）と ki（入次数）の間に相関を作る。
     *
     * <p>
     * 正相関の場合 (swapNum > 0):
     * 
     * <pre>
     *   ko[i] > ko[j]  かつ  ki[i] < ki[j]  → 交換（i に高 ki を渡す）
     *   ko[i] < ko[j]  かつ  ki[i] > ki[j]  → 交換（j に高 ki を渡す）
     * </pre>
     * 
     * つまり「ko の大小関係と ki の大小関係が不一致」のときに交換し、一致に近づける。
     * 繰り返すほど「高 ko ↔ 高 ki」の正相関が強まる。
     *
     * <p>
     * 負相関の場合 (swapNum < 0):
     * 「ko の大小関係と ki の大小関係が一致」しているときに交換し、不一致に近づける。
     *
     * <p>
     * 試行回数 = |swapNum|。swapNum の絶対値がスワップ回数となる。
     *
     * @param ko 各ノードの出次数配列（変更しない）
     * @param ki 各ノードの入次数配列（並び替える）
     * @param swapNum スワップ回数（正: 正相関、負: 負相関、0 は呼び出し側でスキップ）
     * @param seed 乱数シード
     */
    private static void correlateBySwap(int[] ko, int[] ki, int swapNum, long seed) {
        int n = ko.length;
        int nTrials = Math.abs(swapNum);
        boolean positive = (swapNum > 0);
        Random rng = new Random(seed);

        for (int t = 0; t < nTrials; t++) {
            int i = rng.nextInt(n);
            int j = rng.nextInt(n);
            if (i == j)
                continue;

            boolean koiGtkoj = ko[i] > ko[j]; // i の方が出次数が高い
            boolean kiiGtkij = ki[i] > ki[j]; // i の方が入次数が高い

            if (positive) {
                // 正相関：大小関係が不一致 → 交換して一致に近づける
                if (koiGtkoj != kiiGtkij) {
                    int tmp = ki[i];
                    ki[i] = ki[j];
                    ki[j] = tmp;
                }
            } else {
                // 負相関：大小関係が一致 かつ 値が異なる → 交換して不一致に近づける
                if (koiGtkoj == kiiGtkij && ki[i] != ki[j]) {
                    int tmp = ki[i];
                    ki[i] = ki[j];
                    ki[j] = tmp;
                }
            }
        }
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
