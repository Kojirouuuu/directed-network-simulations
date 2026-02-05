package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;
import sim.utils.RandomUtils;

import java.util.Random;
import java.util.Arrays;

/**
 * 入次数をパワーロー分布で、出次数をランダムに配分する。
 *
 * 生成モデルの要点:
 * - 各ノードに (k_i, k_o, k_n) を割り当てる。
 * - directed: out-stub と in-stub をランダムにマッチングして u->v を作る。
 * - nondirected: nondirected stub をランダムにペアリングして u--v を作る（DirectedGraph では両向き2アークに展開）。
 *
 * 注意:
 * - 「単純グラフ（多重辺なし・自己ループなし）」が必要なら、生成後にリジェクト/リワイヤリングが必要。
 */
public final class DirectedCMInPow {

    private static final long SEED_MIX_DIRECTED = 0x9E3779B97F4A7C15L;
    private static final long SEED_MIX_UNDIRECTED = 0x517CC1B727220A95L;

    private DirectedCMInPow() {}

    /**
     * 一般の次数列 (k_i, k_o, k_n) から Directed CM を生成する。
     *
     * @param name グラフ名
     * @param ki   各ノードの in-directed 次数
     * @param ko   各ノードの out-directed 次数
     * @param kn   各ノードの nondirected 次数
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generateFromDegreeSequence(String name, int[] ki, int[] ko, int[] kn, long seed) {
        validateDegreeSequence(name, ki, ko, kn);

        final int n = ki.length;
        final int mDir = (int) sum(ko);
        final int mUnd = (int) (sum(kn) / 2);
        final int mBase = mDir + mUnd;

        int[][] directedStubs = buildDirectedStubs(n, ki, ko, mDir);
        int[] outStubs = directedStubs[0];
        int[] inStubs = directedStubs[1];

        ArrayUtils.shuffle(outStubs, seed);
        ArrayUtils.shuffle(inStubs, seed ^ SEED_MIX_DIRECTED);

        int[] undStubs = buildUndirectedStubs(n, kn, mUnd);
        ArrayUtils.shuffle(undStubs, seed ^ SEED_MIX_UNDIRECTED);

        int[] srcs = new int[mBase];
        int[] dsts = new int[mBase];
        boolean[] isUndirected = new boolean[mBase];
        fillBaseEdgeArrays(srcs, dsts, isUndirected, outStubs, inStubs, undStubs, mDir, mUnd);

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, n, srcs, dsts, isUndirected);
    }

    /**
     * 指定されたパラメータから Directed + Nondirected CM を生成する。
     * 入次数はパワーロー分布、出次数は総数をランダムに配分する。
     *
     * @param name    グラフ名
     * @param n       頂点数
     * @param kInMin  最小入次数
     * @param kInMax  最大入次数
     * @param kuAve   平均無向次数
     * @param gamma   パワーロー指数
     * @param seed    乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(String name, int n, int kInMin, int kInMax, double kuAve, double gamma, long seed) {
        validateGenerateParams(name, n, kInMin, kuAve, gamma);

        Random random = new Random(seed);

        int[] ki = sampleInDegreesPowerLaw(n, kInMin, kInMax, gamma, random);
        int[] ko = assignOutDegreesRandomly(n, (int) sum(ki), random);
        int[] kn = sampleNondirectedDegrees(n, kuAve, random);
        ensureEvenNondirectedSum(kn, random);

        return generateFromDegreeSequence(name, ki, ko, kn, seed ^ SEED_MIX_DIRECTED);
    }

    // --- validation ---

    private static void validateDegreeSequence(String name, int[] ki, int[] ko, int[] kn) {
        if (name == null) {
            throw new IllegalArgumentException("name must be non-null");
        }
        if (ki == null || ko == null || kn == null) {
            throw new IllegalArgumentException("Degree arrays must be non-null");
        }
        if (ki.length != ko.length || ki.length != kn.length) {
            throw new IllegalArgumentException("Degree arrays must have the same length");
        }
        int n = ki.length;
        long sumKi = 0, sumKo = 0, sumKn = 0;
        for (int u = 0; u < n; u++) {
            if (ki[u] < 0 || ko[u] < 0 || kn[u] < 0) {
                throw new IllegalArgumentException("degree values must be non-negative");
            }
            sumKi += ki[u];
            sumKo += ko[u];
            sumKn += kn[u];
        }
        if (sumKi != sumKo) {
            throw new IllegalArgumentException("Sum of in-degrees must equal sum of out-degrees");
        }
        if ((sumKn & 1L) != 0L) {
            throw new IllegalArgumentException("Sum of nondirected degrees must be even");
        }
        if (sumKo > Integer.MAX_VALUE || sumKn > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Sum of degrees must be less than or equal to Integer.MAX_VALUE");
        }
    }

    private static void validateGenerateParams(String name, int n, int kInMin, double kuAve, double gamma) {
        if (name == null) throw new IllegalArgumentException("name must be non-null");
        if (n < 0) throw new IllegalArgumentException("n must be non-negative");
        if (kInMin < 0) throw new IllegalArgumentException("kInMin must be non-negative");
        if (kuAve < 0) throw new IllegalArgumentException("kuAve must be non-negative");
        if (gamma <= 1.0) {
            throw new IllegalArgumentException("gamma must be greater than 1.0 for power-law distribution");
        }
    }

    // --- degree sequence helpers ---

    private static long sum(int[] arr) {
        long s = 0;
        for (int x : arr) s += x;
        return s;
    }

    /**
     * パワーロー分布に従う入次数列をサンプルする。
     */
    private static int[] sampleInDegreesPowerLaw(int n, int kInMin, int kInMax, double gamma, Random random) {
        double[] cdf = buildPowerLawCdf(kInMin, kInMax, gamma);
        int[] ki = new int[n];
        for (int u = 0; u < n; u++) {
            double r = random.nextDouble();
            int idx = Arrays.binarySearch(cdf, r);
            ki[u] = idx >= 0 ? idx + kInMin : -idx - 1 + kInMin;
        }
        return ki;
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
        cdf[len - 1] = 1.0;
        return cdf;
    }

    /**
     * 総出次数 totalOut を n ノードにランダムに配分する。
     */
    private static int[] assignOutDegreesRandomly(int n, int totalOut, Random random) {
        int[] ko = new int[n];
        for (int i = 0; i < totalOut; i++) {
            ko[random.nextInt(n)]++;
        }
        return ko;
    }

    /**
     * 平均 kuAve のポアソン分布で無向次数列をサンプルする。
     */
    private static int[] sampleNondirectedDegrees(int n, double kuAve, Random random) {
        int[] kn = new int[n];
        for (int u = 0; u < n; u++) {
            kn[u] = RandomUtils.getPoisson(kuAve, random);
        }
        return kn;
    }

    /**
     * 無向次数の合計を偶数にする（nondirected のペアリングのため）。kn を破壊的に更新する。
     */
    private static void ensureEvenNondirectedSum(int[] kn, Random random) {
        long sumKn = sum(kn);
        if ((sumKn & 1L) == 0L) return;
        int n = kn.length;
        for (int attempt = 0; attempt < n * 2; attempt++) {
            int u = random.nextInt(n);
            if (kn[u] > 0) {
                kn[u]--;
                return;
            }
        }
        kn[random.nextInt(n)]++;
    }

    // --- stub and edge construction ---

    /**
     * 有向用の out-stub / in-stub 配列を構築する。返り値は { outStubs, inStubs }。
     */
    private static int[][] buildDirectedStubs(int n, int[] ki, int[] ko, int mDir) {
        int[] outStubs = new int[mDir];
        int[] inStubs = new int[mDir];
        int pOut = 0, pIn = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < ko[u]; t++) outStubs[pOut++] = u;
            for (int t = 0; t < ki[u]; t++) inStubs[pIn++] = u;
        }
        if (pOut != mDir || pIn != mDir) {
            throw new IllegalStateException("Stub count mismatch: pOut=" + pOut + ", pIn=" + pIn + ", mDir=" + mDir);
        }
        return new int[][]{ outStubs, inStubs };
    }

    /**
     * 無向用の stub 配列を構築する（長さ mUnd*2）。ペア (0,1), (2,3), ... でエッジになる。
     */
    private static int[] buildUndirectedStubs(int n, int[] kn, int mUnd) {
        int totalStubs = mUnd * 2;
        int[] undStubs = new int[totalStubs];
        int p = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < kn[u]; t++) undStubs[p++] = u;
        }
        if (p != totalStubs) {
            throw new IllegalStateException("Undirected stub count mismatch: p=" + p + ", expected=" + totalStubs);
        }
        return undStubs;
    }

    /**
     * 有向・無向スタブから base edge 配列 (srcs, dsts, isUndirected) を埋める。
     */
    private static void fillBaseEdgeArrays(int[] srcs, int[] dsts, boolean[] isUndirected,
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
