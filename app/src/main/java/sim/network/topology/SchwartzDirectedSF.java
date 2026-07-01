package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;

import java.util.Arrays;
import java.util.Random;

/**
 * Schwartz et al. 2002 "Percolation in Directed Scale-Free Networks" で提案された
 * 相関を持つ有向スケールフリーネットワークを生成するクラス。
 *
 * <p>
 * λ_in と λ_out の異なるべき則指数を持ち、確率 A で
 * K = round(J^{(λ_in-1)/(λ_out-1)}) という決定論的写像で入出次数を相関させる。
 * これは PowPow の事後スワップ相関とは異なる相関注入方式である。
 *
 * <p>
 * 本実装は論文 Eq. (8) の相関構造を、有限サイズ・有限カットオフの
 * configuration model として近似的に実装するものである。
 * 具体的には mIn ≥ 1 を要請するため j = 0 のマスを扱わず、
 * また Σ K = Σ J を満たすため K を ±1 単位で後調整する点で
 * Eq. (8) の P(j, k) と厳密には異なる。
 *
 * <p>
 * Configuration Model のセマンティクス通り、自己ループおよび多重辺は許容される。
 */
public final class SchwartzDirectedSF {

    private static final long SEED_MIX_IN = 0x9E3779B97F4A7C15L;
    private static final long SEED_MIX_OUT_SAMPLE = 0x517CC1B727220A95L;
    private static final long SEED_MIX_OUT_CORR = 0xC2B2AE3D27D4EB4FL;
    private static final long SEED_MIX_BALANCE = 0xDEADBEEFCAFEBABEL;
    private static final long SEED_MIX_SHUFFLE_OUT = 0x1234567890ABCDEFL;
    private static final long SEED_MIX_SHUFFLE_IN = 0xFEDCBA0987654321L;

    private SchwartzDirectedSF() {
    }

    // =========================================================
    // Public API
    // =========================================================

    /**
     * Schwartz らの相関有向スケールフリーネットワークを生成する。
     *
     * <p>
     * アルゴリズム概要:
     * <ol>
     * <li>入次数 J[i] をべき則 P(j) ∝ j^{-λ_in}, j ∈ [mIn, kInMax] からサンプリング。</li>
     * <li>出次数 K[i] を確率 A で K = clamp(round(J^{(λ_in-1)/(λ_out-1)}), [mOut,
     * kOutMax]) と決定論的に与え、
     * 確率 1-A で P(k) ∝ k^{-λ_out}, k ∈ [mOut, kOutMax] から独立にサンプリング。</li>
     * <li>Σ K = Σ J になるよう K を ±1 単位で調整（balanceSums）。</li>
     * <li>スタブをシャッフルしてペアリング（自己ループ・多重辺は許容）。</li>
     * </ol>
     *
     * @param name グラフ名
     * @param n 頂点数 (>= 1)
     * @param mIn 入次数の最小値 (>= 1)
     * @param kInMax 入次数の最大値 (>= mIn)
     * @param mOut 出次数の最小値 (>= 1)
     * @param kOutMax 出次数の最大値 (>= mOut)
     * @param gammaIn 入次数べき則指数 λ_in (> 1.0)
     * @param gammaOut 出次数べき則指数 λ_out (> 1.0)
     * @param corrA 入出次数の相関確率 A (0.0 <= A <= 1.0)
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(
            String name, int n,
            int mIn, int kInMax, int mOut, int kOutMax,
            double gammaIn, double gammaOut, double corrA,
            long seed) {

        validate(name, n, mIn, kInMax, mOut, kOutMax, gammaIn, gammaOut, corrA);

        // Step 1: 入次数 J を λ_in のべき則からサンプリング
        int[] ji = samplePowerLaw(n, mIn, kInMax, gammaIn, seed ^ SEED_MIX_IN);

        // Step 2: 出次数 K を確率 A で J と相関、そうでなければ λ_out のべき則から独立にサンプリング
        int[] ki = sampleOutDegreesWithCorrelation(
                ji, mOut, kOutMax, gammaIn, gammaOut, corrA,
                seed ^ SEED_MIX_OUT_SAMPLE, seed ^ SEED_MIX_OUT_CORR);

        // Step 3: Σ K = Σ J となるよう K を調整（PowPow.balanceSums と同形だが K 側を調整）
        balanceOutSumToInSum(ki, ji, mOut, kOutMax, seed ^ SEED_MIX_BALANCE);

        // Step 4: スタブをシャッフルしてペアリング → 辺配列を構築
        return buildGraphFromDegreeSequence(name, n, ji, ki, seed);
    }

    // =========================================================
    // Validation
    // =========================================================

    private static void validate(
            String name, int n, int mIn, int kInMax, int mOut, int kOutMax,
            double gammaIn, double gammaOut, double corrA) {
        if (name == null)
            throw new IllegalArgumentException("name must be non-null");
        if (n < 1)
            throw new IllegalArgumentException("n must be >= 1");
        if (mIn < 1)
            throw new IllegalArgumentException("mIn must be >= 1");
        if (mOut < 1)
            throw new IllegalArgumentException("mOut must be >= 1");
        if (kInMax < mIn)
            throw new IllegalArgumentException("kInMax must be >= mIn");
        if (kOutMax < mOut)
            throw new IllegalArgumentException("kOutMax must be >= mOut");
        if (gammaIn <= 1.0)
            throw new IllegalArgumentException("gammaIn must be > 1.0");
        if (gammaOut <= 1.0)
            throw new IllegalArgumentException("gammaOut must be > 1.0");
        if (corrA < 0.0 || corrA > 1.0)
            throw new IllegalArgumentException("corrA must be in [0.0, 1.0]");
    }

    // =========================================================
    // Degree sequence helpers
    // =========================================================

    /**
     * べき則分布 P(k) ∝ k^{-gamma} に従う次数列を [kMin, kMax] からサンプリングする。
     * PowPow.samplePowerLawDistribution と同一の CDF 逆変換による実装。
     */
    private static int[] samplePowerLaw(int n, int kMin, int kMax, double gamma, long seed) {
        Random rng = new Random(seed);
        double[] cdf = buildPowerLawCdf(kMin, kMax, gamma);
        int[] k = new int[n];
        for (int u = 0; u < n; u++) {
            double r = rng.nextDouble();
            int idx = Arrays.binarySearch(cdf, r);
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
        cdf[len - 1] = 1.0;
        return cdf;
    }

    /**
     * 各ノードの出次数を決定する。
     * 入次数 ji[u] > 0 かつ U(0,1) < corrA のとき K[u] = clamp(round(ji[u]^exponent),
     * [mOut, kOutMax])、
     * そうでなければ λ_out のべき則から独立にサンプリングする。
     *
     * <p>
     * 効率のため、独立サンプリング用の CDF は corrA < 1.0 の場合のみ事前構築する。
     *
     * @param ji 入次数列（変更しない）
     * @param mOut 出次数の最小値
     * @param kOutMax 出次数の最大値
     * @param gammaIn λ_in
     * @param gammaOut λ_out
     * @param corrA 相関確率 A
     * @param sampleSeed 独立サンプリング用シード
     * @param decisionSeed 相関判定用シード
     * @return 出次数列
     */
    private static int[] sampleOutDegreesWithCorrelation(
            int[] ji, int mOut, int kOutMax,
            double gammaIn, double gammaOut, double corrA,
            long sampleSeed, long decisionSeed) {
        int n = ji.length;
        int[] ki = new int[n];

        // 相関写像の指数: K^(λ_out - 1) = J^(λ_in - 1) より K = J^{(λ_in - 1)/(λ_out - 1)}
        double exponent = (gammaIn - 1.0) / (gammaOut - 1.0);

        // 独立サンプリング用 CDF（必要なときのみ構築）
        double[] outCdf = (corrA < 1.0) ? buildPowerLawCdf(mOut, kOutMax, gammaOut) : null;

        Random decisionRng = new Random(decisionSeed);
        Random sampleRng = new Random(sampleSeed);

        for (int u = 0; u < n; u++) {
            boolean correlated = (ji[u] > 0) && (decisionRng.nextDouble() < corrA);
            if (correlated) {
                long mapped = Math.round(Math.pow(ji[u], exponent));
                if (mapped < mOut) {
                    mapped = mOut;
                } else if (mapped > kOutMax) {
                    mapped = kOutMax;
                }
                ki[u] = (int) mapped;
            } else {
                double r = sampleRng.nextDouble();
                int idx = Arrays.binarySearch(outCdf, r);
                ki[u] = (idx >= 0 ? idx : -idx - 1) + mOut;
            }
        }
        return ki;
    }

    /**
     * Σ K[i] == Σ J[i] となるよう K（出次数列）を ±1 単位で調整する。
     * 擬似コード prompt.md L41-L51 に一致。
     *
     * <p>
     * まず sum(ji) が K の実現可能範囲 [n·mOut, n·kOutMax] に入っているかを検査し、
     * 範囲外なら {@link IllegalArgumentException} を投げる。
     *
     * <p>
     * 調整自体は、現時点で +1（または -1）可能なノードのプールを最初に走査して作成し、
     * そこからランダムに選んで増減する。境界に達したノードは swap-and-pop でプールから外す。
     * 計算量は O(n + |sum(ji) - sum(ki)|) で、有効候補が少ないケースでも線形時間で完了する。
     *
     * @throws IllegalArgumentException sum(ji) が [n·mOut, n·kOutMax] の外側のとき
     */
    private static void balanceOutSumToInSum(
            int[] ki, int[] ji, int mOut, int kOutMax, long seed) {
        final int n = ki.length;
        final long sumJi = sum(ji);
        long sumKi = sum(ki);

        // 実現可能性チェック: K[i] ∈ [mOut, kOutMax] という制約のもとで sum(ki) = sum(ji) にできるか
        // long で乗算することで int オーバーフローを回避する
        final long minPossible = (long) n * mOut;
        final long maxPossible = (long) n * kOutMax;
        if (sumJi < minPossible || sumJi > maxPossible) {
            throw new IllegalArgumentException(
                    "sum(ji)=" + sumJi + " is outside feasible range ["
                            + minPossible + ", " + maxPossible + "]"
                            + " for K with n=" + n + ", mOut=" + mOut + ", kOutMax=" + kOutMax);
        }

        Random rng = new Random(seed);

        if (sumKi < sumJi) {
            // K が不足 → +1 可能なノードを集めて、そこからランダムに選んで増分
            int[] pool = new int[n];
            int poolSize = 0;
            for (int u = 0; u < n; u++) {
                if (ki[u] < kOutMax)
                    pool[poolSize++] = u;
            }
            while (sumKi < sumJi) {
                if (poolSize == 0) {
                    // 実現可能性チェック後は到達しないはず（防御的）
                    throw new IllegalStateException(
                            "No nodes available to increment K (should not happen after feasibility check)");
                }
                int p = rng.nextInt(poolSize);
                int u = pool[p];
                ki[u]++;
                sumKi++;
                if (ki[u] >= kOutMax) {
                    // 上限に達したら swap-and-pop で候補から外す
                    pool[p] = pool[--poolSize];
                }
            }
        } else if (sumKi > sumJi) {
            // K が過剰 → -1 可能なノードを集めて、そこからランダムに選んで減分
            int[] pool = new int[n];
            int poolSize = 0;
            for (int u = 0; u < n; u++) {
                if (ki[u] > mOut)
                    pool[poolSize++] = u;
            }
            while (sumKi > sumJi) {
                if (poolSize == 0) {
                    throw new IllegalStateException(
                            "No nodes available to decrement K (should not happen after feasibility check)");
                }
                int p = rng.nextInt(poolSize);
                int u = pool[p];
                ki[u]--;
                sumKi--;
                if (ki[u] <= mOut) {
                    pool[p] = pool[--poolSize];
                }
            }
        }
        // sumKi == sumJi のときは何もしない
    }

    private static long sum(int[] arr) {
        long s = 0;
        for (int x : arr)
            s += x;
        return s;
    }

    // =========================================================
    // Edge construction
    // =========================================================

    /**
     * (ji, ki) からスタブを構築し、シャッフル後にペアリングして DirectedGraph を返す。
     * 自己ループ・多重辺は許容する。
     *
     * @param name グラフ名
     * @param n 頂点数
     * @param ji 入次数列
     * @param ki 出次数列（Σ ki == Σ ji が成立していること）
     * @param seed 乱数シード（内部で複数の用途にミックスする）
     * @return DirectedGraph
     */
    private static DirectedGraph buildGraphFromDegreeSequence(
            String name, int n, int[] ji, int[] ki, long seed) {
        // int キャスト前にオーバーフロー検査（極端に密なグラフで sum(ki) が Integer.MAX_VALUE を超えうる）
        long edgeCount = sum(ki);
        if (edgeCount > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Too many directed edges: " + edgeCount);
        }
        final int mDir = (int) edgeCount;
        if (sum(ji) != edgeCount) {
            throw new IllegalStateException(
                    "Internal error: sum(ji) != sum(ki) after balancing: sumJi="
                            + sum(ji) + ", sumKi=" + edgeCount);
        }

        // 出スタブ・入スタブを構築
        int[] outStubs = new int[mDir];
        int[] inStubs = new int[mDir];
        int pOut = 0, pIn = 0;
        for (int u = 0; u < n; u++) {
            for (int t = 0; t < ki[u]; t++)
                outStubs[pOut++] = u;
            for (int t = 0; t < ji[u]; t++)
                inStubs[pIn++] = u;
        }
        if (pOut != mDir || pIn != mDir) {
            throw new IllegalStateException(
                    "Directed stub count mismatch: pOut=" + pOut + ", pIn=" + pIn + ", mDir=" + mDir);
        }

        ArrayUtils.shuffle(outStubs, seed ^ SEED_MIX_SHUFFLE_OUT);
        ArrayUtils.shuffle(inStubs, seed ^ SEED_MIX_SHUFFLE_IN);

        // 辺配列を生成（無向辺なし）
        int[] srcs = new int[mDir];
        int[] dsts = new int[mDir];
        boolean[] isUndirected = new boolean[mDir];
        for (int i = 0; i < mDir; i++) {
            srcs[i] = outStubs[i];
            dsts[i] = inStubs[i];
            isUndirected[i] = false;
        }

        return DirectedGraph.fromEdgeListWithUndirectedFlag(name, n, srcs, dsts, isUndirected);
    }
}
