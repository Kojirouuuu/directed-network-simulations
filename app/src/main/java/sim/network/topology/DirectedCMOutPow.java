package sim.network.topology;

import sim.network.DirectedGraph;
import sim.utils.ArrayUtils;
import sim.utils.RandomUtils;

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

public class DirectedCMOutPow {
    private DirectedCMOutPow() {}

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
        ArrayUtils.shuffle(inStubs, seed ^ 0x9E3779B97F4A7C15L);

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

        ArrayUtils.shuffle(undStubs, seed ^ 0x517CC1B727220A95L);

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
    public static DirectedGraph generate(String name, int n, int kHat, double gamma, long seed) {
        if (name == null) throw new IllegalArgumentException("name must be non-null");
        if (n < 0) throw new IllegalArgumentException("n must be non-negative");
        if (kHat < 0) throw new IllegalArgumentException("kHat must be non-negative");
        if (gamma <= 1.0) {
            throw new IllegalArgumentException("gamma must be greater than 1.0 for power-law distribution");
        }

        Random random = new Random(seed);

        int[] ki = new int[n];
        int[] ko = new int[n];
        int[] kn = new int[n];

        // ki[u] の生成部分を修正
        // 入次数 (Power Law) の生成
        double kInMin = 1.0; 
        double kInMax = n - 1; 

        // 逆変換法のための定数計算
        double kInMinPow = Math.pow(kInMin, 1.0 - gamma);
        double kInMaxPow = Math.pow(kInMax, 1.0 - gamma);
        double oneOverOneMinusGamma = 1.0 / (1.0 - gamma);

        // ユーザーの式の係数を計算
        double nMinus1Pow2MinusGamma = Math.pow(n - 1, 2.0 - gamma);
        double coefficient = (2.0 - gamma) / (1.0 - nMinus1Pow2MinusGamma);

        // まず、正規化前の値を生成
        double[] koRaw = new double[n];
        double sumKoRaw = 0.0;
        for (int u = 0; u < n; u++) {
            double r = random.nextDouble();
            
            // 逆関数法による Power Law 乱数 x の生成
            double x = Math.pow((kInMaxPow - kInMinPow) * r + kInMinPow, oneOverOneMinusGamma);
            x = Math.max(1.0, Math.min(n - 1, x));
            
            // ユーザーの式に従って ki[u] を計算
            koRaw[u] = coefficient * kHat * Math.pow(x, -gamma);
            sumKoRaw += koRaw[u];
        }

        // 平均が kHat になるように正規化（目標合計: n * kHat）
        long targetSum = (long) n * kHat;
        double scale = targetSum / sumKoRaw;
        
        long sumKo = 0;
        for (int u = 0; u < n; u++) {
            ko[u] = (int) Math.round(koRaw[u] * scale);
            if (ko[u] < 1) ko[u] = 1;
            if (ko[u] > n - 1) ko[u] = n - 1;
            sumKo += ko[u];
        }

        // 丸め誤差を調整して、合計が targetSum になるようにする
        long deltaKo = targetSum - sumKo;
        while (deltaKo != 0) {
            int u = random.nextInt(n);
            if (deltaKo > 0) {
                if (ko[u] < n - 1) {
                    ko[u]++;
                    deltaKo--;
                    sumKo++;
                }
            } else {
                if (ko[u] > 1) {
                    ko[u]--;
                    deltaKo++;
                    sumKo--;
                }
            }
        }

        // sum(ko) を sum(ki) に合わせる（分布の歪みは O(sqrt(n)) 程度で、n が大きいほど無視できる）
        long sumKi = 0;
        for (int u = 0; u < n; u++) {
            ki[u] = RandomUtils.getPoisson(kHat, random);
            sumKi += ki[u];
        }

        long delta = sumKi - sumKo;

        // delta>0: ko を減らす必要がある / delta<0: ko を増やす必要がある
        while (delta != 0) {
            int u = random.nextInt(n);
            if (delta > 0) {
                if (ki[u] > 0) {
                    ki[u]--;
                    delta--;
                }
            } else {
                ki[u]++;
                delta++;
            }
        }

        // 無向次数 (Poisson) の生成
        long sumKn = 0;
        for (int u = 0; u < n; u++) {
            kn[u] = RandomUtils.getPoisson(kHat, random);
            sumKn += kn[u];
        }

        // sumKn が偶数になるように調整（nondirected はペアリングするために偶数が必要）
        if ((sumKn & 1L) != 0L) {
            // 奇数なので1つ減らす（0より大きいノードを探す）
            boolean adjusted = false;
            for (int attempt = 0; attempt < n * 2 && !adjusted; attempt++) {
                int u = random.nextInt(n);
                if (kn[u] > 0) {
                    kn[u]--;
                    sumKn--;
                    adjusted = true;
                }
            }
            // すべてのノードが0の場合（理論上は起こりにくいが）
            if (!adjusted) {
                int u = random.nextInt(n);
                kn[u]++;
                sumKn++;
            }
        }

        return generateFromDegreeSequence(name, ki, ko, kn, seed ^ 0x9E3779B97F4A7C15L);
    }
}
