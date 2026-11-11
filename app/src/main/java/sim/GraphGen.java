package sim;

import sim.network.DirectedGraph;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;
import java.util.stream.IntStream;

public class GraphGen {
    // 実行制御
    /** true: batchSize 個のワーカーで並列生成、false: 逐次生成 */
    private static final boolean parallelGeneration = false;
    private static final boolean randomizeByEdgeSwaps = true; // 書き出し前に次数保存ランダム辺スワップを行うか
    /** 並列ワーカー数（コア数） */
    private static final int batchSize = 10;
    /** 各コアが作成するネットワーク数 */
    private static final int itr = 1;

    // ネットワークの基本設定
    private static final String networkType = "rev-higgs-social"; // ネットワークタイプ
    private static final int n = 500_000;

    // 次数パラメータ
    private static final int kMin = 5;
    private static final int kMax = (int) Math.sqrt(n);
    private static final int kInMin = 5;
    private static final int kInMax = (int) Math.sqrt(n);
    private static final int kOutMin = 5;
    private static final int kOutMax = (int) Math.sqrt(n);
    private static final double kuAve = 10.0;
    private static final int kuMin = 5;
    private static final int kuMax = (int) Math.sqrt(n);

    // トポロジー生成パラメータ
    private static final int m0 = 0;
    private static final int m = 0;
    private static final double gamma = 2.5;
    private static final int swapNum = 0;

    // SchwartzDirectedSF 用パラメータ
    private static final double gammaIn = 2.5;
    private static final double gammaOut = 2.5;
    private static final double corrA = 0.0;

    // TwoPointDirectedCM: set all four when selecting this topology.
    private static final Integer k1 = null, k2 = null;
    private static final Double p = null, r = null;

    // リワイヤリング・乱数パラメータ
    private static final double targetReciprocity = 0.00;
    private static final long seed = 42L;
    private static final long randomSwapSeed = 4242L;

    public static void main(String[] args) throws IOException {
        if (parallelGeneration) {
            IntStream.range(0, batchSize).parallel().forEach(GraphGen::runWorker);
        } else {
            for (int workerId = 0; workerId < batchSize; workerId++) {
                runWorker(workerId);
            }
        }
    }

    private static void runWorker(int workerId) {
        for (int k = 0; k < itr; k++) {
            int runIndex = workerId + k * batchSize;
            generateAndWrite(workerId, runIndex);
        }
    }

    private static void generateAndWrite(int workerId, int runIndex) {
        try {
            var parameters = SwitchUtils.GraphGeneratorParams.builder(networkType, n)
                    .kdMin(kMin)
                    .kdMax(kMax)
                    .kInMin(kInMin)
                    .kInMax(kInMax)
                    .kOutMin(kOutMin)
                    .kOutMax(kOutMax)
                    .kuMin(kuMin)
                    .kuMax(kuMax)
                    .kuAve(kuAve)
                    .gamma(gamma)
                    .m0(m0)
                    .m(m)
                    .swapNum(swapNum)
                    .gammaIn(gammaIn)
                    .gammaOut(gammaOut)
                    .corrA(corrA).k1(k1).k2(k2).p(p).r(r).build();
            DirectedGraph g = SwitchUtils.generateGraph(parameters, seed + runIndex);

            if (randomizeByEdgeSwaps) {
                g = g.randomizeByEdgeSwaps(randomSwapSeed + runIndex);
            }

            // パス構成: out/edgelist/{NetworkPath}/randomization={mode}/{filename}
            Path networkPath = SwitchUtils.buildNetworkPath(parameters.withN(g.n));
            Path outputDir = Paths.get("out/edgelist")
                    .resolve(SwitchUtils.appendRandomizationPath(networkPath, randomizeByEdgeSwaps))
                    .resolve(String.format(Locale.ROOT, "reciprocity=%.6f", targetReciprocity));
            String fileName = String.format("%d.csv", runIndex);
            Path outputPath = outputDir.resolve(fileName);

            g.writeEdgeList(outputPath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
