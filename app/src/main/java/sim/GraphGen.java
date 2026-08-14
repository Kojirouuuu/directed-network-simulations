package sim;

import sim.network.DirectedGraph;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;
import java.util.stream.IntStream;

public class GraphGen {
    /** true: batchSize 個のワーカーで並列生成、false: 逐次生成 */
    private static final boolean parallelGeneration = false;
    /** 並列ワーカー数（コア数） */
    private static final int batchSize = 10;

    private static final String networkType = "rev-higgs-social"; // ネットワークタイプ
    private static final int n = 500_000;
    private static final int kMin = 5;
    private static final int kMax = (int) Math.sqrt(n);
    private static final double gamma = 2.5;
    private static final int swapNum = 0;
    /** 各コアが作成するネットワーク数 */
    private static final int itr = 1;
    private static final long seed = 42L;
    private static final double targetReciprocity = 0.01;
    private static final int maxIncreaseAttemptsPerEdge = 200;
    private static final int neutralSwapAttemptsPerEdge = 1;

    private static final int kInMin = 5;
    private static final int kInMax = (int) Math.sqrt(n);
    private static final int kOutMin = 5;
    private static final int kOutMax = (int) Math.sqrt(n);
    private static final int kuMin = 5;
    private static final int kuMax = (int) Math.sqrt(n);
    private static final double kuAve = 10.0;
    private static final double gammaIn = 2.5;
    private static final double gammaOut = 2.5;
    private static final double corrA = 0.0;
    private static final int m0 = 0;
    private static final int m = 0;

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
            DirectedGraph g = SwitchUtils.generateGraph(networkType, n,
                    null, kMin, kMax, kInMin, kInMax, kOutMin,
                    kOutMax,
                    kuMin, kuMax,
                    kuAve, gamma, m0, m, swapNum,
                    gammaIn, gammaOut, corrA,
                    seed + runIndex);

            long maxIncreaseAttempts = (long) maxIncreaseAttemptsPerEdge * g.m;
            long neutralSwapAttempts = (long) neutralSwapAttemptsPerEdge * g.m;
            DirectedGraph g2 = g.increaseReciprocity(
                    targetReciprocity, maxIncreaseAttempts, neutralSwapAttempts, seed + runIndex);
            synchronized (System.out) {
                System.out.println("");
                System.out.println("--------------------------------");
                System.out.println("[worker " + workerId + " run " + runIndex + "]");
                g2.printInfo();
                System.out.println("--------------------------------");
                System.out.println("");
            }

            // パス構成: out/edgelist/{NetworkPath}/{filename}
            Path networkPath = SwitchUtils.buildNetworkPath(networkType, g2.n,
                    null, null, null, null, null, null, kMin, kMax, null, null, null, null,
                    gamma, swapNum,
                    null, null, null);
            Path outputDir = Paths.get("out/edgelist")
                    .resolve(networkPath)
                    .resolve(String.format(Locale.ROOT, "reciprocity=%.6f", targetReciprocity));
            String fileName = String.format("%d.csv", runIndex);
            Path outputPath = outputDir.resolve(fileName);

            g2.writeEdgeList(outputPath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
