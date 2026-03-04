package sim;

import sim.network.DirectedGraph;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.IntStream;

public class GraphGen {
    /** 並列ワーカー数（コア数） */
    private static final int batchSize = 20;

    /** 各コアが作成するネットワーク数 */
    private static final int n = 500_000;
    private static final int kdMin = 5;
    private static final int kdMax = 1000;
    private static final double gamma = 2.5;
    private static final int swapNum = 0;
    private static final int itr = 1;
    private static final long seed = 42L;

    public static void main(String[] args) throws IOException {
        IntStream.range(0, batchSize).parallel().forEach(workerId -> {
            for (int k = 0; k < itr; k++) {
                int runIndex = workerId + k * batchSize;
                try {
                    DirectedGraph g = SwitchUtils.generateGraph("PowPow", n,
                        null, kdMin, kdMax, null, null, null, null, null, null, null, gamma, null, null, swapNum, seed + runIndex);
                    synchronized (System.out) {
                        System.out.println("");
                        System.out.println("--------------------------------");
                        System.out.println("[worker " + workerId + " run " + runIndex + "]");
                        g.printInfo();
                        System.out.println("--------------------------------");
                        System.out.println("");
                    }

                    // パス構成: out/edgelist/{NetworkPath}/{filename}
                    Path networkPath = SwitchUtils.buildNetworkPath(g.name, g.n,
                            null, null, null, null, null, null,
                            kdMin, kdMax, null, null, null, null,
                            gamma, swapNum);
                    Path outputDir = Paths.get("out/edgelist").resolve(networkPath);
                    String fileName = String.format("%d.csv", runIndex);
                    Path outputPath = outputDir.resolve(fileName);

                    g.writeEdgeList(outputPath);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        });
    }
}
