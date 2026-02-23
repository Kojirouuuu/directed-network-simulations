package sim;

import sim.network.DirectedGraph;
import sim.network.topology.PowPow;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class GraphGen {
    public static void main(String[] args) {
        int n = 500_000;
        int kMin = 5;
        int kMax = (int) Math.pow(n, 0.5);
        double gamma = 2.5;
        int itrs = 20;
        int swapNum = 10000;

        for (int itr = 0; itr < itrs; itr++) {

            long seed = 1234567890 + itr * 100;
            DirectedGraph g = PowPow.generate("PowPow", n, kMin, kMax, gamma, swapNum, seed);
            System.out.println("");
            System.out.println("--------------------------------");
            g.printInfo();
            System.out.println("--------------------------------");
            System.out.println("");

            // パス構成: out/edgelist/{NetworkPath}/{filename}
            Path networkPath = SwitchUtils.buildNetworkPath("PowPow", n,
                    null, null, null, null, null, null,
                    kMin, kMax, null, null,
                    gamma, swapNum);
            Path outputDir = Paths.get("out/edgelist").resolve(networkPath);
            String fileName = String.format("%s_%d.csv", g.name, itr);
            Path outputPath = outputDir.resolve(fileName);

            try {
                g.writeEdgeList(outputPath);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
