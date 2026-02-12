package sim;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedCM;
import sim.network.topology.DirectedCMInPow;
import sim.network.topology.DirectedCMOutPow;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;

public class GraphGen {
    public static void main(String[] args) {
        int n = 1_000_0000;
        int kOutMin = 3;
        int kOutMax = n - 1;
        double kuAve = 6.2;
        double[] gammaList = {2.7};
        int itrs = 10;

        for (int itr = 0; itr < itrs; itr++) {
            // String outputDir = String.format("out/edgelist/n=%d/kHat=%d/", n, kHat);

            long seed = 1234567890 + itr * 100;
            for (double gamma : gammaList) {
                DirectedGraph g = DirectedCMInPow.generate("DirectedCMInPow", n, kOutMin, kOutMax, kuAve, gamma, seed);
                System.out.println("");
                System.out.println("--------------------------------");
                g.printInfo();
                System.out.println("--------------------------------");
                System.out.println("");
                // String fileName = String.format("%s_%d.txt", g.name, itr);
                // Path outputPath = Paths.get(outputDir + fileName);
                // try {
                //     g.writeEdgeList(outputPath);
                // } catch (IOException e) {
                //     e.printStackTrace();
                // }
            }
        }
    }
}
