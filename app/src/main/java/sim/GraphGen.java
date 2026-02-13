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
        int n = 50_000;
        int kOutMin = 5;
        int kOutMax = (int) Math.pow(n, 0.5);
        double kuAve = 0;
        double[] gammaList = {2.5};
        int itrs = 10;

        for (int itr = 0; itr < itrs; itr++) {

            long seed = 1234567890 + itr * 100;
            for (double gamma : gammaList) {
                String outputDir = String.format("out/edgelist/n=%d/gamma=%.2f/kOutMin=%d/kOutMax=%d/", n, gamma, kOutMin, kOutMax);
                DirectedGraph g = DirectedCMOutPow.generate("DirectedCMOutPow", n, kOutMin, kOutMax, kuAve, gamma, seed);
                System.out.println("");
                System.out.println("--------------------------------");
                g.printInfo();
                System.out.println("--------------------------------");
                System.out.println("");
                String fileName = String.format("%s_%d.txt", g.name, itr);
                Path outputPath = Paths.get(outputDir + fileName);
                try {
                    g.writeEdgeList(outputPath);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
