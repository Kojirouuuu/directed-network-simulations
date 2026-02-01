package sim;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedCM;
import sim.network.topology.DirectedCMInPow;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;

public class GraphGen {
    public static void main(String[] args) {
        int n = 1_000_0000;
        int[] kHatList = {4};
        double[] gammaList = {2.4};
        int itrs = 100;

        for (int itr = 0; itr < itrs; itr++) {
            for (int kHat : kHatList) {
                // String outputDir = String.format("out/edgelist/n=%d/kHat=%d/", n, kHat);

                long seed = 1234567890 + itr * 100 + kHat;
                for (double gamma : gammaList) {
                    DirectedGraph g = DirectedCMInPow.generate("DirectedCMInPow", n, kHat, gamma, seed);
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
}
