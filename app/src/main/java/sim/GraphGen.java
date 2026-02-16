package sim;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedCM;
import sim.network.topology.DirectedCMInPow;
import sim.network.topology.DirectedCMOutPow;
import sim.network.topology.BA;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.IOException;

public class GraphGen {
    public static void main(String[] args) {
        int n = 1_000_000;
        int m0 = 6;
        int m = 5;
        boolean isDirected = true;
        int itrs = 10;

        for (int itr = 0; itr < itrs; itr++) {

        long seed = 1234567890 + itr * 100;
            String outputDir = String.format("out/edgelist/n=%d/m0=%d/m=%d/isDirected=%b/", n, m0, m, isDirected);
            DirectedGraph g = BA.generate("BA", n, m0, m, isDirected, seed);
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
