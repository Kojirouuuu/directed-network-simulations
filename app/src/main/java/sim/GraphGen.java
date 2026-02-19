package sim;

import sim.network.DirectedGraph;
import sim.network.topology.undirected.BA;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class GraphGen {
    public static void main(String[] args) {
        int n = 1_000_000;
        int m0 = 6;
        int m = 5;
        int itrs = 10;

        for (int itr = 0; itr < itrs; itr++) {

            long seed = 1234567890 + itr * 100;
            DirectedGraph g = BA.generate("BA", n, m0, m, seed);
            System.out.println("");
            System.out.println("--------------------------------");
            g.printInfo();
            System.out.println("--------------------------------");
            System.out.println("");

            // パス構成: out/edgelist/{NetworkPath}/{filename}
            Path networkPath = SwitchUtils.buildNetworkPath("BA", n, null, null, null, null, m0, m);
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
