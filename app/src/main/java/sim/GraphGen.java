package sim;

import sim.network.DirectedGraph;
import sim.network.topology.EgoTwitter;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class GraphGen {
    public static void main(String[] args) throws IOException {
        int n = 500_000;
        int kMin = 7;
        int kMax = 1000;
        double gamma = 4.0;
        int itrs = 1;

        for (int itr = 0; itr < itrs; itr++) {

            DirectedGraph g = EgoTwitter.loadFromDefaultEdgeList();
            System.out.println("");
            System.out.println("--------------------------------");
            g.printInfo();
            System.out.println("--------------------------------");
            System.out.println("");

            // パス構成: out/edgelist/{NetworkPath}/{filename}
            Path networkPath = SwitchUtils.buildNetworkPath("ego-Twitter", g.n,
                    null, null, null, null, null, null, null, null,
                    kMin, kMax, kMin, kMax,
                    gamma, null);
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
