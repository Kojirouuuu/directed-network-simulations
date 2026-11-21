package sim.simulation;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.network.DirectedGraph;

class BondPercolationResultTest {
    @TempDir Path temp;

    @Test
    void writesSortedFrequenciesIncludingLargestTiesAndIsolatesAndAppendsOneHeader() throws IOException {
        var graph = DirectedGraph.fromEdgeListWithUndirectedFlag("undirected", 7,
                new int[] {0, 1, 3, 4}, new int[] {1, 2, 4, 5}, new boolean[] {true, true, true, true});
        var simulator = new BondPercolationSimulator(graph);
        Path csv = temp.resolve("nested/sizes.csv");
        simulator.run(1, new Random(42)).writeComponentSizeDistributionCsv(csv, 3, 1, true);
        simulator.run(0, new Random(42)).writeComponentSizeDistributionCsv(csv, 4, 0, true);
        assertEquals(List.of("itr,p,component_type,size,count",
                "3,1.000000000,wcc,1,1", "3,1.000000000,wcc,3,2",
                "3,1.000000000,scc,1,1", "3,1.000000000,scc,3,2",
                "4,0.000000000,wcc,1,7", "4,0.000000000,scc,1,7"), Files.readAllLines(csv));
    }

    @Test
    void mixedAndDirectedDistributionsAgreeWithSummaryStatistics() throws IOException {
        for (boolean mixed : new boolean[] {false, true}) {
            var graph = DirectedGraph.fromEdgeListWithUndirectedFlag("fixture", 5,
                    new int[] {0, 1, 2}, new int[] {1, 2, 3}, new boolean[] {mixed, false, mixed});
            for (double p : new double[] {0, 0.5, 1}) {
                var result = new BondPercolationSimulator(graph).run(p, new Random(42));
                Path csv = temp.resolve("sizes.csv");
                Files.deleteIfExists(csv);
                result.writeComponentSizeDistributionCsv(csv, 0, p, true);
                List<String> rows = Files.readAllLines(csv);
                for (String type : new String[] {"wcc", "scc"}) {
                    int vertices = 0, components = 0, largest = 0, second = 0;
                    for (String line : rows.subList(1, rows.size())) {
                        String[] fields = line.split(",");
                        if (!fields[2].equals(type)) continue;
                        int size = Integer.parseInt(fields[3]), count = Integer.parseInt(fields[4]);
                        vertices += size * count;
                        components += count;
                        second = count > 1 ? size : largest;
                        largest = size;
                    }
                    assertEquals(5, vertices);
                    assertEquals(type.equals("wcc") ? result.gwccSize : result.gsccSize, largest);
                    assertEquals(type.equals("wcc") ? result.gwccSecondSize : result.gsccSecondSize, second);
                    assertEquals(type.equals("wcc") ? result.gwccMeanSize : result.gsccMeanSize,
                            (double) vertices / components, 1e-12);
                }
            }
        }
        var mixed = DirectedGraph.fromEdgeListWithUndirectedFlag("fixture", 5,
                new int[] {0, 1, 2}, new int[] {1, 2, 3}, new boolean[] {true, false, true});
        Path known = temp.resolve("known.csv");
        new BondPercolationSimulator(mixed).run(1, new Random(42))
                .writeComponentSizeDistributionCsv(known, 0, 1, true);
        assertEquals(List.of("itr,p,component_type,size,count", "0,1.000000000,wcc,1,1",
                "0,1.000000000,wcc,4,1", "0,1.000000000,scc,1,1", "0,1.000000000,scc,2,2"),
                Files.readAllLines(known));
    }

    @Test
    void emptyGraphWritesOnlyHeaderAndNonAppendPreservesExistingFile() throws IOException {
        var graph = DirectedGraph.fromEdgeListWithUndirectedFlag("empty", 0,
                new int[0], new int[0], new boolean[0]);
        var result = new BondPercolationSimulator(graph).run(0.5, new Random(42));
        Path csv = temp.resolve("sizes.csv");
        Files.writeString(csv, "original");
        result.writeComponentSizeDistributionCsv(csv, 0, 0.5, false);
        assertEquals("original", Files.readString(csv));
        assertEquals(List.of("itr,p,component_type,size,count"),
                Files.readAllLines(temp.resolve("sizes (1).csv")));
    }
}
