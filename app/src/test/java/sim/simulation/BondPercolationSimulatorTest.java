package sim.simulation;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.Random;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;

class BondPercolationSimulatorTest {
    @Test
    void undirectedComponentsIncludeIsolatedVerticesAndTiedLargestComponents() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("arbitrary-name", 7,
                new int[] {0, 1, 3, 4}, new int[] {1, 2, 4, 5},
                new boolean[] {true, true, true, true});
        BondPercolationSimulator simulator = new BondPercolationSimulator(graph);
        assertStatistics(simulator.run(1, new Random(42)), 3, 3, 7.0 / 3, 3, 3, 7.0 / 3);
        assertStatistics(simulator.run(0, new Random(42)), 1, 1, 1, 1, 1, 1);
    }

    @Test
    void seededUndirectedSamplingMatchesDirectedReferenceAtIntermediateProbabilities() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("undirected", 8,
                new int[] {0, 1, 2, 2, 4, 5, 6}, new int[] {1, 2, 0, 3, 5, 6, 4},
                new boolean[] {true, true, true, true, true, true, true});
        BondPercolationSimulator simulator = new BondPercolationSimulator(graph);
        for (long seed = 0; seed < 20; seed++) {
            for (double p : new double[] {0.2, 0.5, 0.8}) {
                Random sampling = new Random(seed);
                var sources = new ArrayList<Integer>();
                var targets = new ArrayList<Integer>();
                for (int u = 0; u < graph.n; u++) {
                    var range = graph.outNeighborRange(u);
                    for (int e = range.start; e < range.end; e++) {
                        int v = graph.getOutNeighbor(e);
                        if (u < v && sampling.nextDouble() < p) {
                            sources.add(u);
                            targets.add(v);
                            sources.add(v);
                            targets.add(u);
                        }
                    }
                }
                // Mark reciprocal arcs as directed, forcing the full SCC algorithm.
                DirectedGraph reference = DirectedGraph.fromEdgeListWithUndirectedFlag("reference", graph.n,
                        sources.stream().mapToInt(Integer::intValue).toArray(),
                        targets.stream().mapToInt(Integer::intValue).toArray(), new boolean[sources.size()]);
                BondPercolationResult expected = new BondPercolationSimulator(reference).run(1, new Random(1));
                assertStatistics(simulator.run(p, new Random(seed)),
                        expected.gwccSize, expected.gwccSecondSize, expected.gwccMeanSize,
                        expected.gsccSize, expected.gsccSecondSize, expected.gsccMeanSize);
            }
        }
    }

    @Test
    void mixedGraphKeepsDistinctWeakAndStrongComponentsRegardlessOfName() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("ER", 5,
                new int[] {0, 1, 2}, new int[] {1, 2, 3}, new boolean[] {true, false, true});
        assertStatistics(new BondPercolationSimulator(graph).run(1, new Random(42)),
                4, 1, 2.5, 2, 2, 5.0 / 3);
    }

    @Test
    void reciprocalDirectedEdgesAreSampledIndependently() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("ER", 2,
                new int[] {0, 1}, new int[] {1, 0}, new boolean[] {false, false});
        Random sampling = new Random(42) {
            private int calls;

            @Override
            public double nextDouble() {
                return calls++ == 0 ? 0.25 : 0.75;
            }
        };
        assertStatistics(new BondPercolationSimulator(graph).run(0.5, sampling), 2, 0, 2, 1, 1, 1);
    }

    @Test
    void emptyGraphHasZeroComponentStatistics() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("empty", 0,
                new int[0], new int[0], new boolean[0]);
        assertStatistics(new BondPercolationSimulator(graph).run(0.5, new Random(42)), 0, 0, 0, 0, 0, 0);
    }

    private static void assertStatistics(BondPercolationResult result,
            int weak, int weakSecond, double weakMean, int strong, int strongSecond, double strongMean) {
        assertEquals(weak, result.gwccSize);
        assertEquals(weakSecond, result.gwccSecondSize);
        assertEquals(weakMean, result.gwccMeanSize, 1e-12);
        assertEquals(strong, result.gsccSize);
        assertEquals(strongSecond, result.gsccSecondSize);
        assertEquals(strongMean, result.gsccMeanSize, 1e-12);
    }
}
