package sim.network.topology;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.IntRange;

class DirectedERTest {

    @Test
    void probabilityZeroGeneratesEmptyGraph() {
        DirectedGraph graph = DirectedER.generateFromP(10, 0.0, 42L);

        assertEquals("DirectedER", graph.name);
        assertEquals(10, graph.n);
        assertEquals(0, graph.m);
    }

    @Test
    void probabilityOneGeneratesEveryNonLoopArc() {
        DirectedGraph graph = DirectedER.generateFromP(4, 1.0, 42L);

        assertEquals(12, graph.m);
        assertSimplePurelyDirectedGraph(graph);
        for (int vertex = 0; vertex < graph.n; vertex++) {
            assertEquals(3, degree(graph.outNeighborRange(vertex)));
            assertEquals(3, degree(graph.inNeighborRange(vertex)));
        }
    }

    @Test
    void fixedSeedIsReproducible() {
        DirectedGraph first = DirectedER.generateFromP(1_000, 0.01, 7L);
        DirectedGraph second = DirectedER.generateFromP(1_000, 0.01, 7L);

        assertEquals(edgeSet(first), edgeSet(second));
    }

    @Test
    void meanDegreeAndProbabilityApisAreEquivalent() {
        int n = 1_000;
        double meanDegree = 8.0;
        long seed = 123L;

        DirectedGraph fromMeanDegree = DirectedER.generateFromMeanDegree(n, meanDegree, seed);
        DirectedGraph fromProbability = DirectedER.generateFromP(n, meanDegree / (n - 1.0), seed);

        assertEquals(edgeSet(fromProbability), edgeSet(fromMeanDegree));
    }

    @Test
    void sparseGraphHasExpectedPoissonLikeDegreeStatistics() {
        int n = 20_000;
        double meanDegree = 4.0;
        DirectedGraph graph = DirectedER.generateFromMeanDegree(n, meanDegree, 2026L);

        assertSimplePurelyDirectedGraph(graph);
        assertEquals(meanDegree, graph.m / (double) n, 0.06);

        int zeroInDegree = 0;
        int zeroOutDegree = 0;
        for (int vertex = 0; vertex < n; vertex++) {
            if (degree(graph.inNeighborRange(vertex)) == 0) {
                zeroInDegree++;
            }
            if (degree(graph.outNeighborRange(vertex)) == 0) {
                zeroOutDegree++;
            }
        }
        double expectedZeroFraction = Math.pow(1.0 - meanDegree / (n - 1.0), n - 1.0);
        assertEquals(expectedZeroFraction, zeroInDegree / (double) n, 0.006);
        assertEquals(expectedZeroFraction, zeroOutDegree / (double) n, 0.006);
    }

    @Test
    void validatesArguments() {
        assertThrows(IllegalArgumentException.class, () -> DirectedER.generateFromP(0, 0.1, 1L));
        assertThrows(IllegalArgumentException.class, () -> DirectedER.generateFromP(10, -0.1, 1L));
        assertThrows(IllegalArgumentException.class, () -> DirectedER.generateFromP(10, 1.1, 1L));
        assertThrows(IllegalArgumentException.class, () -> DirectedER.generateFromP(10, Double.NaN, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DirectedER.generateFromP(10, Double.POSITIVE_INFINITY, 1L));

        assertThrows(IllegalArgumentException.class,
                () -> DirectedER.generateFromMeanDegree(10, -1.0, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DirectedER.generateFromMeanDegree(10, 10.0, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DirectedER.generateFromMeanDegree(10, Double.NaN, 1L));

        assertEquals(0, DirectedER.generateFromMeanDegree(1, 0.0, 1L).m);
        assertThrows(IllegalArgumentException.class,
                () -> DirectedER.generateFromMeanDegree(1, 0.1, 1L));
    }

    private static void assertSimplePurelyDirectedGraph(DirectedGraph graph) {
        Set<Long> edges = new HashSet<>();
        for (int source = 0; source < graph.n; source++) {
            IntRange range = graph.outNeighborRange(source);
            for (int index = range.start; index < range.end; index++) {
                int destination = graph.getOutNeighbor(index);
                assertTrue(source != destination, "self-loop: " + source + " -> " + destination);
                assertFalse(graph.isOutUndirected(index));
                assertTrue(edges.add(edgeKey(source, destination)),
                        "duplicate edge: " + source + " -> " + destination);
            }
        }
        assertEquals(graph.m, edges.size());
    }

    private static Set<Long> edgeSet(DirectedGraph graph) {
        Set<Long> edges = new HashSet<>();
        for (int source = 0; source < graph.n; source++) {
            IntRange range = graph.outNeighborRange(source);
            for (int index = range.start; index < range.end; index++) {
                edges.add(edgeKey(source, graph.getOutNeighbor(index)));
            }
        }
        return edges;
    }

    private static long edgeKey(int source, int destination) {
        return ((long) source << 32) | (destination & 0xffffffffL);
    }

    private static int degree(IntRange range) {
        return range.end - range.start;
    }
}
