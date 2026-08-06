package sim.network;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

class DirectedGraphTest {

    @Test
    void increasesReciprocityWhilePreservingDegrees() {
        DirectedGraph graph = graph(4,
                new int[] { 0, 1, 2, 3 },
                new int[] { 1, 2, 3, 0 });
        int[] inDegrees = inDegrees(graph);
        int[] outDegrees = outDegrees(graph);

        DirectedGraph rewired = graph.increaseReciprocity(1.0, 10_000, 100, 42L);

        assertEquals(1.0, rewired.reciprocity());
        assertEquals(graph.n, rewired.n);
        assertEquals(graph.m, rewired.m);
        assertArrayEquals(inDegrees, inDegrees(rewired));
        assertArrayEquals(outDegrees, outDegrees(rewired));
        assertNoSelfLoopsOrParallelEdges(rewired);
    }

    @Test
    void fixedSeedProducesSameGraph() {
        DirectedGraph graph = graph(6,
                new int[] { 0, 1, 2, 3, 4, 5 },
                new int[] { 1, 2, 3, 4, 5, 0 });

        DirectedGraph first = graph.increaseReciprocity(2.0 / 3.0, 10_000, 200, 7L);
        DirectedGraph second = graph.increaseReciprocity(2.0 / 3.0, 10_000, 200, 7L);

        assertEquals(edgeSet(first), edgeSet(second));
    }

    @Test
    void neutralSwapsPreserveReciprocalArcCount() {
        DirectedGraph graph = graph(6,
                new int[] { 0, 1, 2, 3, 4, 5 },
                new int[] { 1, 0, 3, 2, 5, 4 });

        DirectedGraph rewired = graph.increaseReciprocity(1.0, 0, 1_000, 123L);

        assertEquals(graph.countReciprocalArcs(), rewired.countReciprocalArcs());
        assertArrayEquals(inDegrees(graph), inDegrees(rewired));
        assertArrayEquals(outDegrees(graph), outDegrees(rewired));
    }

    @Test
    void countsReciprocalParallelArcsByMatchedMultiplicityAndIgnoresLoops() {
        DirectedGraph graph = graph(3,
                new int[] { 0, 0, 0, 1, 2 },
                new int[] { 1, 1, 1, 0, 2 });

        assertEquals(2, graph.countReciprocalArcs());
        assertEquals(0.4, graph.reciprocity());

        DirectedGraph unchanged = graph.increaseReciprocity(0.4, 0, 0, 1L);
        assertEquals(2, unchanged.countReciprocalArcs());
    }

    @Test
    void rejectsUndirectedEdgesAndInvalidArguments() {
        DirectedGraph undirected = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "undirected", 2, new int[] { 0 }, new int[] { 1 }, new boolean[] { true });
        DirectedGraph directed = graph(2, new int[] { 0 }, new int[] { 1 });

        assertThrows(IllegalArgumentException.class,
                () -> undirected.increaseReciprocity(0.5, 10, 10, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> directed.increaseReciprocity(-0.1, 10, 10, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> directed.increaseReciprocity(Double.NaN, 10, 10, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> directed.increaseReciprocity(0.5, -1, 10, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> directed.increaseReciprocity(0.5, 10, -1, 1L));
    }

    @Test
    void throwsWhenTargetCannotBeReachedWithinAttemptLimit() {
        DirectedGraph graph = graph(2, new int[] { 0 }, new int[] { 1 });

        IllegalStateException exception = assertThrows(IllegalStateException.class,
                () -> graph.increaseReciprocity(1.0, 100, 0, 1L));

        assertTrue(exception.getMessage().contains("achieved 0.0"));
        assertEquals(0.0, graph.reciprocity());
    }

    @Test
    void emptyGraphHasZeroReciprocity() {
        DirectedGraph graph = graph(3, new int[0], new int[0]);

        assertEquals(0, graph.countReciprocalArcs());
        assertEquals(0.0, graph.reciprocity());
        assertEquals(0.0, graph.increaseReciprocity(0.0, 0, 10, 1L).reciprocity());
    }

    private static DirectedGraph graph(int n, int[] sources, int[] destinations) {
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", n, sources, destinations, new boolean[sources.length]);
    }

    private static int[] inDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int u = 0; u < graph.n; u++) {
            DirectedGraph.IntRange range = graph.inNeighborRange(u);
            degrees[u] = range.end - range.start;
        }
        return degrees;
    }

    private static int[] outDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int u = 0; u < graph.n; u++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(u);
            degrees[u] = range.end - range.start;
        }
        return degrees;
    }

    private static Set<Long> edgeSet(DirectedGraph graph) {
        Set<Long> edges = new HashSet<>();
        for (int u = 0; u < graph.n; u++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(u);
            for (int i = range.start; i < range.end; i++) {
                edges.add(((long) u << 32) | (graph.getOutNeighbor(i) & 0xffffffffL));
            }
        }
        return edges;
    }

    private static void assertNoSelfLoopsOrParallelEdges(DirectedGraph graph) {
        Set<Long> edges = new HashSet<>();
        for (int u = 0; u < graph.n; u++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(u);
            for (int i = range.start; i < range.end; i++) {
                int v = graph.getOutNeighbor(i);
                assertFalse(u == v, "self-loop: " + u + " -> " + v);
                long key = ((long) u << 32) | (v & 0xffffffffL);
                assertTrue(edges.add(key), "parallel edge: " + u + " -> " + v);
            }
        }
    }
}
