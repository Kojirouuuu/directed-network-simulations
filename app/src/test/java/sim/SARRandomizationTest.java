package sim;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Path;
import java.util.Arrays;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;

class SARRandomizationTest {

    @Test
    void recordsEachRandomizationModeInItsOutputPath() {
        Path networkPath = Path.of("DirectedER/N=1000/kAve=6.00");

        assertEquals(networkPath.resolve("randomization=none"),
                SAR.appendRandomizationPath(networkPath, SAR.RandomizationMode.NONE));
        assertEquals(networkPath.resolve("randomization=edge-swap"),
                SAR.appendRandomizationPath(networkPath, SAR.RandomizationMode.EDGE_SWAP));
        assertEquals(networkPath.resolve("randomization=in-degree-shuffle"),
                SAR.appendRandomizationPath(networkPath, SAR.RandomizationMode.SHUFFLE_IN_DEGREES));
        assertEquals(networkPath.resolve("randomization=out-degree-shuffle"),
                SAR.appendRandomizationPath(networkPath, SAR.RandomizationMode.SHUFFLE_OUT_DEGREES));
    }

    @Test
    void degreeShuffleModesUseEdgeSwappedInputFiles() {
        assertTrue(SAR.RandomizationMode.EDGE_SWAP.usesEdgeSwappedInput());
        assertTrue(SAR.RandomizationMode.SHUFFLE_IN_DEGREES.usesEdgeSwappedInput());
        assertTrue(SAR.RandomizationMode.SHUFFLE_OUT_DEGREES.usesEdgeSwappedInput());
        assertFalse(SAR.RandomizationMode.NONE.usesEdgeSwappedInput());
    }

    @Test
    void dispatchesLoadedGraphsAccordingToMode() {
        DirectedGraph graph = graphWithDistinctDegreeSequences();

        assertSame(graph, SAR.applyRandomization(graph, SAR.RandomizationMode.NONE, 7L, true));
        assertSame(graph, SAR.applyRandomization(graph, SAR.RandomizationMode.EDGE_SWAP, 7L, true));

        DirectedGraph inShuffled = SAR.applyRandomization(
                graph, SAR.RandomizationMode.SHUFFLE_IN_DEGREES, 7L, true);
        DirectedGraph outShuffled = SAR.applyRandomization(
                graph, SAR.RandomizationMode.SHUFFLE_OUT_DEGREES, 7L, true);

        assertArrayEquals(outDegrees(graph), outDegrees(inShuffled));
        assertArrayEquals(sorted(inDegrees(graph)), sorted(inDegrees(inShuffled)));
        assertArrayEquals(inDegrees(graph), inDegrees(outShuffled));
        assertArrayEquals(sorted(outDegrees(graph)), sorted(outDegrees(outShuffled)));
    }

    @Test
    void appliesEdgeSwapsToGeneratedGraphs() {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 8,
                new int[] { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7 },
                new int[] { 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0, 0, 1 },
                new boolean[16]);

        DirectedGraph randomized = SAR.applyRandomization(
                graph, SAR.RandomizationMode.EDGE_SWAP, 0L, false);

        assertEquals("test_randomized", randomized.name);
        assertArrayEquals(inDegrees(graph), inDegrees(randomized));
        assertArrayEquals(outDegrees(graph), outDegrees(randomized));
    }

    private static DirectedGraph graphWithDistinctDegreeSequences() {
        int[] sources = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 3 };
        int[] destinations = { 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 };
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 5, sources, destinations, new boolean[sources.length]);
    }

    private static int[] inDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int vertex = 0; vertex < graph.n; vertex++) {
            DirectedGraph.IntRange range = graph.inNeighborRange(vertex);
            degrees[vertex] = range.end - range.start;
        }
        return degrees;
    }

    private static int[] outDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int vertex = 0; vertex < graph.n; vertex++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(vertex);
            degrees[vertex] = range.end - range.start;
        }
        return degrees;
    }

    private static int[] sorted(int[] values) {
        int[] copy = Arrays.copyOf(values, values.length);
        Arrays.sort(copy);
        return copy;
    }
}
