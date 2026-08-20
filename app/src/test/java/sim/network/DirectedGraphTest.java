package sim.network;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.junit.jupiter.api.Test;

class DirectedGraphTest {

    @Test
    void countsFeedForwardLoopsAsOrderedNonInducedTriples() {
        DirectedGraph completeDirectedTriangle = graph(3,
                new int[] { 0, 1, 1, 2, 0, 2 },
                new int[] { 1, 0, 2, 1, 2, 0 });
        DirectedGraph parallelEdgesAndLoop = graph(3,
                new int[] { 0, 0, 0, 1, 1, 0 },
                new int[] { 0, 1, 1, 2, 2, 2 });

        assertEquals(6L, completeDirectedTriangle.countFeedForwardLoops());
        assertEquals(1L, parallelEdgesAndLoop.countFeedForwardLoops());
    }

    @Test
    void feedForwardLoopCountMatchesBruteForceOnSmallMultigraphs() {
        for (long seed = 0L; seed < 20L; seed++) {
            Random random = new Random(seed);
            int[] sources = new int[24];
            int[] destinations = new int[24];
            for (int edge = 0; edge < sources.length; edge++) {
                sources[edge] = random.nextInt(6);
                destinations[edge] = random.nextInt(6);
            }
            DirectedGraph graph = graph(6, sources, destinations);

            assertEquals(bruteForceFeedForwardLoopCount(graph), graph.countFeedForwardLoops(),
                    "seed=" + seed);
        }
    }

    @Test
    void increasesFeedForwardLoopsWhilePreservingDegrees() {
        DirectedGraph graph = graph(5,
                new int[] { 0, 1, 0, 4 },
                new int[] { 1, 2, 3, 2 });
        int[] inDegrees = inDegrees(graph);
        int[] outDegrees = outDegrees(graph);
        Set<Long> originalEdges = edgeSet(graph);

        DirectedGraph rewired = graph.increaseFeedForwardLoops(1L, 1_000L, 42L);

        assertEquals(0L, graph.countFeedForwardLoops());
        assertEquals(1L, rewired.countFeedForwardLoops());
        assertEquals(rewired.countFeedForwardLoops(), copyGraph(rewired).countFeedForwardLoops());
        assertEquals(originalEdges, edgeSet(graph));
        assertEquals(graph.n, rewired.n);
        assertEquals(graph.m, rewired.m);
        assertArrayEquals(inDegrees, inDegrees(rewired));
        assertArrayEquals(outDegrees, outDegrees(rewired));
        assertNoSelfLoopsOrParallelEdges(rewired);
        assertTrue(edgeSet(rewired).contains(edgeKey(0, 2)));
        assertTrue(edgeSet(rewired).contains(edgeKey(4, 3)));
    }

    @Test
    void feedForwardLoopRewiringLeavesExistingLoopsAndParallelEdgesUntouched() {
        DirectedGraph graph = graph(7,
                new int[] { 0, 1, 0, 4, 5, 5, 5 },
                new int[] { 1, 2, 3, 2, 5, 6, 6 });

        DirectedGraph rewired = graph.increaseFeedForwardLoops(1L, 10_000L, 42L);

        assertEquals(1, edgeMultiplicity(rewired, 5, 5));
        assertEquals(2, edgeMultiplicity(rewired, 5, 6));
        assertEquals(rewired.countFeedForwardLoops(), copyGraph(rewired).countFeedForwardLoops());
    }

    @Test
    void feedForwardLoopRewiringIsDeterministicForFixedSeed() {
        DirectedGraph graph = graph(5,
                new int[] { 0, 1, 0, 4 },
                new int[] { 1, 2, 3, 2 });

        DirectedGraph first = graph.increaseFeedForwardLoops(1L, 1_000L, 7L);
        DirectedGraph second = graph.increaseFeedForwardLoops(1L, 1_000L, 7L);

        assertEquals(edgeSet(first), edgeSet(second));
        assertEquals(first.countFeedForwardLoops(), second.countFeedForwardLoops());
    }

    @Test
    void repeatedFeedForwardLoopSwapsKeepIncrementalCountExact() {
        DirectedGraph graph = graph(15,
                new int[] { 0, 1, 0, 4, 5, 6, 5, 9, 10, 11, 10, 14 },
                new int[] { 1, 2, 3, 2, 6, 7, 8, 7, 11, 12, 13, 12 });

        DirectedGraph rewired = graph.increaseFeedForwardLoops(3L, 100_000L, 99L);

        assertEquals(3L, rewired.countFeedForwardLoops());
        assertEquals(bruteForceFeedForwardLoopCount(rewired), rewired.countFeedForwardLoops());
        assertArrayEquals(inDegrees(graph), inDegrees(rewired));
        assertArrayEquals(outDegrees(graph), outDegrees(rewired));
    }

    @Test
    void feedForwardLoopRewiringValidatesArgumentsAndReportsUnreachableTarget() {
        DirectedGraph graph = graph(3,
                new int[] { 0, 1 },
                new int[] { 1, 2 });
        DirectedGraph undirected = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "undirected", 2, new int[] { 0 }, new int[] { 1 }, new boolean[] { true });

        assertThrows(IllegalArgumentException.class,
                () -> graph.increaseFeedForwardLoops(-1L, 10L, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> graph.increaseFeedForwardLoops(1L, -1L, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> undirected.increaseFeedForwardLoops(1L, 10L, 1L));

        IllegalStateException exception = assertThrows(IllegalStateException.class,
                () -> graph.increaseFeedForwardLoops(1L, 100L, 1L));
        assertTrue(exception.getMessage().contains("target NFFL 1"));
        assertTrue(exception.getMessage().contains("achieved 0"));
    }

    @Test
    void feedForwardLoopTargetAlreadyReachedNeedsNoAttempts() {
        DirectedGraph graph = graph(3,
                new int[] { 0, 1, 0 },
                new int[] { 1, 2, 2 });

        DirectedGraph unchanged = graph.increaseFeedForwardLoops(1L, 0L, 1L);

        assertEquals(edgeSet(graph), edgeSet(unchanged));
        assertEquals(1L, unchanged.countFeedForwardLoops());
    }

    @Test
    void randomizesByEdgeSwapsWhilePreservingDegrees() {
        DirectedGraph graph = graph(8,
                new int[] { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7 },
                new int[] { 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0, 0, 1 });
        Set<Long> originalEdges = edgeSet(graph);
        int[] originalInDegrees = inDegrees(graph);
        int[] originalOutDegrees = outDegrees(graph);

        DirectedGraph randomized = graph.randomizeByEdgeSwaps(0L);

        assertEquals("test_randomized", randomized.name);
        assertEquals(graph.n, randomized.n);
        assertEquals(graph.m, randomized.m);
        assertArrayEquals(originalInDegrees, inDegrees(randomized));
        assertArrayEquals(originalOutDegrees, outDegrees(randomized));
        assertEquals(originalEdges, edgeSet(graph));
        assertFalse(originalEdges.equals(edgeSet(randomized)));
        assertTrue(randomized.reciprocity() > graph.reciprocity());
        assertNoSelfLoopsOrParallelEdges(randomized);
    }

    @Test
    void randomEdgeSwapsAreDeterministicForFixedSeed() {
        DirectedGraph graph = graph(8,
                new int[] { 0, 1, 2, 3, 4, 5, 6, 7 },
                new int[] { 1, 2, 3, 4, 5, 6, 7, 0 });

        DirectedGraph first = graph.randomizeByEdgeSwaps(7L);
        DirectedGraph second = graph.randomizeByEdgeSwaps(7L);

        assertEquals(edgeSet(first), edgeSet(second));
    }

    @Test
    void randomEdgeSwapsRejectUndirectedEdgesAndReportUnreachableTarget() {
        DirectedGraph undirected = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "undirected", 2, new int[] { 0 }, new int[] { 1 }, new boolean[] { true });
        DirectedGraph complete = graph(3,
                new int[] { 0, 0, 1, 1, 2, 2 },
                new int[] { 1, 2, 0, 2, 0, 1 });

        assertThrows(IllegalArgumentException.class, () -> undirected.randomizeByEdgeSwaps(1L));

        IllegalStateException exception = assertThrows(
                IllegalStateException.class, () -> complete.randomizeByEdgeSwaps(1L));
        assertTrue(exception.getMessage().contains("60 random edge swaps"));
        assertTrue(exception.getMessage().contains("600 attempts"));
        assertTrue(exception.getMessage().contains("accepted 0"));
    }

    @Test
    void randomEdgeSwapsReturnNamedCopyForEmptyGraph() {
        DirectedGraph empty = graph(3, new int[0], new int[0]);

        DirectedGraph randomized = empty.randomizeByEdgeSwaps(1L);

        assertEquals("test_randomized", randomized.name);
        assertEquals(0, randomized.m);
        assertArrayEquals(inDegrees(empty), inDegrees(randomized));
        assertArrayEquals(outDegrees(empty), outDegrees(randomized));
    }

    @Test
    void shufflesInDegreeSequenceWhilePreservingBothDegreeDistributions() {
        DirectedGraph graph = graphWithDistinctDegreeSequences();
        int[] originalInDegrees = inDegrees(graph);
        int[] originalOutDegrees = outDegrees(graph);

        DirectedGraph randomized = graph.randomizeByShuffledDegreeSequence(
                DirectedGraph.DegreeSide.IN, 7L);

        assertEquals("test_in_degree_shuffled", randomized.name);
        assertEquals(graph.n, randomized.n);
        assertEquals(graph.m, randomized.m);
        assertArrayEquals(originalOutDegrees, outDegrees(randomized));
        assertArrayEquals(sorted(originalInDegrees), sorted(inDegrees(randomized)));
        assertFalse(Arrays.equals(originalInDegrees, inDegrees(randomized)));
        assertArrayEquals(originalInDegrees, inDegrees(graph));
        assertArrayEquals(originalOutDegrees, outDegrees(graph));
    }

    @Test
    void shufflesOutDegreeSequenceWhilePreservingBothDegreeDistributions() {
        DirectedGraph graph = graphWithDistinctDegreeSequences();
        int[] originalInDegrees = inDegrees(graph);
        int[] originalOutDegrees = outDegrees(graph);

        DirectedGraph randomized = graph.randomizeByShuffledDegreeSequence(
                DirectedGraph.DegreeSide.OUT, 7L);

        assertEquals("test_out_degree_shuffled", randomized.name);
        assertEquals(graph.n, randomized.n);
        assertEquals(graph.m, randomized.m);
        assertArrayEquals(originalInDegrees, inDegrees(randomized));
        assertArrayEquals(sorted(originalOutDegrees), sorted(outDegrees(randomized)));
        assertFalse(Arrays.equals(originalOutDegrees, outDegrees(randomized)));
        assertArrayEquals(originalInDegrees, inDegrees(graph));
        assertArrayEquals(originalOutDegrees, outDegrees(graph));
    }

    @Test
    void shuffledDegreeRandomizationIsDeterministicAndAllowsConfigurationModelEdges() {
        DirectedGraph graph = graphWithDistinctDegreeSequences();

        DirectedGraph first = graph.randomizeByShuffledDegreeSequence(
                DirectedGraph.DegreeSide.IN, 123L);
        DirectedGraph second = graph.randomizeByShuffledDegreeSequence(
                DirectedGraph.DegreeSide.IN, 123L);
        DirectedGraph loops = graph(1, new int[] { 0, 0 }, new int[] { 0, 0 })
                .randomizeByShuffledDegreeSequence(DirectedGraph.DegreeSide.OUT, 1L);

        assertEquals(edgeMultiset(first), edgeMultiset(second));
        assertEquals(2, edgeMultiplicity(loops, 0, 0));
    }

    @Test
    void shuffledDegreeRandomizationHandlesEmptyGraphAndRejectsInvalidInputs() {
        DirectedGraph empty = graph(3, new int[0], new int[0]);
        DirectedGraph randomized = empty.randomizeByShuffledDegreeSequence(
                DirectedGraph.DegreeSide.IN, 1L);
        DirectedGraph undirected = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "undirected", 2, new int[] { 0 }, new int[] { 1 }, new boolean[] { true });

        assertEquals("test_in_degree_shuffled", randomized.name);
        assertEquals(0, randomized.m);
        assertThrows(IllegalArgumentException.class,
                () -> empty.randomizeByShuffledDegreeSequence(null, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> undirected.randomizeByShuffledDegreeSequence(DirectedGraph.DegreeSide.IN, 1L));
    }

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
    void acceptsIncreaseAttemptBudgetLargerThanIntegerMaxValue() {
        DirectedGraph graph = graph(2, new int[] { 0, 1 }, new int[] { 1, 0 });

        DirectedGraph unchanged = graph.increaseReciprocity(
                1.0, (long) Integer.MAX_VALUE + 1L, 0L, 1L);

        assertEquals(1.0, unchanged.reciprocity());
    }

    @Test
    void rejectsBudgetBelowTheoreticalMinimumBeforeTryingSwaps() {
        DirectedGraph graph = graph(4,
                new int[] { 0, 1, 2, 3 },
                new int[] { 1, 2, 3, 0 });

        IllegalStateException exception = assertThrows(IllegalStateException.class,
                () -> graph.increaseReciprocity(1.0, 0L, 0L, 1L));

        assertTrue(exception.getMessage().contains("theoretically required"));
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

    private static DirectedGraph graphWithDistinctDegreeSequences() {
        return graph(5,
                new int[] { 0, 0, 0, 0, 1, 1, 1, 2, 2, 3 },
                new int[] { 1, 2, 2, 3, 3, 3, 4, 4, 4, 4 });
    }

    private static int[] sorted(int[] values) {
        int[] copy = Arrays.copyOf(values, values.length);
        Arrays.sort(copy);
        return copy;
    }

    private static List<Long> edgeMultiset(DirectedGraph graph) {
        List<Long> edges = new ArrayList<>(graph.m);
        for (int source = 0; source < graph.n; source++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(source);
            for (int edge = range.start; edge < range.end; edge++) {
                edges.add(edgeKey(source, graph.getOutNeighbor(edge)));
            }
        }
        edges.sort(Long::compareTo);
        return edges;
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

    private static long edgeKey(int source, int destination) {
        return ((long) source << 32) | (destination & 0xffffffffL);
    }

    private static int edgeMultiplicity(DirectedGraph graph, int source, int destination) {
        int count = 0;
        DirectedGraph.IntRange range = graph.outNeighborRange(source);
        for (int edge = range.start; edge < range.end; edge++) {
            if (graph.getOutNeighbor(edge) == destination) {
                count++;
            }
        }
        return count;
    }

    private static DirectedGraph copyGraph(DirectedGraph graph) {
        int[] sources = new int[graph.m];
        int[] destinations = new int[graph.m];
        int edge = 0;
        for (int source = 0; source < graph.n; source++) {
            DirectedGraph.IntRange range = graph.outNeighborRange(source);
            for (int index = range.start; index < range.end; index++) {
                sources[edge] = source;
                destinations[edge] = graph.getOutNeighbor(index);
                edge++;
            }
        }
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                graph.name + "_copy", graph.n, sources, destinations, new boolean[graph.m]);
    }

    private static long bruteForceFeedForwardLoopCount(DirectedGraph graph) {
        Set<Long> edges = edgeSet(graph);
        long count = 0L;
        for (int source = 0; source < graph.n; source++) {
            for (int middle = 0; middle < graph.n; middle++) {
                if (middle == source || !edges.contains(edgeKey(source, middle))) {
                    continue;
                }
                for (int destination = 0; destination < graph.n; destination++) {
                    if (destination != source && destination != middle
                            && edges.contains(edgeKey(middle, destination))
                            && edges.contains(edgeKey(source, destination))) {
                        count++;
                    }
                }
            }
        }
        return count;
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
