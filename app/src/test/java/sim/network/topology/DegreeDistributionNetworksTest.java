package sim.network.topology;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.IntRange;
import sim.network.topology.DegreeDistributionNetworks.Distribution;
import sim.network.topology.DegreeDistributionNetworks.Variant;

class DegreeDistributionNetworksTest {

    @Test
    void generatesFourCoupledDegreeDistributionVariants() {
        int kdMin = 2;
        int kdMax = 40;
        List<Variant> variants = DegreeDistributionNetworks.generate(
                "comparison", 5_000, kdMin, kdMax, 2.5, 42L);

        assertEquals(4, variants.size());
        assertVariant(variants.get(0), Distribution.POW, Distribution.POW);
        assertVariant(variants.get(1), Distribution.POW, Distribution.POI);
        assertVariant(variants.get(2), Distribution.POI, Distribution.POW);
        assertVariant(variants.get(3), Distribution.POI, Distribution.POI);

        DirectedGraph powPow = variants.get(0).graph();
        DirectedGraph powPoi = variants.get(1).graph();
        DirectedGraph poiPow = variants.get(2).graph();
        DirectedGraph poiPoi = variants.get(3).graph();

        int[] k1 = inDegrees(powPow);
        int[] k2 = outDegrees(powPow);
        int[] p1 = outDegrees(powPoi);

        assertArrayEquals(k1, inDegrees(powPoi), "Pow/Pow and Pow/Poi must share k1");
        assertArrayEquals(p1, inDegrees(poiPow), "reversed graph must use p1 as in-degrees");
        assertArrayEquals(k1, outDegrees(poiPow), "reversed graph must use k1 as out-degrees");
        assertArrayEquals(p1, inDegrees(poiPoi), "Poi/Poi must reuse p1 as in-degrees");

        for (int degree : k1) {
            assertTrue(degree >= kdMin && degree <= kdMax);
        }
        for (int degree : k2) {
            assertTrue(degree >= kdMin && degree <= kdMax);
        }

        for (Variant variant : variants) {
            assertEquals(sum(inDegrees(variant.graph())), sum(outDegrees(variant.graph())));
            assertEquals(powPow.m, variant.graph().m);
        }

        assertIsExactReverse(powPoi, poiPow);
    }

    @Test
    void fixedSeedIsReproducible() {
        List<Variant> first = DegreeDistributionNetworks.generate("comparison", 1_000, 1, 20, 2.8, 7L);
        List<Variant> second = DegreeDistributionNetworks.generate("comparison", 1_000, 1, 20, 2.8, 7L);

        for (int i = 0; i < first.size(); i++) {
            assertArrayEquals(inDegrees(first.get(i).graph()), inDegrees(second.get(i).graph()));
            assertArrayEquals(outDegrees(first.get(i).graph()), outDegrees(second.get(i).graph()));
            assertSameOutgoingEdges(first.get(i).graph(), second.get(i).graph());
        }
    }

    @Test
    void invalidArgumentsThrow() {
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate(null, 100, 1, 10, 2.5, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate("test", 0, 1, 10, 2.5, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate("test", 100, 0, 10, 2.5, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate("test", 100, 10, 9, 2.5, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate("test", 100, 1, 10, 1.0, 1L));
        assertThrows(IllegalArgumentException.class,
                () -> DegreeDistributionNetworks.generate("test", 100, 1, 10, Double.NaN, 1L));
    }

    private static void assertVariant(Variant variant, Distribution in, Distribution out) {
        assertEquals(in, variant.inDistribution());
        assertEquals(out, variant.outDistribution());
    }

    private static int[] inDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int u = 0; u < graph.n; u++) {
            IntRange range = graph.inNeighborRange(u);
            degrees[u] = range.end - range.start;
        }
        return degrees;
    }

    private static int[] outDegrees(DirectedGraph graph) {
        int[] degrees = new int[graph.n];
        for (int u = 0; u < graph.n; u++) {
            IntRange range = graph.outNeighborRange(u);
            degrees[u] = range.end - range.start;
        }
        return degrees;
    }

    private static long sum(int[] values) {
        long total = 0L;
        for (int value : values) {
            total += value;
        }
        return total;
    }

    private static void assertIsExactReverse(DirectedGraph graph, DirectedGraph reversed) {
        for (int u = 0; u < graph.n; u++) {
            IntRange out = graph.outNeighborRange(u);
            IntRange reversedIn = reversed.inNeighborRange(u);
            assertEquals(out.end - out.start, reversedIn.end - reversedIn.start);
            for (int offset = 0; offset < out.end - out.start; offset++) {
                assertEquals(graph.getOutNeighbor(out.start + offset),
                        reversed.getInNeighbor(reversedIn.start + offset));
            }
        }
    }

    private static void assertSameOutgoingEdges(DirectedGraph first, DirectedGraph second) {
        assertEquals(first.n, second.n);
        assertEquals(first.m, second.m);
        for (int u = 0; u < first.n; u++) {
            IntRange firstRange = first.outNeighborRange(u);
            IntRange secondRange = second.outNeighborRange(u);
            assertEquals(firstRange.end - firstRange.start, secondRange.end - secondRange.start);
            for (int offset = 0; offset < firstRange.end - firstRange.start; offset++) {
                assertEquals(first.getOutNeighbor(firstRange.start + offset),
                        second.getOutNeighbor(secondRange.start + offset));
            }
        }
    }
}
