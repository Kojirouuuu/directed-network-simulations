package sim.network.topology;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class TwoPointDirectedCMTest {
    @Test
    void countsRepresentCorrelationEndpointsAndHalfUpRounding() {
        assertEquals(new TwoPointDirectedCM.DegreePairCounts(375, 125, 125, 375),
                TwoPointDirectedCM.degreePairCounts(1000, 0.5, 0.5));
        assertEquals(new TwoPointDirectedCM.DegreePairCounts(20, 0, 0, 80),
                TwoPointDirectedCM.degreePairCounts(100, 0.2, 1));
        assertEquals(new TwoPointDirectedCM.DegreePairCounts(0, 20, 20, 60),
                TwoPointDirectedCM.degreePairCounts(100, 0.2, -0.25));
        var rounded = TwoPointDirectedCM.degreePairCounts(5, 0.5, -1);
        assertEquals(new TwoPointDirectedCM.DegreePairCounts(1, 2, 2, 0), rounded);
        assertEquals(0.6, rounded.realizedP(), 1e-12);
        assertEquals(-2.0 / 3, rounded.realizedR(), 1e-12);
    }

    @Test
    void generatedGraphPreservesJointDegreesAcrossSeeds() {
        for (long seed : new long[] {42, 43}) {
            var graph = TwoPointDirectedCM.generate("test", 1000, 2, 8, 0.5, 0.5, seed);
            int[] counts = new int[4];
            for (int vertex = 0; vertex < graph.n; vertex++) {
                var in = graph.inNeighborRange(vertex);
                var out = graph.outNeighborRange(vertex);
                int ki = in.end - in.start;
                int ko = out.end - out.start;
                assertTrue(ki == 2 || ki == 8);
                assertTrue(ko == 2 || ko == 8);
                counts[(ki == 2 ? 0 : 2) + (ko == 2 ? 0 : 1)]++;
            }
            assertArrayEquals(new int[] {375, 125, 125, 375}, counts);
            assertEquals(5000, graph.m);
            for (int edge = 0; edge < graph.m; edge++) {
                assertFalse(graph.isOutUndirected(edge));
            }
        }
    }

    @Test
    void seedReproducesConnectionsAndDifferentSeedChangesThem() {
        var first = TwoPointDirectedCM.generate("test", 100, 2, 8, 0.5, 0.5, 42);
        var same = TwoPointDirectedCM.generate("test", 100, 2, 8, 0.5, 0.5, 42);
        var other = TwoPointDirectedCM.generate("test", 100, 2, 8, 0.5, 0.5, 43);
        boolean changed = false;
        for (int edge = 0; edge < first.m; edge++) {
            assertEquals(first.getOutNeighbor(edge), same.getOutNeighbor(edge));
            assertEquals(first.getInNeighbor(edge), same.getInNeighbor(edge));
            changed |= first.getOutNeighbor(edge) != other.getOutNeighbor(edge);
        }
        assertTrue(changed);
    }

    @Test
    void rejectsInvalidParametersBeforeAllocatingGraph() {
        for (double p : new double[] {0, 1, Double.NaN, Double.POSITIVE_INFINITY}) {
            assertThrows(IllegalArgumentException.class,
                    () -> TwoPointDirectedCM.degreePairCounts(100, p, 0));
        }
        for (double r : new double[] {-0.251, 1.01, Double.NaN, Double.NEGATIVE_INFINITY}) {
            assertThrows(IllegalArgumentException.class,
                    () -> TwoPointDirectedCM.degreePairCounts(100, 0.2, r));
        }
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.degreePairCounts(1, 0.5, 0));
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.degreePairCounts(2, 0.01, 0));
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.generate("test", 10, 0, 8, 0.5, 0, 42));
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.generate("test", 10, 8, 8, 0.5, 0, 42));
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.generate(null, 10, 2, 8, 0.5, 0, 42));
        assertThrows(IllegalArgumentException.class,
                () -> TwoPointDirectedCM.generate("test", 2, 1, Integer.MAX_VALUE, 0.5, 0, 42));
    }
}
