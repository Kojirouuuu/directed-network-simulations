package sim.simulation;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;

class SARFixedIntervalRecordingTest {

    private static final double EPSILON = 1.0e-12;

    @Test
    void recordsFixedGridAndFillsTerminalStateThroughLastGridPoint() {
        DirectedGraph graph = graph(1, new int[0], new int[0]);

        for (boolean useGillespie : new boolean[] { false, true }) {
            SARResult result = simulate(useGillespie, graph, 0.0, 0.0,
                    0.35, 0.1, new int[] { 2 }, new int[] { 0 }, 42L);

            assertGrid(result, 0.0, 0.1, 0.2, 0.3);
            assertEquals(4, result.A.size());
            for (int i = 0; i < result.times.size(); i++) {
                assertEquals(1, result.A.get(i));
                assertEquals(0, result.R.get(i));
                assertEquals(0, result.Phi.get(i));
            }
        }
    }

    @Test
    void includesTMaxWhenItIsOnTheGridDespiteFloatingPointRounding() {
        DirectedGraph graph = graph(1, new int[0], new int[0]);

        for (boolean useGillespie : new boolean[] { false, true }) {
            SARResult result = simulate(useGillespie, graph, 0.0, 0.0,
                    0.3, 0.1, new int[] { 2 }, new int[] { 0 }, 7L);
            assertGrid(result, 0.0, 0.1, 0.2, 0.3);
        }
    }

    @Test
    void sampledAdoptedAndRecoveredCountsReflectRecoveryTime() {
        DirectedGraph graph = graph(1, new int[0], new int[0]);

        for (boolean useGillespie : new boolean[] { false, true }) {
            SARResult result = simulate(useGillespie, graph, 0.0, 1.0,
                    3.0, 0.25, new int[] { 2 }, new int[] { 0 }, 42L);
            double recoveryTime = result.tRecover[0];
            assertTrue(Double.isFinite(recoveryTime));

            for (int i = 0; i < result.times.size(); i++) {
                boolean recovered = result.times.get(i) >= recoveryTime;
                assertEquals(recovered ? 0 : 1, result.A.get(i));
                assertEquals(recovered ? 1 : 0, result.R.get(i));
            }
        }
    }

    @Test
    void capturesPhiOnlyTransmissionChanges() {
        DirectedGraph graph = graph(2, new int[] { 0 }, new int[] { 1 });

        for (boolean useGillespie : new boolean[] { false, true }) {
            SARResult result = simulate(useGillespie, graph, 10.0, 0.0,
                    2.0, 0.1, new int[] { 2, 2 }, new int[] { 0 }, 123L);

            assertEquals(0, result.Phi.get(0));
            assertEquals(1, result.Phi.get(result.Phi.size() - 1));
            assertTrue(result.Phi.stream().anyMatch(value -> value == 1));
            assertTrue(result.A.stream().allMatch(value -> value == 1));
            assertTrue(result.R.stream().allMatch(value -> value == 0));
        }
    }

    @Test
    void rejectsInvalidOrImpracticallySmallDt() {
        DirectedGraph graph = graph(1, new int[0], new int[0]);
        for (boolean useGillespie : new boolean[] { false, true }) {
            for (double dt : new double[] { 0.0, -0.1, Double.NaN, Double.POSITIVE_INFINITY,
                    Double.MIN_VALUE }) {
                assertThrows(IllegalArgumentException.class,
                        () -> simulate(useGillespie, graph, 0.0, 0.0,
                                1.0, dt, new int[] { 2 }, new int[] { 0 }, 1L));
            }
        }
    }

    @Test
    void existingApiStillRecordsEventTimes() {
        DirectedGraph graph = graph(1, new int[0], new int[0]);

        SARResult eventDriven = SARSimulator.simulate(
                graph, 0.0, 0.0, 0.0, 1.0, new int[] { 2 }, new int[] { 0 }, 1L);
        SARResult gillespie = SARGillespieSimulator.simulate(
                graph, 0.0, 0.0, 0.0, 1.0, new int[] { 2 }, new int[] { 0 }, 1L);

        assertGrid(eventDriven, 0.0);
        assertGrid(gillespie, 0.0);
    }

    private static SARResult simulate(boolean useGillespie, DirectedGraph graph,
            double lambdaDirected, double mu, double tMax, double dt,
            int[] thresholds, int[] initialInfecteds, long seed) {
        return useGillespie
                ? SARGillespieSimulator.simulate(
                        graph, lambdaDirected, 0.0, mu, tMax, dt, thresholds, initialInfecteds, seed)
                : SARSimulator.simulate(
                        graph, lambdaDirected, 0.0, mu, tMax, dt, thresholds, initialInfecteds, seed);
    }

    private static DirectedGraph graph(int n, int[] sources, int[] destinations) {
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", n, sources, destinations, new boolean[sources.length]);
    }

    private static void assertGrid(SARResult result, double... expectedTimes) {
        assertEquals(expectedTimes.length, result.times.size());
        for (int i = 0; i < expectedTimes.length; i++) {
            assertEquals(expectedTimes[i], result.times.get(i), EPSILON);
        }
    }
}
