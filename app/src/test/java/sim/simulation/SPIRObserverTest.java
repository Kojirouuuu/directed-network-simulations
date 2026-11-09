package sim.simulation;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.network.DirectedGraph;
import sim.simulation.SPIRGillespieSimulator.ModelMode;

class SPIRObserverTest extends SPIRTestSupport {
    @TempDir Path temp;

    private DirectedGraph graph(int n, boolean undirected, int... pairs) {
        int[] from = new int[pairs.length / 2], to = new int[from.length];
        boolean[] flags = new boolean[from.length];
        for (int k = 0; k < from.length; k++) {
            from[k] = pairs[2 * k];
            to[k] = pairs[2 * k + 1];
        }
        Arrays.fill(flags, undirected);
        return DirectedGraph.fromEdgeListWithUndirectedFlag("test", n, from, to, flags);
    }

    private SPIRResult simulate(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, int[] initial, int[] observers, long seed) {
        return simulateApi(
                graph, mode, lambda1, lambda2, mu1, mu2, 100, initial, observers, seed);
    }

    @Test
    void detectsIImmediatelyAndIsolatesObserverNeighborhoodButNotTriggeringI() {
        var graph = graph(4, true, 0, 1, 0, 2, 1, 2, 1, 3);
        var result = simulate(graph, ModelMode.SPIR, 0, 100, 1, 1, new int[] {1}, new int[] {0}, 42);
        assertEquals(result.tOnset[1], result.tQuarantine[0]);
        assertEquals(result.tOnset[1], result.tQuarantine[2]);
        assertFalse(result.quarantined[1]);
        assertFalse(result.quarantined[3]);
        assertTrue(Double.isNaN(result.tInfect[0]));
        assertTrue(Double.isNaN(result.tInfect[2]));
        assertTrue(result.tInfect[3] > result.tOnset[1]);
        assertEquals(2, result.quarantinedCount);
        for (int k = 0; k < result.times.size(); k++) {
            assertEquals(4, result.S.get(k) + result.P.get(k) + result.I.get(k) + result.R.get(k));
        }
    }

    @Test
    void presymptomaticObserverCanTransmitAndCannotDetectP() {
        var graph = graph(3, false, 0, 1, 1, 2);
        var result = simulate(graph, ModelMode.SPIR, 100, 0, 0, 0, new int[] {0}, new int[] {1}, 42);
        assertTrue(Double.isFinite(result.tInfect[2]));
        assertEquals(3, result.finalP);
        assertEquals(0, result.quarantinedCount);
        assertTrue(Double.isNaN(result.tOnset[1]));
    }

    @Test
    void observerSelfOnsetQuarantinesSelfAndNeighborsAndRecoveryContinues() {
        var graph = graph(3, true, 0, 1, 1, 2);
        int selfDetections = 0;
        for (int seed = 0; seed < 100; seed++) {
            var result = simulate(graph, ModelMode.SPIR, 10, 0, 1, 1, new int[] {0}, new int[] {1}, seed);
            if (Double.isFinite(result.tOnset[1]) && result.tOnset[1] < result.tOnset[0]
                    && (Double.isNaN(result.tOnset[2]) || result.tOnset[1] < result.tOnset[2])) {
                selfDetections++;
                for (int node = 0; node < 3; node++) {
                    assertEquals(result.tOnset[1], result.tQuarantine[node]);
                    if (Double.isFinite(result.tInfect[node])) {
                        assertTrue(result.tRecover[node] > result.tQuarantine[node]);
                    }
                }
                assertEquals(0, result.finalI);
                assertEquals(0, result.finalP);
            }
        }
        assertTrue(selfDetections > 0);
    }

    @Test
    void quarantinedPresymptomaticObserverStopsMonitoringButStillProgresses() {
        var graph = graph(5, true, 0, 1, 1, 2, 2, 3, 4, 2);
        int cases = 0;
        for (int seed = 0; seed < 300; seed++) {
            var result =
                    simulate(graph, ModelMode.SPIR, 5, 0, 1, 1, new int[] {0, 4}, new int[] {1, 2}, seed);
            if (Double.isFinite(result.tQuarantine[2]) && result.tQuarantine[2] == result.tOnset[0]
                    && result.tInfect[2] < result.tQuarantine[2] && result.tQuarantine[2] < result.tOnset[2]
                    && result.tOnset[4] > result.tOnset[0]) {
                cases++;
                assertTrue(result.tRecover[2] > result.tOnset[2]);
                assertFalse(result.quarantined[3]);
                assertTrue(Double.isNaN(result.tInfect[3]) || result.tInfect[3] <= result.tQuarantine[2]);
            }
        }
        assertTrue(cases > 0);
    }

    @Test
    void isolatedNonObserverOnsetIsNotVisibleAcrossDisconnectedContacts() {
        var graph = graph(5, true, 0, 1, 1, 2, 2, 3, 3, 4);
        int cases = 0;
        for (int seed = 0; seed < 100; seed++) {
            var result =
                    simulate(graph, ModelMode.SPIR, 0, 0, 1, 1, new int[] {0, 2}, new int[] {1, 3}, seed);
            if (result.tOnset[0] < result.tOnset[2]) {
                cases++;
                assertEquals(result.tOnset[0], result.tQuarantine[2]);
                assertTrue(result.tRecover[2] > result.tOnset[2]);
                assertFalse(result.quarantined[3]);
                assertFalse(result.quarantined[4]);
            }
        }
        assertTrue(cases > 0);
    }

    @Test
    void sirInitializationExcludesAllInitialIFromIsolation() {
        var graph = graph(4, true, 0, 1, 0, 2, 0, 3);
        var result = simulate(graph, ModelMode.SIR, 0, 100, 0, 1, new int[] {1, 2}, new int[] {0}, 42);
        assertTrue(result.quarantined[0]);
        assertTrue(result.quarantined[3]);
        assertFalse(result.quarantined[1]);
        assertFalse(result.quarantined[2]);
        assertEquals(0, result.tQuarantine[0]);
        assertEquals(2, result.recovered);
        assertEquals(2, result.finalS);
    }

    @Test
    void newlyInfectedSirNodeTriggersDetectionBeforeItCanTransmitAgain() {
        var graph = graph(4, true, 0, 1, 1, 2, 2, 3);
        var result = simulate(graph, ModelMode.SIR, 0, 100, 0, 1, new int[] {0}, new int[] {2}, 42);
        assertTrue(result.tInfect[1] > 0);
        assertEquals(result.tInfect[1], result.tQuarantine[2]);
        assertTrue(Double.isNaN(result.tInfect[2]));
        assertTrue(Double.isNaN(result.tInfect[3]));
        assertFalse(result.quarantined[1]);
    }

    @Test
    void allObserversActiveAtOneOnsetReactBeforeMonitoringStops() {
        var graph = graph(5, true, 0, 1, 0, 2, 1, 2, 1, 3, 2, 4);
        var result = simulate(graph, ModelMode.SIR, 0, 0, 0, 1, new int[] {0}, new int[] {1, 2}, 42);
        assertEquals(4, result.quarantinedCount);
        assertFalse(result.quarantined[0]);
        assertTrue(result.quarantined[3]);
        assertTrue(result.quarantined[4]);
    }

    @Test
    void monitoringUsesBothIncomingAndOutgoingAdjacency() {
        for (var graph : new DirectedGraph[] {graph(3, false, 1, 0, 1, 2), graph(3, false, 0, 1, 2, 1)}) {
            var result = simulate(graph, ModelMode.SIR, 0, 0, 0, 0, new int[] {0}, new int[] {1}, 42);
            assertTrue(result.quarantined[1]);
            assertTrue(result.quarantined[2]);
            assertFalse(result.quarantined[0]);
        }
    }

    @Test
    void emptyObserverSetPreservesTrajectoryAndInvalidInitialObserverIsRejected() {
        var graph = graph(3, true, 0, 1, 1, 2);
        var old = simulateApi(graph, ModelMode.SPIR, 1, 1, 1, 1, 100, new int[] {0}, 42);
        var current = simulate(graph, ModelMode.SPIR, 1, 1, 1, 1, new int[] {0}, new int[0], 42);
        assertEquals(old.times, current.times);
        assertEquals(old.P, current.P);
        assertArrayEquals(old.tRecover, current.tRecover);
        assertEquals(0, current.quarantinedCount);
        assertThrows(IllegalArgumentException.class,
                () -> simulate(graph, ModelMode.SPIR, 1, 1, 1, 1, new int[] {0}, new int[] {0}, 42));
        assertThrows(IllegalArgumentException.class,
                () -> simulate(graph, ModelMode.SPIR, 1, 1, 1, 1, new int[] {0}, new int[] {3}, 42));
    }

    @Test
    void fixedGridAndCsvReportObserverCountAndCumulativeIsolation() throws Exception {
        var graph = graph(3, true, 0, 1, 1, 2);
        var result = simulateApi(
                graph, ModelMode.SPIR, 0, 0, 1, 1, 100, 0.5, new int[] {0}, new int[] {1}, 42);
        assertEquals(201, result.Q.size());
        for (int k = 0; k < result.Q.size(); k++) {
            assertEquals(result.times.get(k) < result.tOnset[0] ? 0 : 2, result.Q.get(k));
        }
        Path csv = temp.resolve("observer.csv");
        result.writeCsv(csv, 0, 1.0 / 3, 1.0 / 3, ModelMode.SPIR, 0, 0, 1, 1, true, true);
        var lines = Files.readAllLines(csv);
        assertTrue(lines.get(0).endsWith(",phi,observer_count,Q"));
        String[] cells = lines.get(1).split(",", -1);
        assertEquals(15, cells.length);
        assertEquals("1", cells[13]);
        assertEquals("2", cells[14]);
        assertThrows(java.io.IOException.class,
                () -> result.writeCsv(csv, 1, 1.0 / 3, ModelMode.SPIR, 0, 0, 1, 1, true, true));
    }
}
