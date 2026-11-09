package sim.simulation;

import static org.junit.jupiter.api.Assertions.*;

import java.nio.file.*;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.network.DirectedGraph;
import sim.simulation.SPIRGillespieSimulator.ModelMode;

class SPIRSimulatorTest extends SPIRTestSupport {
    @TempDir Path temp;

    private DirectedGraph graph(boolean undirected) {
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 2, new int[] {0}, new int[] {1}, new boolean[] {undirected});
    }
    private SPIRResult simulate(
            ModelMode mode, double l1, double l2, double m1, double m2, int[] initial, long seed) {
        return simulateApi(graph(false), mode, l1, l2, m1, m2, 100, initial, seed);
    }
    @Test
    void transitionsConservePopulationAndAreOrdered() {
        SPIRResult r = simulate(ModelMode.SPIR, 100, 100, 1, 1, new int[] {0, 0}, 42);
        assertEquals(1, r.P.get(0));
        for (int k = 0; k < r.times.size(); k++) {
            assertEquals(2, r.S.get(k) + r.P.get(k) + r.I.get(k) + r.R.get(k));
            assertTrue(r.S.get(k) >= 0 && r.P.get(k) >= 0 && r.I.get(k) >= 0 && r.R.get(k) >= 0);
        }
        assertEquals(2, r.recovered);
        for (int u = 0; u < 2; u++) {
            assertTrue(r.tInfect[u] <= r.tOnset[u]);
            assertTrue(r.tOnset[u] < r.tRecover[u]);
        }
        assertEquals(100, r.finalTime);
    }
    @Test
    void respectsEdgeDirectionAndUndirectedFlag() {
        SPIRResult forward = simulate(ModelMode.SPIR, 1, 0, 0, 0, new int[] {1}, 42);
        assertTrue(Double.isNaN(forward.tInfect[0]));
        SPIRResult reverse = simulateApi(
                graph(true), ModelMode.SPIR, 1, 0, 0, 0, 100, new int[] {1}, 42);
        assertTrue(Double.isFinite(reverse.tInfect[0]));
    }
    @Test
    void zeroRatesAndEmptyInitialStateRemainStationary() {
        SPIRResult r = simulate(ModelMode.SPIR, 0, 0, 0, 0, new int[] {0}, 42);
        assertEquals(1, r.finalP);
        assertEquals(0, r.finalI);
        r = simulate(ModelMode.SPIR, 1, 1, 1, 1, new int[0], 42);
        assertEquals(2, r.finalS);
        assertTrue(Double.isNaN(r.tInfect[0]));
    }
    @Test
    void switchingToIActivatesLambda2() {
        SPIRResult r = simulate(ModelMode.SPIR, 0, 100, 1, 0, new int[] {0}, 42);
        assertTrue(r.tInfect[1] > r.tOnset[0]);
        assertTrue(Double.isNaN(r.tRecover[0]));
    }
    @Test
    void sirBypassesPAndIgnoresUnusedRates() {
        SPIRResult a = simulate(ModelMode.SIR, 0, 1, 0, 1, new int[] {0}, 42);
        SPIRResult b = simulate(ModelMode.SIR, 99, 1, 99, 1, new int[] {0}, 42);
        assertEquals(a.times, b.times);
        assertEquals(a.I, b.I);
        assertTrue(a.P.stream().allMatch(x -> x == 0));
        assertEquals(0, a.tOnset[0]);
    }
    @Test
    void reproducesSeedsAndFixedGridIncludingTerminalState() {
        SPIRResult a = simulate(ModelMode.SPIR, 1, 1, 1, 1, new int[] {0}, 42);
        SPIRResult b = simulate(ModelMode.SPIR, 1, 1, 1, 1, new int[] {0}, 42);
        assertEquals(a.times, b.times);
        assertArrayEquals(a.tRecover, b.tRecover);
        SPIRResult grid = simulateApi(
                graph(false), ModelMode.SPIR, 0, 0, 1, 1, 3, 0.1, new int[] {0}, 42);
        assertEquals(31, grid.times.size());
        for (int k = 0; k < grid.times.size(); k++) {
            double t = grid.times.get(k);
            int p = t < grid.tOnset[0] ? 1 : 0;
            int i = t >= grid.tOnset[0] && (Double.isNaN(grid.tRecover[0]) || t < grid.tRecover[0]) ? 1 : 0;
            assertEquals(p, grid.P.get(k));
            assertEquals(i, grid.I.get(k));
        }
        SPIRResult shortGrid = simulateApi(
                graph(false), ModelMode.SPIR, 0, 0, 0, 0, 0.35, 0.1, new int[] {0}, 42);
        assertEquals(4, shortGrid.times.size());
    }
    @Test
    void validatesArguments() {
        for (double value : new double[] {-1, Double.NaN, Double.POSITIVE_INFINITY}) {
            assertThrows(IllegalArgumentException.class,
                    () -> simulate(ModelMode.SPIR, value, 1, 1, 1, new int[] {0}, 1));
            assertThrows(IllegalArgumentException.class,
                    () -> simulate(ModelMode.SPIR, 1, value, 1, 1, new int[] {0}, 1));
            assertThrows(IllegalArgumentException.class,
                    () -> simulate(ModelMode.SPIR, 1, 1, value, 1, new int[] {0}, 1));
            assertThrows(IllegalArgumentException.class,
                    () -> simulate(ModelMode.SPIR, 1, 1, 1, value, new int[] {0}, 1));
        }
        assertThrows(
                IllegalArgumentException.class, () -> simulate(ModelMode.SPIR, 1, 1, 1, 1, new int[] {2}, 1));
        assertThrows(IllegalArgumentException.class,
                ()
                        -> simulateApi(
                                graph(false), ModelMode.SPIR, 1, 1, 1, 1, 0, new int[] {0}, 1));
        assertThrows(IllegalArgumentException.class,
                ()
                        -> simulateApi(
                                graph(false), ModelMode.SPIR, 1, 1, 1, 1, 1, 0, new int[] {0}, 1));
    }
    @Test
    void oneEdgeMatchesAnalyticInfectionProbability() {
        for (ModelMode mode : ModelMode.values()) {
            int infected = 0, samples = 10000;
            Random seeds = new Random(7331);
            for (int k = 0; k < samples; k++) {
                SPIRResult r = simulate(mode, 1, 1, 1, 1, new int[] {0}, seeds.nextLong());
                if (Double.isFinite(r.tInfect[1]))
                    infected++;
            }
            assertEquals(mode == ModelMode.SIR ? 0.5 : 0.75, (double) infected / samples, 0.025);
        }
    }
    @Test
    void csvPreservesSchemaAndBlankUnusedRates() throws Exception {
        Path path = temp.resolve("results.csv");
        SPIRResult r = simulate(ModelMode.SIR, 0, 1, 0, 1, new int[] {0}, 42);
        for (int itr = 0; itr < 2; itr++) r.writeCsv(path, itr, 0.5, ModelMode.SIR, 0, 1, 0, 1, true, true);
        var lines = Files.readAllLines(path);
        assertEquals(3, lines.size());
        assertEquals("model,itr,rho_0,lambda1,lambda2,mu1,mu2,time,S,P,I,R", lines.get(0));
        String[] cells = lines.get(1).split(",", -1);
        assertEquals(12, cells.length);
        assertEquals("", cells[3]);
        assertEquals("", cells[5]);
        assertEquals("100.000000000", cells[7]);
        assertEquals("0", cells[9]);
        Path trajectory = temp.resolve("times.csv");
        r.writeCsv(trajectory, 0, 0.5, ModelMode.SIR, 0, 1, 0, 1, false, true);
        assertEquals(r.times.size() + 1, Files.readAllLines(trajectory).size());
    }
}
