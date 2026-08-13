package sim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import sim.network.DirectedGraph;
import sim.simulation.SARResult;
import sim.simulation.SARSimulator;

class NetworkSweepSARTest {

    @TempDir
    Path tempDir;

    @Test
    void validatesApplicableSweepAxes() {
        NetworkSweepSAR.validateSweepAxes(
                "PowPow", new double[] { 2.3, 2.7 }, new int[] { 0, 100 },
                NetworkSweepSAR.RewiringMode.RECIPROCITY,
                new double[] { 0.0, 0.2 }, null);
        NetworkSweepSAR.validateSweepAxes(
                "PowPow", new double[] { 2.3, 2.7 }, new int[] { 0, 100 },
                NetworkSweepSAR.RewiringMode.FFL,
                null, new long[] { 0L, 1_000L });

        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "ER", new double[] { 2.3, 2.7 }, new int[] { 0 },
                        NetworkSweepSAR.RewiringMode.RECIPROCITY, new double[] { 0.0 }, null));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "SameInOut", new double[] { 2.5 }, new int[] { 0, 1 },
                        NetworkSweepSAR.RewiringMode.RECIPROCITY, new double[] { 0.0 }, null));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "PowPow", new double[0], new int[] { 0 },
                        NetworkSweepSAR.RewiringMode.RECIPROCITY, new double[] { 0.0 }, null));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "PowPow", new double[] { 2.5 }, new int[] { 0 },
                        NetworkSweepSAR.RewiringMode.RECIPROCITY, new double[] { 1.1 }, null));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "PowPow", new double[] { 2.5 }, new int[] { 0 },
                        NetworkSweepSAR.RewiringMode.FFL, null, new long[] { -1L }));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.validateSweepAxes(
                        "PowPow", new double[] { 2.5 }, new int[] { 0 },
                        NetworkSweepSAR.RewiringMode.FFL, null, new long[] { 1L, 1L }));

        assertEquals(2, NetworkSweepSAR.rewiringTargetCount(
                NetworkSweepSAR.RewiringMode.RECIPROCITY,
                new double[] { 0.0, 0.2 }, new long[] { 1L }));
        assertEquals(3, NetworkSweepSAR.rewiringTargetCount(
                NetworkSweepSAR.RewiringMode.FFL,
                new double[] { 0.0 }, new long[] { 1L, 2L, 3L }));
    }

    @Test
    void scalesAttemptBudgetsBeyondIntegerRange() {
        assertEquals(1_000L, NetworkSweepSAR.scaledAttempts(100, 10, "testAttempts"));
        assertEquals(2L * Integer.MAX_VALUE,
                NetworkSweepSAR.scaledAttempts(Integer.MAX_VALUE, 2, "testAttempts"));
        assertThrows(IllegalArgumentException.class,
                () -> NetworkSweepSAR.scaledAttempts(100, -1, "testAttempts"));
    }

    @Test
    void writesNetworkParametersAndMeasuredReciprocityToCsv() throws Exception {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 2, new int[] { 0 }, new int[] { 1 }, new boolean[] { false });
        SARResult result = SARSimulator.simulate(
                graph, 0.0, 0.0, 1.0, 1.0,
                new int[] { 1, 1 }, new int[] { 0 }, 123L);
        Path output = tempDir.resolve("results.csv");

        NetworkSweepSAR.writeResultCsv(
                output, result, true, 3, "PowPow", 2.5, 100,
                0.2, 0.25, null, null,
                0.01, 2.0, 0.0, 1.0, true);

        List<String> lines = Files.readAllLines(output);
        assertEquals(2, lines.size());
        assertEquals(
                "itr,network_type,gamma,swap_num,target_r,actual_r,target_n_ffl,actual_n_ffl,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi",
                lines.get(0));
        assertTrue(lines.get(1).startsWith(
                "3,PowPow,2.500000000,100,0.200000000,0.250000000,,,0.010000000,2.000000000,0.000000000,1.000000000,"));
    }

    @Test
    void writesTargetAndMeasuredFflToCsv() throws Exception {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 3, new int[] { 0, 1, 0 }, new int[] { 1, 2, 2 }, new boolean[3]);
        SARResult result = SARSimulator.simulate(
                graph, 0.0, 0.0, 1.0, 1.0,
                new int[] { 1, 1, 1 }, new int[] { 0 }, 123L);
        Path output = tempDir.resolve("ffl-results.csv");

        NetworkSweepSAR.writeResultCsv(
                output, result, true, 4, "PowPow", 2.5, 0,
                null, 0.125, 10L, 12L,
                0.01, 2.0, 0.0, 1.0, true);

        List<String> lines = Files.readAllLines(output);
        assertEquals(2, lines.size());
        assertTrue(lines.get(1).startsWith(
                "4,PowPow,2.500000000,0,,0.125000000,10,12,0.010000000,2.000000000,0.000000000,1.000000000,"));
    }
}
