package sim.simulation;
import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import sim.network.topology.undirected.ER;
class SPIRRefactorTest {
    @Test
    void finalOnlyPreservesTerminalStateWithoutIntermediateHistory() {
        var prepared = new SPIRSimulator.PreparedGraph(ER.generateERFromKAve(100, 4, 42));
        for (var engine : SPIRSimulator.Engine.values()) {
            var history = SPIRSimulator.simulate(prepared, SPIRSimulator.ModelMode.SPIR, engine,
                    SPIRSimulator.RecordingMode.EVENT_TIMES, 1, 1, 1, 1, 100, .1, new int[] {0},
                    new int[] {1}, 42);
            var terminal = SPIRSimulator.simulate(prepared, SPIRSimulator.ModelMode.SPIR, engine,
                    SPIRSimulator.RecordingMode.FINAL_ONLY, 1, 1, 1, 1, 100, .1, new int[] {0}, new int[] {1},
                    42);
            assertEquals(1, terminal.times.size());
            assertEquals(history.finalS, terminal.finalS);
            assertEquals(history.recovered, terminal.recovered);
            assertArrayEquals(history.tRecover, terminal.tRecover);
            double[] retained = terminal.tInfect.clone();
            SPIRSimulator.simulate(prepared, SPIRSimulator.ModelMode.SPIR, engine,
                    SPIRSimulator.RecordingMode.FINAL_ONLY, 4, 4, 1, 1, 100, .1, new int[] {2}, new int[] {3},
                    1);
            assertArrayEquals(retained, terminal.tInfect);
        }
    }
    @Test void legacyApiRetainsEventHistoryAndModeCompatibility() {
        var graph = ER.generateERFromKAve(20,4,42);
        var legacy = SPIRGillespieSimulator.simulate(graph,SPIRGillespieSimulator.ModelMode.SPIR,
                1,1,1,1,100,new int[]{0},42);
        var modern = SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(graph),
                SPIRSimulator.ModelMode.SPIR,SPIRSimulator.Engine.GILLESPIE,
                SPIRSimulator.RecordingMode.EVENT_TIMES,1,1,1,1,100,.1,new int[]{0},new int[0],42);
        assertEquals(legacy.times,modern.times);
        assertArrayEquals(legacy.tInfect,modern.tInfect);
        assertThrows(IllegalArgumentException.class, () -> SPIRGillespieSimulator.simulate(graph,null,
                1,1,1,1,100,new int[]{0},42));
    }

}
