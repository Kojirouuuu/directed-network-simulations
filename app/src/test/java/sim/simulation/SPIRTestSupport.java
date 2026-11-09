package sim.simulation;
import sim.network.DirectedGraph;
import sim.simulation.SPIRGillespieSimulator.ModelMode;
abstract class SPIRTestSupport {
    protected SPIRSimulator.Engine engine() {
        return SPIRSimulator.Engine.GILLESPIE;
    }
    protected SPIRResult simulateApi(DirectedGraph g, ModelMode mode, double l1, double l2, double m1,
            double m2, double end, int[] initial, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(g),
                SPIRSimulator.ModelMode.valueOf(mode.name()), engine(),
                SPIRSimulator.RecordingMode.EVENT_TIMES, l1, l2, m1, m2, end, 0.1, initial, new int[0], seed);
    }
    protected SPIRResult simulateApi(DirectedGraph g, ModelMode mode, double l1, double l2, double m1,
            double m2, double end, double dt, int[] initial, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(g),
                SPIRSimulator.ModelMode.valueOf(mode.name()), engine(),
                SPIRSimulator.RecordingMode.FIXED_INTERVAL, l1, l2, m1, m2, end, dt, initial, new int[0],
                seed);
    }
    protected SPIRResult simulateApi(DirectedGraph g, ModelMode mode, double l1, double l2, double m1,
            double m2, double end, int[] initial, int[] observers, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(g),
                SPIRSimulator.ModelMode.valueOf(mode.name()), engine(),
                SPIRSimulator.RecordingMode.EVENT_TIMES, l1, l2, m1, m2, end, 0.1, initial, observers, seed);
    }
    protected SPIRResult simulateApi(DirectedGraph g, ModelMode mode, double l1, double l2, double m1,
            double m2, double end, double dt, int[] initial, int[] observers, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(g),
                SPIRSimulator.ModelMode.valueOf(mode.name()), engine(),
                SPIRSimulator.RecordingMode.FIXED_INTERVAL, l1, l2, m1, m2, end, dt, initial, observers,
                seed);
    }
}
