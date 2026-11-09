package sim.simulation;
import sim.network.DirectedGraph;
/** Compatibility entry point for Gillespie event-history and fixed-grid recording. */
public final class SPIRGillespieSimulator {
    private SPIRGillespieSimulator() {}
    public enum ModelMode { SPIR, SIR }
    private static SPIRSimulator.ModelMode canonical(ModelMode mode) {
        if (mode == null) throw new IllegalArgumentException("mode is required");
        return SPIRSimulator.ModelMode.valueOf(mode.name());
    }

    public static SPIRResult simulate(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, double tMax, int[] initial, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(graph),
                canonical(mode), SPIRSimulator.Engine.GILLESPIE,
                SPIRSimulator.RecordingMode.EVENT_TIMES, lambda1, lambda2, mu1, mu2, tMax, 0.1, initial,
                new int[0], seed);
    }

    public static SPIRResult simulate(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, double tMax, double dt, int[] initial, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(graph),
                canonical(mode), SPIRSimulator.Engine.GILLESPIE,
                SPIRSimulator.RecordingMode.FIXED_INTERVAL, lambda1, lambda2, mu1, mu2, tMax, dt, initial,
                new int[0], seed);
    }

    public static SPIRResult simulate(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, double tMax, int[] initial, int[] observers, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(graph),
                canonical(mode), SPIRSimulator.Engine.GILLESPIE,
                SPIRSimulator.RecordingMode.EVENT_TIMES, lambda1, lambda2, mu1, mu2, tMax, 0.1, initial,
                observers, seed);
    }

    public static SPIRResult simulate(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, double tMax, double dt, int[] initial, int[] observers, long seed) {
        return SPIRSimulator.simulate(new SPIRSimulator.PreparedGraph(graph),
                canonical(mode), SPIRSimulator.Engine.GILLESPIE,
                SPIRSimulator.RecordingMode.FIXED_INTERVAL, lambda1, lambda2, mu1, mu2, tMax, dt, initial,
                observers, seed);
    }
}
