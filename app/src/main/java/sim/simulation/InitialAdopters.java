package sim.simulation;

import java.util.Arrays;
import java.util.stream.IntStream;
import sim.utils.ArrayUtils;

/** Deterministic initial adopter selection shared by experiment runners. */
public final class InitialAdopters {
    private InitialAdopters() {}

    public static int[] sample(int n, double rho0, long seed) {
        int[] nodes = IntStream.range(0, n).toArray();
        ArrayUtils.shuffle(nodes, seed);
        return Arrays.copyOfRange(nodes, 0, Math.max(1, (int) (n * rho0)));
    }
}
