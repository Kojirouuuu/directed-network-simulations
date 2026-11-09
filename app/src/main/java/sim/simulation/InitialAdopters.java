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

    /** Sample from non-observers, retaining the count relative to the full population. */
    public static int[] sample(int n, double rho0, boolean[] observers, long seed) {
        if (n <= 0 || observers == null || observers.length != n
                || !Double.isFinite(rho0) || rho0 < 0 || rho0 > 1) {
            throw new IllegalArgumentException("positive n, observer mask of length n and rho0 in [0,1] required");
        }
        int count = Math.max(1, (int) (n * rho0));
        int[] nodes = IntStream.range(0, n).filter(node -> !observers[node]).toArray();
        if (count > nodes.length) {
            throw new IllegalArgumentException("not enough non-observers for " + count + " initial infections");
        }
        ArrayUtils.shuffle(nodes, seed);
        return Arrays.copyOf(nodes, count);
    }
}
