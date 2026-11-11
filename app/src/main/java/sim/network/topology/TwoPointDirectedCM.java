package sim.network.topology;

import java.util.Random;
import sim.network.DirectedGraph;

/** Two-point joint-degree configuration model; loops and parallel arcs are allowed. */
public final class TwoPointDirectedCM {
    private TwoPointDirectedCM() {}

    public record DegreePairCounts(int n11, int n12, int n21, int n22) {
        public double realizedP() {
            return (double) (n11 + n12) / (n11 + n12 + n21 + n22);
        }
        public double realizedR() {
            double n = (double) n11 + n12 + n21 + n22;
            return 1.0 - n * n12 / ((double) (n11 + n12) * (n21 + n22));
        }
    }

    public static void validateProbabilities(double p, double r) {
        if (!Double.isFinite(p) || p <= 0 || p >= 1)
            throw new IllegalArgumentException("p must be finite and strictly between 0 and 1");
        double lower = -Math.min(p / (1 - p), (1 - p) / p);
        if (!Double.isFinite(r) || r < lower || r > 1)
            throw new IllegalArgumentException("r must be finite and in [" + lower + ", 1]");
    }

    /** Fixed counts, half-up rounding; off-diagonal counts are equal to balance stubs. */
    public static DegreePairCounts degreePairCounts(int n, double p, double r) {
        validateProbabilities(p, r);
        if (n < 2) throw new IllegalArgumentException("n must be >= 2");
        int n1 = (int) Math.floor(n * p + 0.5);
        int n2 = n - n1;
        if (n1 == 0 || n2 == 0)
            throw new IllegalArgumentException("n is too small to represent both degree classes");
        double pn = (double) n1 / n;
        int b = (int) Math.floor(n * pn * (1 - pn) * (1 - r) + 0.5);
        b = Math.max(0, Math.min(Math.min(n1, n2), b));
        return new DegreePairCounts(n1 - b, b, b, n2 - b);
    }

    public static DirectedGraph generate(String name, int n, int k1, int k2,
            double p, double r, long seed) {
        if (name == null) throw new IllegalArgumentException("name must be non-null");
        if (k1 < 1 || k2 <= k1) throw new IllegalArgumentException("require 1 <= k1 < k2");
        var c = degreePairCounts(n, p, r);
        long edges = (long) (c.n11 + c.n12) * k1 + (long) (c.n21 + c.n22) * k2;
        if (edges > Integer.MAX_VALUE) throw new IllegalArgumentException("edge count exceeds int capacity");
        int[] ki = new int[n];
        int[] ko = new int[n];
        int[] counts = {c.n11, c.n12, c.n21, c.n22};
        int pos = 0;
        for (int type = 0; type < 4; type++) {
            for (int t = 0; t < counts[type]; t++, pos++) {
                ki[pos] = type < 2 ? k1 : k2;
                ko[pos] = type % 2 == 0 ? k1 : k2;
            }
        }
        Random rng = new Random(seed);
        for (int i = n - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            int tmp = ki[i]; ki[i] = ki[j]; ki[j] = tmp;
            tmp = ko[i]; ko[i] = ko[j]; ko[j] = tmp;
        }
        return DirectedCM.generateFromDegreeSequence(name, ki, ko, new int[n],
                seed ^ 0x9E3779B97F4A7C15L);
    }
}
