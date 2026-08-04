package sim.network.topology;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Consumer;

import sim.network.DirectedGraph;

/**
 * 同じ実現平均次数を持つ Pow/Pow、Pow/Poi、Poi/Pow、Poi/Poi の
 * 有向コンフィギュレーションモデルを一組として生成する。
 */
public final class DegreeDistributionNetworks {

    private static final long SEED_MIX_K2 = 0x9E3779B97F4A7C15L;
    private static final long SEED_MIX_BALANCE = 0xDEADBEEFCAFEBABEL;
    private static final long SEED_MIX_POISSON_1 = 0x517CC1B727220A95L;
    private static final long SEED_MIX_POISSON_2 = 0x1234567890ABCDEFL;
    private static final long SEED_MIX_GRAPH = 0xD1B54A32D192ED03L;

    private DegreeDistributionNetworks() {
    }

    /** 次数分布の表示名。CSVにもこの値を書き出す。 */
    public enum Distribution {
        POW("Pow"),
        POI("Poi");

        private final String label;

        Distribution(String label) {
            this.label = label;
        }

        public String label() {
            return label;
        }
    }

    /** 一つの入次数・出次数分布の組み合わせと、そのグラフ。 */
    public record Variant(
            Distribution inDistribution,
            Distribution outDistribution,
            DirectedGraph graph) {
    }

    /**
     * 4種類のネットワークを Pow/Pow、Pow/Poi、Poi/Pow、Poi/Poi の順で返す。
     * Poi/Pow は Pow/Poi の辺をすべて反転したグラフである。
     */
    public static List<Variant> generate(
            String name, int n, int kdMin, int kdMax, double gamma, long seed) {
        List<Variant> variants = new ArrayList<>(4);
        forEachVariant(name, n, kdMin, kdMax, gamma, seed, variants::add);
        return List.copyOf(variants);
    }

    /**
     * 4種類のネットワークを順番に生成して処理する。
     * 大規模実験で4グラフを同時に保持しないためのAPI。
     */
    public static void forEachVariant(
            String name, int n, int kdMin, int kdMax, double gamma, long seed,
            Consumer<Variant> consumer) {
        validateParams(name, n, kdMin, kdMax, gamma);
        if (consumer == null) {
            throw new IllegalArgumentException("consumer must be non-null");
        }

        int[] k1 = samplePowerLaw(n, kdMin, kdMax, gamma, seed);
        int[] k2 = samplePowerLaw(n, kdMin, kdMax, gamma, seed ^ SEED_MIX_K2);
        long edgeCount = sum(k1);
        if (edgeCount > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Sum of directed degrees must be <= Integer.MAX_VALUE");
        }
        balanceToTarget(k2, edgeCount, kdMin, kdMax, seed ^ SEED_MIX_BALANCE);

        int[] p1 = assignRandomly(n, (int) edgeCount, seed ^ SEED_MIX_POISSON_1);
        int[] p2 = assignRandomly(n, (int) edgeCount, seed ^ SEED_MIX_POISSON_2);
        int[] noUndirectedEdges = new int[n];

        DirectedGraph powPow = DirectedCMOutPow.generateFromDegreeSequence(
                name + "_Pow_Pow", k1, k2, noUndirectedEdges, seed ^ SEED_MIX_GRAPH);
        consumer.accept(new Variant(Distribution.POW, Distribution.POW, powPow));

        DirectedGraph powPoi = DirectedCMOutPow.generateFromDegreeSequence(
                name + "_Pow_Poi", k1, p1, noUndirectedEdges, seed ^ SEED_MIX_GRAPH ^ SEED_MIX_POISSON_1);
        consumer.accept(new Variant(Distribution.POW, Distribution.POI, powPoi));

        DirectedGraph poiPow = powPoi.reverseDirectedEdges(name + "_Poi_Pow");
        consumer.accept(new Variant(Distribution.POI, Distribution.POW, poiPow));

        DirectedGraph poiPoi = DirectedCMOutPow.generateFromDegreeSequence(
                name + "_Poi_Poi", p1, p2, noUndirectedEdges, seed ^ SEED_MIX_GRAPH ^ SEED_MIX_POISSON_2);
        consumer.accept(new Variant(Distribution.POI, Distribution.POI, poiPoi));
    }

    private static void validateParams(String name, int n, int kdMin, int kdMax, double gamma) {
        if (name == null) {
            throw new IllegalArgumentException("name must be non-null");
        }
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1");
        }
        if (kdMin < 1) {
            throw new IllegalArgumentException("kdMin must be >= 1");
        }
        if (kdMax < kdMin) {
            throw new IllegalArgumentException("kdMax must be >= kdMin");
        }
        if (!Double.isFinite(gamma) || gamma <= 1.0) {
            throw new IllegalArgumentException("gamma must be finite and > 1.0");
        }
    }

    private static int[] samplePowerLaw(int n, int kdMin, int kdMax, double gamma, long seed) {
        double[] cdf = buildPowerLawCdf(kdMin, kdMax, gamma);
        Random random = new Random(seed);
        int[] degrees = new int[n];
        for (int u = 0; u < n; u++) {
            int index = Arrays.binarySearch(cdf, random.nextDouble());
            degrees[u] = (index >= 0 ? index : -index - 1) + kdMin;
        }
        return degrees;
    }

    private static double[] buildPowerLawCdf(int kdMin, int kdMax, double gamma) {
        int length = kdMax - kdMin + 1;
        double[] cdf = new double[length];
        double total = 0.0;
        for (int k = kdMin; k <= kdMax; k++) {
            total += Math.pow(k, -gamma);
            cdf[k - kdMin] = total;
        }
        for (int i = 0; i < length; i++) {
            cdf[i] /= total;
        }
        cdf[length - 1] = 1.0;
        return cdf;
    }

    /** k1 は固定したまま、k2 のみを1ずつ調整して次数和を一致させる。 */
    private static void balanceToTarget(
            int[] degrees, long target, int kdMin, int kdMax, long seed) {
        Random random = new Random(seed);
        long current = sum(degrees);
        while (current != target) {
            boolean increment = current < target;
            int start = random.nextInt(degrees.length);
            int selected = -1;
            for (int offset = 0; offset < degrees.length; offset++) {
                int u = (start + offset) % degrees.length;
                if ((increment && degrees[u] < kdMax) || (!increment && degrees[u] > kdMin)) {
                    selected = u;
                    break;
                }
            }
            if (selected < 0) {
                throw new IllegalStateException("Unable to balance power-law degree sequence");
            }
            degrees[selected] += increment ? 1 : -1;
            current += increment ? 1 : -1;
        }
    }

    /** 固定された総次数を頂点へ一様ランダムに配分する。 */
    private static int[] assignRandomly(int n, int totalDegree, long seed) {
        Random random = new Random(seed);
        int[] degrees = new int[n];
        for (int i = 0; i < totalDegree; i++) {
            degrees[random.nextInt(n)]++;
        }
        return degrees;
    }

    private static long sum(int[] values) {
        long total = 0L;
        for (int value : values) {
            total += value;
        }
        return total;
    }
}
