package sim.network.topology;

import java.util.Arrays;
import java.util.Random;

import sim.network.DirectedGraph;

/**
 * Erdős–Rényi 型の単純有向グラフを生成する。
 *
 * <p>自己ループを除く各順序付き頂点対 {@code u -> v} を、他の辺とは独立に
 * 確率 {@code p} で生成する。したがって各頂点の入次数と出次数は
 * {@code Binomial(n - 1, p)} に従い、疎なグラフでは平均 {@code (n - 1) * p} の
 * ポアソン分布に近似する。</p>
 */
public final class DirectedER {

    private static final String GRAPH_NAME = "DirectedER";
    private static final int MAX_ARRAY_LENGTH = Integer.MAX_VALUE - 8;

    private DirectedER() {
    }

    /**
     * 各順序付き頂点対を確率 {@code p} で結ぶ。
     *
     * @param n 頂点数
     * @param p 有向辺の生成確率
     * @param seed 乱数シード
     * @return 自己ループと多重辺を持たない有向グラフ
     */
    public static DirectedGraph generateFromP(int n, double p, long seed) {
        validateProbabilityParameters(n, p);

        long possibleEdgeCount = (long) n * (n - 1L);
        if (possibleEdgeCount == 0L || p == 0.0) {
            return buildGraph(n, new int[0], new int[0]);
        }
        if (p == 1.0) {
            if (possibleEdgeCount > MAX_ARRAY_LENGTH) {
                throw new IllegalArgumentException("Generated edge count exceeds the supported int array size");
            }
            int edgeCount = (int) possibleEdgeCount;
            int[] sources = new int[edgeCount];
            int[] destinations = new int[edgeCount];
            int edge = 0;
            for (int source = 0; source < n; source++) {
                for (int destination = 0; destination < n; destination++) {
                    if (source != destination) {
                        sources[edge] = source;
                        destinations[edge] = destination;
                        edge++;
                    }
                }
            }
            return buildGraph(n, sources, destinations);
        }

        EdgeBuffer edges = new EdgeBuffer(initialCapacity(possibleEdgeCount, p));
        Random random = new Random(seed);
        double logFailureProbability = Math.log1p(-p);
        long position = -1L;
        int destinationsPerSource = n - 1;

        while (position < possibleEdgeCount - 1L) {
            double failures = Math.floor(Math.log1p(-random.nextDouble()) / logFailureProbability);
            long remainingPositions = possibleEdgeCount - position - 1L;
            if (failures >= remainingPositions) {
                break;
            }
            position += (long) failures + 1L;

            int source = (int) (position / destinationsPerSource);
            int destinationIndex = (int) (position % destinationsPerSource);
            int destination = destinationIndex >= source ? destinationIndex + 1 : destinationIndex;
            edges.add(source, destination);
        }

        return buildGraph(n, edges.sources(), edges.destinations());
    }

    /** 乱数シードを現在時刻から与える確率指定版。 */
    public static DirectedGraph generateFromP(int n, double p) {
        return generateFromP(n, p, System.currentTimeMillis());
    }

    /**
     * 期待平均入次数および出次数を指定して生成する。
     *
     * @param n 頂点数
     * @param meanDegree 期待平均入次数および出次数
     * @param seed 乱数シード
     * @return directed ER グラフ
     */
    public static DirectedGraph generateFromMeanDegree(int n, double meanDegree, long seed) {
        validateMeanDegreeParameters(n, meanDegree);
        double p = n == 1 ? 0.0 : meanDegree / (n - 1.0);
        return generateFromP(n, p, seed);
    }

    /** 乱数シードを現在時刻から与える平均次数指定版。 */
    public static DirectedGraph generateFromMeanDegree(int n, double meanDegree) {
        return generateFromMeanDegree(n, meanDegree, System.currentTimeMillis());
    }

    private static void validateProbabilityParameters(int n, double p) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1");
        }
        if (!Double.isFinite(p) || p < 0.0 || p > 1.0) {
            throw new IllegalArgumentException("p must be finite and in [0, 1]");
        }
    }

    private static void validateMeanDegreeParameters(int n, double meanDegree) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be >= 1");
        }
        if (!Double.isFinite(meanDegree) || meanDegree < 0.0 || meanDegree > n - 1.0) {
            throw new IllegalArgumentException("meanDegree must be finite and in [0, n - 1]");
        }
    }

    private static int initialCapacity(long possibleEdgeCount, double p) {
        double expectedEdgeCount = possibleEdgeCount * p;
        return (int) Math.min(1_000_000.0, Math.ceil(expectedEdgeCount));
    }

    private static DirectedGraph buildGraph(int n, int[] sources, int[] destinations) {
        return DirectedGraph.fromEdgeListWithUndirectedFlag(
                GRAPH_NAME, n, sources, destinations, new boolean[sources.length]);
    }

    private static final class EdgeBuffer {
        private int[] sources;
        private int[] destinations;
        private int size;

        private EdgeBuffer(int initialCapacity) {
            sources = new int[initialCapacity];
            destinations = new int[initialCapacity];
        }

        private void add(int source, int destination) {
            if (size == sources.length) {
                grow();
            }
            sources[size] = source;
            destinations[size] = destination;
            size++;
        }

        private void grow() {
            if (size == MAX_ARRAY_LENGTH) {
                throw new IllegalArgumentException("Generated edge count exceeds the supported int array size");
            }
            int newLength = size == 0 ? 1 : size + Math.max(1, size / 2);
            if (newLength < 0 || newLength > MAX_ARRAY_LENGTH) {
                newLength = MAX_ARRAY_LENGTH;
            }
            sources = Arrays.copyOf(sources, newLength);
            destinations = Arrays.copyOf(destinations, newLength);
        }

        private int[] sources() {
            return Arrays.copyOf(sources, size);
        }

        private int[] destinations() {
            return Arrays.copyOf(destinations, size);
        }
    }
}
