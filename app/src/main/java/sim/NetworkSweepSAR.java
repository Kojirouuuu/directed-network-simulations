package sim;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;

import sim.network.DirectedGraph;
import sim.simulation.SARGillespieSimulator;
import sim.simulation.SARResult;
import sim.simulation.SARSimulator;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;
import sim.utils.SwitchUtils;

/**
 * ネットワーク生成パラメータと相互辺割合を変えながら、ネットワーク生成とSAR実験を行う。
 * 設定は {@link SimulationConfig} の配列をソース上で編集する。
 */
public final class NetworkSweepSAR {

    private static final int PROGRESS_BAR_LENGTH = 100;
    private static final int PROGRESS_UPDATE_INTERVAL_MS = 100;

    private static final long RNG_BASE_SEED = 7L;
    private static final long SIM_BASE_SEED = 12345L;
    private static final long GRAPH_BASE_SEED = 42L;
    private static final long RECIPROCITY_BASE_SEED = 98765L;
    private static final long SEED_OFFSET_NODES = 2000L;

    private NetworkSweepSAR() {
    }

    public static void main(String[] args) throws Exception {
        SimulationConfig config = new SimulationConfig();
        validateSweepAxes(
                config.networkType, config.gammaList, config.swapNumList, config.targetReciprocityList);
        validateSimulationConfig(config);

        long totalTasks = (long) config.batchSize * config.gammaList.length * config.swapNumList.length
                * config.targetReciprocityList.length * config.itrs * config.rho0List.length
                * config.lambdaDirectedList.length * config.lambdaNondirectedList.length;

        System.out.println("Total tasks: " + totalTasks);
        System.out.println(config.networkType + ": N=" + config.N + ", batches=" + config.batchSize);

        int parallelism = Runtime.getRuntime().availableProcessors();
        System.out.println("Parallelism: " + parallelism + " (available processors)");

        AtomicLong done = new AtomicLong(0L);
        AtomicBoolean running = new AtomicBoolean(true);
        Thread renderer = createTotalProgressRenderer(done, totalTasks, running);
        renderer.start();

        try (ForkJoinPool pool = new ForkJoinPool(parallelism)) {
            Future<?> future = pool.submit(() -> IntStream.range(0, config.batchSize).parallel()
                    .forEach(batchIndex -> processBatch(batchIndex, config, done)));
            future.get();
        } finally {
            running.set(false);
            renderer.join();
        }

        System.out.println("All tasks completed");
    }

    private static void processBatch(int batchIndex, SimulationConfig config, AtomicLong done) {
        Path resultsPath = null;
        boolean gammaApplicable = usesGamma(config.networkType);
        boolean swapApplicable = usesSwapNum(config.networkType);

        for (double gamma : config.gammaList) {
            for (int swapNum : config.swapNumList) {
                DirectedGraph baseGraph = generateBaseGraph(config, gamma, swapNum, batchIndex);
                if (resultsPath == null) {
                    resultsPath = prepareOutputPath(baseGraph.n, batchIndex, config);
                }

                long maxIncreaseAttempts = scaledAttempts(
                        baseGraph.m, config.maxIncreaseAttemptsPerEdge, "maxIncreaseAttempts");
                long neutralSwapAttempts = scaledAttempts(
                        baseGraph.m, config.neutralSwapAttemptsPerEdge, "neutralSwapAttempts");
                double baseReciprocity = baseGraph.reciprocity();

                for (double targetReciprocity : config.targetReciprocityList) {
                    long neutralAttemptsForTarget = targetReciprocity > baseReciprocity
                            ? neutralSwapAttempts
                            : 0L;
                    DirectedGraph graph = baseGraph.increaseReciprocity(
                            targetReciprocity,
                            maxIncreaseAttempts,
                            neutralAttemptsForTarget,
                            RECIPROCITY_BASE_SEED + batchIndex);
                    double actualReciprocity = graph.reciprocity();

                    Path edgeListPath = buildEdgeListPath(
                            graph, config, gamma, swapNum, targetReciprocity, batchIndex);
                    try {
                        graph.writeEdgeList(edgeListPath);
                    } catch (IOException e) {
                        throw new RuntimeException(
                                "Edge-list output error (batch=" + batchIndex
                                        + ", gamma=" + gamma + ", swapNum=" + swapNum
                                        + ", targetReciprocity=" + targetReciprocity + ")",
                                e);
                    }

                    synchronized (System.out) {
                        System.out.printf(Locale.ROOT,
                                "%nGenerated %s (target_r=%.6f, actual_r=%.6f)%n",
                                edgeListPath, targetReciprocity, actualReciprocity);
                    }
                    runAllSimulations(
                            graph, config, gamma, swapNum, targetReciprocity, actualReciprocity,
                            gammaApplicable, swapApplicable, batchIndex, resultsPath, done);
                }
            }
        }
    }

    private static Path buildEdgeListPath(
            DirectedGraph graph, SimulationConfig config, double gamma, int swapNum,
            double targetReciprocity, int batchIndex) {
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, graph.n,
                config.kdAve, config.kuAve,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kdMin, config.kdMax, config.kuMin, config.kuMax, config.m0, config.m,
                usesGamma(config.networkType) ? gamma : null,
                usesSwapNum(config.networkType) ? swapNum : null,
                config.gammaIn, config.gammaOut, config.corrA);
        return Paths.get("out/edgelist")
                .resolve(networkPath)
                .resolve(String.format(Locale.ROOT, "reciprocity=%.6f", targetReciprocity))
                .resolve(batchIndex + ".csv");
    }

    private static DirectedGraph generateBaseGraph(
            SimulationConfig config, double gamma, int swapNum, int batchIndex) {
        return SwitchUtils.generateGraph(
                config.networkType, config.N,
                config.kdAve, config.kdMin, config.kdMax,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kuMin, config.kuMax,
                config.kuAve, gamma, config.m0, config.m, swapNum,
                config.gammaIn, config.gammaOut, config.corrA,
                GRAPH_BASE_SEED + batchIndex);
    }

    private static void runAllSimulations(
            DirectedGraph graph, SimulationConfig config, double gamma, int swapNum,
            double targetReciprocity, double actualReciprocity,
            boolean gammaApplicable, boolean swapApplicable,
            int batchIndex, Path resultsPath, AtomicLong done) {
        int[] thresholdList = new int[graph.n];
        Arrays.fill(thresholdList, config.threshold);

        for (int itr = 0; itr < config.itrs; itr++) {
            for (double rho0 : config.rho0List) {
                for (double lambdaDirected : config.lambdaDirectedList) {
                    for (double lambdaNondirected : config.lambdaNondirectedList) {
                        SARResult result = simulate(
                                graph, config, lambdaDirected, lambdaNondirected, rho0,
                                thresholdList, batchIndex, itr);
                        int resultItr = batchIndex * config.itrs + itr;
                        try {
                            writeResultCsv(
                                    resultsPath, result, config.isFinal, resultItr, config.networkType,
                                    gammaApplicable ? gamma : null,
                                    swapApplicable ? swapNum : null,
                                    targetReciprocity, actualReciprocity,
                                    rho0, lambdaDirected, lambdaNondirected, config.mu, true);
                        } catch (IOException e) {
                            throw new RuntimeException(
                                    "CSV output error (batch=" + batchIndex + ", iteration=" + itr
                                            + ", gamma=" + gamma + ", swapNum=" + swapNum
                                            + ", targetReciprocity=" + targetReciprocity + ")",
                                    e);
                        }
                        done.incrementAndGet();
                    }
                }
            }
        }
    }

    private static SARResult simulate(
            DirectedGraph graph, SimulationConfig config,
            double lambdaDirected, double lambdaNondirected, double rho0,
            int[] thresholdList, int batchIndex, int itr) {
        int initialInfectedNum = Math.max(1, (int) (graph.n * rho0));
        int[] nodes = new int[graph.n];
        for (int i = 0; i < graph.n; i++) {
            nodes[i] = i;
        }

        long nodeSeed = RNG_BASE_SEED + (long) batchIndex * config.itrs + itr + SEED_OFFSET_NODES;
        ArrayUtils.shuffle(nodes, nodeSeed);
        int[] init = Arrays.copyOfRange(nodes, 0, initialInfectedNum);
        long simSeed = SIM_BASE_SEED + (long) batchIndex * config.itrs + itr;

        return config.useGillespie
                ? SARGillespieSimulator.simulate(
                        graph, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                        thresholdList, init, simSeed)
                : SARSimulator.simulate(
                        graph, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                        thresholdList, init, simSeed);
    }

    private static Path prepareOutputPath(int graphSize, int batchIndex, SimulationConfig config) {
        String index = String.format("%02d", batchIndex);
        Path outputDir = SwitchUtils.buildSimulationOutputDir(config.optionPath, config.threshold);
        Path basePath = outputDir
                .resolve(config.networkType)
                .resolve("N=" + graphSize)
                .resolve("network-sweep");
        return PathsEx.resolveIndexed(basePath.resolve(String.format("results_%s.csv", index)));
    }

    static void writeResultCsv(
            Path path, SARResult result, boolean finalState, int itr, String networkType,
            Double gamma, Integer swapNum, double targetReciprocity, double actualReciprocity,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        if (!append) {
            path = PathsEx.resolveIndexed(path);
        }
        Files.createDirectories(path.getParent());
        boolean writeHeader = !Files.exists(path) || Files.size(path) == 0L;

        try (BufferedWriter writer = Files.newBufferedWriter(
                path,
                StandardOpenOption.CREATE,
                append ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING);
                PrintWriter out = new PrintWriter(writer)) {
            if (writeHeader) {
                out.println(finalState
                        ? "itr,network_type,gamma,swap_num,target_r,actual_r,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi"
                        : "itr,network_type,gamma,swap_num,target_r,actual_r,rho_0,lambda_d,lambda_u,mu,time,A,R,Phi");
            }

            String gammaValue = gamma == null ? "" : String.format(Locale.ROOT, "%.9f", gamma);
            String swapValue = swapNum == null ? "" : Integer.toString(swapNum);
            if (finalState) {
                int last = result.times.size() - 1;
                out.printf(Locale.ROOT,
                        "%d,%s,%s,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                        itr, networkType, gammaValue, swapValue, targetReciprocity, actualReciprocity,
                        rho0, lambdaDirected, lambdaNondirected, mu, result.times.get(last),
                        result.initialAdoptedTime, result.finalAdoptedTime,
                        result.A.get(last), result.R.get(last), result.Phi.get(last));
            } else {
                for (int i = 0; i < result.times.size(); i++) {
                    out.printf(Locale.ROOT,
                            "%d,%s,%s,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                            itr, networkType, gammaValue, swapValue, targetReciprocity, actualReciprocity,
                            rho0, lambdaDirected, lambdaNondirected, mu, result.times.get(i),
                            result.A.get(i), result.R.get(i), result.Phi.get(i));
                }
            }
        }
    }

    static boolean usesGamma(String networkType) {
        return switch (networkType) {
            case "DirectedCMInPow", "DirectedCMOutPow", "PowPow", "SameInOut", "CM" -> true;
            default -> false;
        };
    }

    static boolean usesSwapNum(String networkType) {
        return "PowPow".equals(networkType);
    }

    static void validateSweepAxes(
            String networkType, double[] gammaList, int[] swapNumList, double[] targetReciprocityList) {
        if (networkType == null || networkType.isBlank()) {
            throw new IllegalArgumentException("networkType must be non-blank");
        }
        if (gammaList == null || gammaList.length == 0) {
            throw new IllegalArgumentException("gammaList must be non-empty");
        }
        if (swapNumList == null || swapNumList.length == 0) {
            throw new IllegalArgumentException("swapNumList must be non-empty");
        }
        if (targetReciprocityList == null || targetReciprocityList.length == 0) {
            throw new IllegalArgumentException("targetReciprocityList must be non-empty");
        }
        if (!usesGamma(networkType) && gammaList.length > 1) {
            throw new IllegalArgumentException(networkType + " does not use gamma; gammaList must have one value");
        }
        if (!usesSwapNum(networkType) && swapNumList.length > 1) {
            throw new IllegalArgumentException(
                    networkType + " does not use swapNum; swapNumList must have one value");
        }

        Set<Double> gammas = new HashSet<>();
        for (double gamma : gammaList) {
            if (!Double.isFinite(gamma)) {
                throw new IllegalArgumentException("gammaList values must be finite");
            }
            if (usesGamma(networkType) && gamma <= 1.0) {
                throw new IllegalArgumentException("gammaList values must be > 1.0");
            }
            if (!gammas.add(gamma)) {
                throw new IllegalArgumentException("gammaList must not contain duplicates");
            }
        }

        Set<Integer> swaps = new HashSet<>();
        for (int swapNum : swapNumList) {
            if (!swaps.add(swapNum)) {
                throw new IllegalArgumentException("swapNumList must not contain duplicates");
            }
        }

        Set<Double> reciprocities = new HashSet<>();
        for (double target : targetReciprocityList) {
            if (!Double.isFinite(target) || target < 0.0 || target > 1.0) {
                throw new IllegalArgumentException("targetReciprocityList values must be in [0, 1]");
            }
            if (!reciprocities.add(target)) {
                throw new IllegalArgumentException("targetReciprocityList must not contain duplicates");
            }
        }
    }

    static long scaledAttempts(int edgeCount, int attemptsPerEdge, String parameterName) {
        if (edgeCount < 0) {
            throw new IllegalArgumentException("edgeCount must be non-negative");
        }
        if (attemptsPerEdge < 0) {
            throw new IllegalArgumentException(parameterName + "PerEdge must be non-negative");
        }
        return Math.multiplyExact((long) edgeCount, attemptsPerEdge);
    }

    private static void validateSimulationConfig(SimulationConfig config) {
        if (config.batchSize < 1 || config.itrs < 1) {
            throw new IllegalArgumentException("batchSize and itrs must be positive");
        }
        if (config.maxIncreaseAttemptsPerEdge < 0 || config.neutralSwapAttemptsPerEdge < 0) {
            throw new IllegalArgumentException("reciprocity attempt multipliers must be non-negative");
        }
        validateNonEmpty(config.rho0List, "rho0List");
        validateNonEmpty(config.lambdaDirectedList, "lambdaDirectedList");
        validateNonEmpty(config.lambdaNondirectedList, "lambdaNondirectedList");
        for (double rho0 : config.rho0List) {
            if (!Double.isFinite(rho0) || rho0 < 0.0 || rho0 > 1.0) {
                throw new IllegalArgumentException("rho0List values must be in [0, 1]");
            }
        }
    }

    private static void validateNonEmpty(double[] values, String name) {
        if (values == null || values.length == 0) {
            throw new IllegalArgumentException(name + " must be non-empty");
        }
        for (double value : values) {
            if (!Double.isFinite(value)) {
                throw new IllegalArgumentException(name + " values must be finite");
            }
        }
    }

    private static Thread createTotalProgressRenderer(
            AtomicLong done, long totalTasks, AtomicBoolean running) {
        return new Thread(() -> {
            long lastPrintedDone = -1L;
            while (running.get() && done.get() < totalTasks) {
                long completed = done.get();
                if (completed != lastPrintedDone) {
                    synchronized (System.out) {
                        renderTotalProgressBar(completed, totalTasks);
                    }
                    lastPrintedDone = completed;
                }
                try {
                    Thread.sleep(PROGRESS_UPDATE_INTERVAL_MS);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    break;
                }
            }
            synchronized (System.out) {
                renderTotalProgressBar(done.get(), totalTasks);
                System.out.println();
            }
        }, "network-sweep-progress-renderer");
    }

    private static void renderTotalProgressBar(long done, long total) {
        long boundedDone = Math.min(done, total);
        int filled = total == 0 ? PROGRESS_BAR_LENGTH
                : (int) (boundedDone * PROGRESS_BAR_LENGTH / total);
        int percent = total == 0 ? 100 : (int) (boundedDone * 100 / total);
        String bar = "#".repeat(filled) + "-".repeat(PROGRESS_BAR_LENGTH - filled);
        System.out.print("\033[2K\rProgress [%s] %3d%% (%d/%d)"
                .formatted(bar, percent, boundedDone, total));
        System.out.flush();
    }

    /** 実験設定。必要に応じてソース上で変更する。 */
    private static final class SimulationConfig {
        final String networkType = "PowPow";
        final String optionPath = "network-sweep";

        final int N = 500_000;
        final Integer kdAve = 5;
        final Integer kdMin = 5;
        final Integer kdMax = (int) Math.sqrt(N);
        final Integer kInMin = 5;
        final Integer kInMax = (int) Math.sqrt(N);
        final Integer kOutMin = 5;
        final Integer kOutMax = (int) Math.sqrt(N);
        final Integer kuMin = 5;
        final Integer kuMax = (int) Math.sqrt(N);
        final Double kuAve = 0.0;
        final Integer m0 = 5;
        final Integer m = 3;
        final Double gammaIn = 2.5;
        final Double gammaOut = 2.5;
        final Double corrA = 0.0;

        final double[] gammaList = { 2.5 };
        final int[] swapNumList = { 0 };
        final double[] targetReciprocityList = { 0.00, 0.0001, 0.001, 0.01 };
        final int maxIncreaseAttemptsPerEdge = 200;
        final int neutralSwapAttemptsPerEdge = 1;

        final boolean isFinal = true;
        final int batchSize = 10;
        final int itrs = 10;
        final double mu = 1.0;
        final double tMax = 200.0;
        final double[] lambdaDirectedList = { 2.0 };
        final double[] lambdaNondirectedList = { 0.0 };
        final double rho0Min = 0.0;
        final double rho0Max = 0.1;
        final double rho0Step = 0.002;
        final double[] rho0List = ArrayUtils.arange(rho0Min, rho0Max, rho0Step);
        final int threshold = 3;
        final boolean useGillespie = false;
    }
}
