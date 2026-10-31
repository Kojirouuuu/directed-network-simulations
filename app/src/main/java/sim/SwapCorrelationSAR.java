package sim;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.UUID;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.stream.IntStream;

import sim.network.DirectedGraph;
import sim.network.topology.PowPow;
import sim.simulation.SARResult;
import sim.simulation.SARSimulator;
import sim.utils.ArrayUtils;
import sim.utils.SwitchUtils;

/** Isolated, reproducible PowPow in/out-degree correlation experiment. */
public final class SwapCorrelationSAR {
    private static final int THRESHOLD = 3;
    private static final int K_MIN = 5;
    private static final int K_MAX = 707;
    private static final double GAMMA = 2.5;
    private static final double RHO_0 = 0.05;
    private static final double MU = 1.0;
    private static final double T_MAX = 200.0;
    private static final long GRAPH_SEED = 42L;
    private static final long NODE_SEED = 2007L;
    private static final long SIM_SEED = 12345L;

    private SwapCorrelationSAR() {
    }

    record Config(int n, int swap, int batches, int itrs, int intervals,
            double lambdaMax, int parallelism, String phase) {
        static Config parse(String[] args) {
            Map<String, String> options = new HashMap<>();
            for (String arg : args) {
                if (!arg.startsWith("--") || !arg.contains("=")) {
                    throw new IllegalArgumentException("Expected --key=value: " + arg);
                }
                String[] pair = arg.substring(2).split("=", 2);
                if (!java.util.Set.of("n", "swap", "batches", "itrs", "intervals", "lambda-max",
                        "parallel", "phase")
                        .contains(pair[0]) || options.putIfAbsent(pair[0], pair[1]) != null) {
                    throw new IllegalArgumentException("Unknown or duplicate option: " + pair[0]);
                }
            }
            if (!options.containsKey("swap")) {
                throw new IllegalArgumentException("--swap is required");
            }
            int n = Integer.parseInt(options.getOrDefault("n", "500000"));
            int swap = Integer.parseInt(options.get("swap"));
            int batches = Integer.parseInt(options.getOrDefault("batches", "10"));
            int itrs = Integer.parseInt(options.getOrDefault("itrs", "20"));
            int intervals = Integer.parseInt(options.getOrDefault("intervals", "100"));
            double lambdaMax = Double.parseDouble(options.getOrDefault("lambda-max", "2.0"));
            int parallelism = Integer.parseInt(options.getOrDefault("parallel",
                    Integer.toString(Math.min(batches, Runtime.getRuntime().availableProcessors()))));
            String phase = options.getOrDefault("phase", "final");
            if (n < K_MIN || batches < 1 || itrs < 1 || intervals < 1
                    || !Double.isFinite(lambdaMax) || lambdaMax <= 0.0
                    || parallelism < 1 || parallelism > batches || swap == Integer.MIN_VALUE
                    || !phase.matches("[A-Za-z0-9_-]+")) {
                throw new IllegalArgumentException("Invalid experiment configuration");
            }
            return new Config(n, swap, batches, itrs, intervals, lambdaMax, parallelism, phase);
        }

        double lambdaAt(int index) {
            return lambdaMax * index / intervals;
        }

        long plannedCount() {
            return (long) batches * itrs * (intervals + 1);
        }
    }

    public static void main(String[] args) throws Exception {
        Config config = Config.parse(args);
        Path networkPath = SwitchUtils.buildNetworkPath("PowPow", config.n,
                null, null, null, null, null, null,
                K_MIN, K_MAX, null, null, null, null,
                GAMMA, config.swap, null, null, null);
        Path outputDir = SwitchUtils.buildSimulationOutputDir(
                "swap-negative/" + config.phase, THRESHOLD)
                .resolve(networkPath)
                .resolve("randomization=none")
                .resolve("run=" + System.currentTimeMillis() + "-"
                        + UUID.randomUUID().toString().substring(0, 8));

        Map<String, Object> network = new LinkedHashMap<>();
        network.put("name", "PowPow");
        network.put("configuredN", config.n);
        network.put("gamma", GAMMA);
        network.put("kdMin", K_MIN);
        network.put("kdMax", K_MAX);
        network.put("swapNum", config.swap);
        network.put("randomizationMode", "none");
        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("threshold", THRESHOLD);
        parameters.put("rho0Values", new double[] { RHO_0 });
        parameters.put("lambdaDirectedValues", IntStream.rangeClosed(0, config.intervals)
                .mapToDouble(config::lambdaAt).toArray());
        parameters.put("lambdaNondirectedValues", new double[] { 0.0 });
        parameters.put("mu", MU);
        parameters.put("tMax", T_MAX);
        parameters.put("graphSeed", GRAPH_SEED);
        parameters.put("nodeSeed", NODE_SEED);
        parameters.put("simulationSeed", SIM_SEED);
        Map<String, Object> execution = new LinkedHashMap<>();
        execution.put("phase", config.phase);
        execution.put("batchSize", config.batches);
        execution.put("iterationsPerBatch", config.itrs);
        execution.put("parallelism", config.parallelism);
        execution.put("finalStateOnly", true);
        execution.put("useGillespie", false);

        SARRunSummary summary = new SARRunSummary(config.plannedCount(), network, parameters, execution);
        Thread shutdownHook = new Thread(summary::finishInterruptedBestEffort,
                "swap-sar-summary-shutdown-hook");
        Runtime.getRuntime().addShutdownHook(shutdownHook);
        Throwable failure = null;
        System.out.printf(Locale.ROOT,
                "PowPow N=%d swap=%d batches=%d itrs=%d lambda points=%d total=%d output=%s%n",
                config.n, config.swap, config.batches, config.itrs,
                config.intervals + 1, config.plannedCount(), outputDir);
        try {
            try (ForkJoinPool pool = new ForkJoinPool(config.parallelism)) {
                Future<?> future = pool.submit(() -> IntStream.range(0, config.batches).parallel()
                        .forEach(batch -> runBatch(config, batch, outputDir, summary)));
                future.get();
            }
        } catch (Exception | Error error) {
            failure = error;
            summary.recordUnhandledErrorIfAbsent(error);
            throw error;
        } finally {
            try {
                if (failure == null) {
                    summary.finishSucceeded();
                } else {
                    summary.finishFailed();
                }
                System.out.println("Run summary: " + summary.summaryPath());
            } finally {
                Runtime.getRuntime().removeShutdownHook(shutdownHook);
            }
        }
    }

    private static void runBatch(Config config, int batch, Path outputDir, SARRunSummary summary) {
        summary.recordBatchStarted(batch, Instant.now());
        String stage = "graph-generation";
        SARRunSummary.ErrorContext context = SARRunSummary.ErrorContext.forBatch(batch);
        try {
            DirectedGraph graph = PowPow.generate("PowPow", config.n, K_MIN, K_MAX,
                    GAMMA, config.swap, GRAPH_SEED + batch);
            summary.registerBatchGraph(batch, graph.n, graph.m, outputDir);
            stage = "degree-pair-output";
            writeDegreePairs(outputDir.resolve(String.format(Locale.ROOT,
                    "degree_pairs_%02d.csv", batch)), graph);

            Path resultPath = outputDir.resolve(String.format(Locale.ROOT, "results_%02d.csv", batch));
            summary.recordResultPath(batch, resultPath);
            int[] thresholds = new int[graph.n];
            Arrays.fill(thresholds, THRESHOLD);
            for (int itr = 0; itr < config.itrs; itr++) {
                int resultItr = batch * config.itrs + itr;
                int[] initial = initialNodes(graph.n, resultItr);
                for (int li = 0; li <= config.intervals; li++) {
                    double lambda = config.lambdaAt(li);
                    context = SARRunSummary.ErrorContext.forSimulation(
                            batch, itr, RHO_0, lambda, 0.0, MU);
                    stage = "simulation";
                    SARResult result = SARSimulator.simulate(graph, lambda, 0.0, MU,
                            T_MAX, thresholds, initial, SIM_SEED + resultItr);
                    stage = "result-output";
                    result.writeFinalStateCsv(resultPath, resultItr, RHO_0, lambda, 0.0, MU, true);
                    summary.recordSimulationCompleted();
                }
            }
            summary.recordBatchCompleted(batch, Instant.now());
            synchronized (System.out) {
                System.out.printf(Locale.ROOT, "Completed swap=%d batch=%d/%d%n",
                        config.swap, batch + 1, config.batches);
            }
        } catch (IOException error) {
            summary.recordError(stage, context, error);
            summary.recordBatchFailed(batch, Instant.now());
            throw new UncheckedIOException(error);
        } catch (RuntimeException | Error error) {
            summary.recordError(stage, context, error);
            summary.recordBatchFailed(batch, Instant.now());
            throw error;
        }
    }

    private static int[] initialNodes(int n, int resultItr) {
        int[] nodes = IntStream.range(0, n).toArray();
        ArrayUtils.shuffle(nodes, NODE_SEED + resultItr);
        return Arrays.copyOf(nodes, Math.max(1, (int) (n * RHO_0)));
    }

    private static void writeDegreePairs(Path path, DirectedGraph graph) throws IOException {
        Map<Long, Integer> counts = new TreeMap<>();
        for (int u = 0; u < graph.n; u++) {
            DirectedGraph.IntRange in = graph.inNeighborRange(u);
            DirectedGraph.IntRange out = graph.outNeighborRange(u);
            int ki = in.end - in.start;
            int ko = out.end - out.start;
            long key = ((long) ki << 32) | Integer.toUnsignedLong(ko);
            counts.merge(key, 1, Integer::sum);
        }
        Files.createDirectories(path.getParent());
        try (BufferedWriter writer = Files.newBufferedWriter(path)) {
            writer.write("ki,ko,count\n");
            for (Map.Entry<Long, Integer> entry : counts.entrySet()) {
                int ki = (int) (entry.getKey() >>> 32);
                int ko = (int) (entry.getKey() & 0xffffffffL);
                writer.write(ki + "," + ko + "," + entry.getValue() + "\n");
            }
        }
    }
}
