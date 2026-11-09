package sim;

import java.nio.file.*;
import java.time.Instant;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;
import sim.SAR.RandomizationMode;
import sim.network.DirectedGraph;
import sim.simulation.InitialAdopters;
import sim.simulation.SPIRGillespieSimulator.ModelMode;
import sim.simulation.SPIRResult;
import sim.simulation.SPIRSimulator;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;
import sim.utils.ProgressBar;
import sim.utils.SwitchUtils;

/**
 * Parallel SPIR/SIR parameter sweeps. Edit SimulationConfig to configure an
 * experiment.
 */
public final class SPIR {
    private static final long GRAPH_SEED = 42, RANDOMIZATION_SEED = 4242, INITIAL_SEED = 2007,
            SIM_SEED = 12345, OBSERVER_SEED = 424242;

    record Population(int[] initial, int[] observers) {
    }

    static Population samplePopulation(int n, double rho0, double phi, long initialSeed, long observerSeed) {
        if (n <= 0 || !Double.isFinite(phi) || phi < 0 || phi > 1) {
            throw new IllegalArgumentException("positive n and phi in [0,1] required");
        }
        int[] nodes = IntStream.range(0, n).toArray();
        ArrayUtils.shuffle(nodes, observerSeed);
        int[] observers = Arrays.copyOf(nodes, (int) (n * phi));
        boolean[] mask = new boolean[n];
        for (int node : observers)
            mask[node] = true;
        return new Population(InitialAdopters.sample(n, rho0, mask, initialSeed), observers);
    }

    record Rates(double lambda1, double lambda2, double mu1, double mu2) {
    }

    static List<Rates> rates(ModelMode mode, double[] l1, double[] l2, double[] m1, double[] m2) {
        return rates(mode, l1, l2, m1, m2, false);
    }

    static List<Rates> rates(ModelMode mode, double[] l1, double[] l2, double[] m1, double[] m2,
            boolean equalLambda) {
        if (mode == null)
            throw new IllegalArgumentException("mode is required");
        for (double[] values : new double[][] { l1, l2, m1, m2 }) {
            if (values == null || values.length == 0)
                throw new IllegalArgumentException("rate lists must not be empty");
            for (double v : values)
                if (!Double.isFinite(v) || v < 0)
                    throw new IllegalArgumentException("rates must be finite and non-negative");
        }
        List<Rates> values = new ArrayList<>();
        for (double a : mode == ModelMode.SIR || equalLambda ? new double[] { 0 } : l1)
            for (double b : l2)
                for (double c : mode == ModelMode.SIR ? new double[] { 0 } : m1)
                    for (double d : m2)
                        values.add(new Rates(mode == ModelMode.SPIR && equalLambda ? b : a, b, c, d));
        return values;
    }

    static Path outputDirectory(DirectedGraph graph, SimulationConfig config) {
        return config.outputRoot
                .resolve(Path.of(
                        "spir", config.optionPath, "model=" + config.mode.name().toLowerCase(Locale.ROOT)))
                .resolve(SARGraphPreparation.appendRandomizationPath(
                        SwitchUtils.buildNetworkPath(config.network.withN(graph.n)), config.randomizationMode,
                        config.sizeMultiplier))
                .resolve("phi=" + Double.toString(config.phi));
    }

    public static void main(String[] args) throws Exception {
        run(new SimulationConfig());
    }

    static void run(SimulationConfig config) throws Exception {
        if (!Double.isFinite(config.phi) || config.phi < 0 || config.phi > 1)
            throw new IllegalArgumentException("phi must be in [0,1]");
        if (config.batchSize <= 0 || config.itrs <= 0 || config.rho0List.length == 0)
            throw new IllegalArgumentException("batch, iteration and rho0 lists must be positive/nonempty");
        for (double rho : config.rho0List)
            if (!Double.isFinite(rho) || rho < 0 || rho > 1)
                throw new IllegalArgumentException("rho0 must be in [0,1]");
        List<Rates> rates = rates(config.mode, config.lambda1List, config.lambda2List,
                config.mu1List, config.mu2List, config.equalLambda);
        long total = Math.multiplyExact(Math.multiplyExact((long) config.batchSize, config.itrs),
                Math.multiplyExact((long) config.rho0List.length, rates.size()));
        Map<String, Object> network = new LinkedHashMap<>();
        network.put("parameters", config.network);
        network.put("randomizationMode", config.randomizationMode);
        network.put("sizeMultiplier", config.sizeMultiplier);
        network.put("loadFromEdgeList", config.loadFromEdgeList);
        network.put("writeEdgeList", config.writeEdgeList);
        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("rho0Values", config.rho0List);
        parameters.put("phi", config.phi);
        parameters.put("engine", config.engine);
        parameters.put("equalLambda", config.equalLambda);
        parameters.put("lambda1Values", config.lambda1List);
        parameters.put("lambda2Values", config.lambda2List);
        parameters.put("mu1Values", config.mu1List);
        parameters.put("mu2Values", config.mu2List);
        parameters.put("tMax", config.tMax);
        parameters.put("dt", config.dt);
        parameters.put("seeds",
                Map.of("graph", GRAPH_SEED, "randomization", RANDOMIZATION_SEED, "initial", INITIAL_SEED,
                        "simulation", SIM_SEED, "observer", OBSERVER_SEED));
        Map<String, Object> execution = Map.of("model", config.mode, "optionPath", config.optionPath,
                "batchSize", config.batchSize, "iterationsPerBatch", config.itrs, "finalStateOnly",
                config.isFinal, "parallelism", Runtime.getRuntime().availableProcessors());
        SPIRRunSummary summary = new SPIRRunSummary(total, network, parameters, execution);
        Thread shutdown = new Thread(summary::finishInterruptedBestEffort, "spir-summary-shutdown");
        Runtime.getRuntime().addShutdownHook(shutdown);
        AtomicLong done = new AtomicLong();
        Throwable failure = null;
        try {
            DirectedGraph shared = config.randomizationMode == RandomizationMode.JOINT_DEGREE_CM
                    ? SARGraphPreparation.loadGraph(config.network, config.loadFromEdgeList,
                            config.randomizationMode, GRAPH_SEED, 0)
                    : null;
            try (ForkJoinPool pool = new ForkJoinPool()) {
                pool.submit(() -> IntStream.range(0, config.batchSize).parallel().forEach(batch -> {
                    summary.recordBatchStarted(batch, Instant.now());
                    SPIRRunSummary.ErrorContext context = SPIRRunSummary.ErrorContext.forBatch(batch);
                    String stage = "graph-preparation";
                    try {
                        DirectedGraph graph = shared == null
                                ? SARGraphPreparation.loadGraph(config.network, config.loadFromEdgeList,
                                        config.randomizationMode, GRAPH_SEED, batch)
                                : shared;
                        graph = SARGraphPreparation.applyRandomization(graph, config.randomizationMode,
                                RANDOMIZATION_SEED + batch, config.loadFromEdgeList,
                                config.sizeMultiplier);
                        var prepared = new SPIRSimulator.PreparedGraph(graph);
                        Path directory = outputDirectory(graph, config);
                        summary.registerBatchGraph(batch, graph.n, graph.m, directory);
                        if (config.writeEdgeList && !config.loadFromEdgeList) {
                            Path edgePath = Path.of("out", "edgelist")
                                    .resolve(SARGraphPreparation.appendRandomizationPath(
                                            SwitchUtils.buildNetworkPath(
                                                    config.network.withN(graph.n)),
                                            config.randomizationMode, config.sizeMultiplier))
                                    .resolve(batch + ".csv");
                            graph.writeEdgeList(edgePath);
                        }
                        Path csv = PathsEx.resolveIndexed(
                                directory.resolve(String.format(Locale.ROOT, "results_%02d.csv", batch)));
                        summary.recordResultPath(batch, csv);
                        try (var csvBatch = SPIRResult.openCsv(csv)) {
                            for (int itr = 0; itr < config.itrs; itr++) {
                                long index = (long) batch * config.itrs + itr;
                                for (double rho : config.rho0List) {
                                    Population population = samplePopulation(graph.n, rho, config.phi,
                                            INITIAL_SEED + index, OBSERVER_SEED + index);
                                    int[] initial = population.initial();
                                    for (Rates r : rates) {
                                        context = SPIRRunSummary.ErrorContext.forSimulation(
                                                batch, itr, rho, r.lambda1, r.lambda2, r.mu1, r.mu2);
                                        stage = "simulation";
                                        var result = SPIRSimulator.simulate(prepared,
                                                SPIRSimulator.ModelMode.valueOf(config.mode.name()),
                                                config.engine,
                                                config.isFinal
                                                        ? SPIRSimulator.RecordingMode.FINAL_ONLY
                                                        : SPIRSimulator.RecordingMode.FIXED_INTERVAL,
                                                r.lambda1, r.lambda2, r.mu1, r.mu2, config.tMax,
                                                config.dt, initial, population.observers(),
                                                SIM_SEED + index);
                                        stage = "csv-output";
                                        csvBatch.write(result, Math.toIntExact(index), rho, config.phi,
                                                config.mode, r.lambda1, r.lambda2, r.mu1, r.mu2,
                                                config.isFinal);
                                        summary.recordSimulationCompleted();
                                        synchronized (System.out) {
                                            ProgressBar.render(done.incrementAndGet(), total);
                                        }
                                    }
                                }
                            }
                        }
                        summary.recordBatchCompleted(batch, Instant.now());
                    } catch (Exception | Error e) {
                        summary.recordError(stage, context, e);
                        if (context.iteration() != null)
                            summary.recordSimulationFailed();
                        summary.recordBatchFailed(batch, Instant.now());
                        throw new RuntimeException(e);
                    }
                })).get();
            }
            summary.finishSucceeded();
        } catch (Exception | Error e) {
            failure = e;
            summary.recordUnhandledErrorIfAbsent(e);
            throw e;
        } finally {
            try {
                if (failure != null)
                    summary.finishFailed();
            } finally {
                try {
                    Runtime.getRuntime().removeShutdownHook(shutdown);
                } catch (IllegalStateException ignored) {
                }
                System.out.println("\nRun summary: " + summary.summaryPath());
            }
        }
    }

    /** シミュレーション設定。実験条件はここで変更する。 */
    static class SimulationConfig {
        Path outputRoot = Path.of("out");
        ModelMode mode = ModelMode.SIR;
        SPIRSimulator.Engine engine = SPIRSimulator.Engine.GILLESPIE;
        String optionPath = "o-spir-test";
        SwitchUtils.GraphGeneratorParams network = SwitchUtils.GraphGeneratorParams.builder("ER", 100_000).kuAve(4.0)
                .build();
        boolean loadFromEdgeList = false, writeEdgeList = false, isFinal = true;
        RandomizationMode randomizationMode = RandomizationMode.NONE;
        int sizeMultiplier = 1, batchSize = 10, itrs = 10;
        double tMax = 100, dt = 0.1;
        double phi = 0.3; // オブサーバーの割合 [0,1]
        double[] rho0List = { 0.001 }; // 初期P: 100,000頂点中100人
        boolean equalLambda = true; // true: lambda2List の各値をlambda1・lambda2共通に使用
        double[] lambda1List = { 0.0 };
        // double[] lambda1List = ArrayUtils.arange(0.0, 2.0, 0.02);
        double[] lambda2List = ArrayUtils.arange(0.0, 1.0, 0.01);
        double[] mu1List = { 1.0 }, mu2List = { 1.0 };
    }
}
