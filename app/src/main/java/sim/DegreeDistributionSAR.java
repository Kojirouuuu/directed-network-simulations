package sim;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;

import sim.network.DirectedGraph;
import sim.network.topology.DegreeDistributionNetworks;
import sim.network.topology.DegreeDistributionNetworks.Variant;
import sim.simulation.SARGillespieSimulator;
import sim.simulation.SARResult;
import sim.simulation.SARSimulator;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;
import sim.utils.SwitchUtils;

/**
 * 入次数・出次数を Pow/Poi の4通りに組み合わせて比較するSAR実験。
 */
public final class DegreeDistributionSAR {

    private static final int NETWORK_VARIANT_COUNT = 4;
    private static final int PROGRESS_BAR_LENGTH = 100;
    private static final int PROGRESS_UPDATE_INTERVAL_MS = 100;

    private static final long RNG_BASE_SEED = 7L;
    private static final long SIM_BASE_SEED = 12345L;
    private static final long GRAPH_BASE_SEED = 42L;
    private static final long SEED_OFFSET_NODES = 2000L;

    private DegreeDistributionSAR() {
    }

    public static void main(String[] args) throws Exception {
        SimulationConfig config = new SimulationConfig();
        long totalTasks = (long) config.batchSize * NETWORK_VARIANT_COUNT * config.itrs
                * config.rho0List.length * config.lambdaDirectedList.length
                * config.lambdaNondirectedList.length;

        System.out.println("Total tasks: " + totalTasks);
        System.out.println(config.networkType + ": N=" + config.N + ", itrs=" + config.itrs);

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
        Path resultsPath = prepareOutputPath(batchIndex, config);
        DegreeDistributionNetworks.forEachVariant(
                config.networkType, config.N, config.kdMin, config.kdMax, config.gamma,
                GRAPH_BASE_SEED + batchIndex,
                variant -> processVariant(variant, batchIndex, config, resultsPath, done));
    }

    private static void processVariant(Variant variant, int batchIndex, SimulationConfig config,
            Path resultsPath, AtomicLong done) {
        DirectedGraph graph = variant.graph();
        for (int itr = 0; itr < config.itrs; itr++) {
            for (double rho0 : config.rho0List) {
                for (double lambdaDirected : config.lambdaDirectedList) {
                    for (double lambdaNondirected : config.lambdaNondirectedList) {
                        int[] thresholdList = new int[graph.n];
                        Arrays.fill(thresholdList, config.threshold);

                        runSimulation(variant, config, lambdaDirected, lambdaNondirected, rho0,
                                thresholdList, batchIndex, itr, resultsPath);
                        done.incrementAndGet();
                    }
                }
            }
        }
    }

    private static Path prepareOutputPath(int batchIndex, SimulationConfig config) {
        String index = String.format("%02d", batchIndex);
        Path outputDir = SwitchUtils.buildSimulationOutputDir(config.optionPath, config.threshold);
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, config.N,
                null, null,
                null, null, null, null,
                config.kdMin, config.kdMax, null, null, null, null,
                config.gamma, null,
                null, null, null);
        return PathsEx.resolveIndexed(
                outputDir.resolve(networkPath).resolve(String.format("results_%s.csv", index)));
    }

    private static void runSimulation(Variant variant, SimulationConfig config,
            double lambdaDirected, double lambdaNondirected, double rho0, int[] thresholdList,
            int batchIndex, int itr, Path resultsPath) {
        DirectedGraph graph = variant.graph();
        int initialInfectedNum = Math.max(1, (int) (graph.n * rho0));

        int[] nodes = new int[graph.n];
        for (int i = 0; i < graph.n; i++) {
            nodes[i] = i;
        }
        long nodeSeed = RNG_BASE_SEED + (long) batchIndex * config.itrs + itr + SEED_OFFSET_NODES;
        ArrayUtils.shuffle(nodes, nodeSeed);
        int[] init = Arrays.copyOfRange(nodes, 0, initialInfectedNum);

        long simSeed = SIM_BASE_SEED + (long) batchIndex * config.itrs + itr;
        SARResult result = config.useGillespie
                ? SARGillespieSimulator.simulate(
                        graph, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                        thresholdList, init, simSeed)
                : SARSimulator.simulate(
                        graph, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                        thresholdList, init, simSeed);

        int resultItr = batchIndex * config.itrs + itr;
        String inDistribution = variant.inDistribution().label();
        String outDistribution = variant.outDistribution().label();
        try {
            if (config.isFinal) {
                result.writeFinalStateCsv(resultsPath, resultItr, inDistribution, outDistribution,
                        rho0, lambdaDirected, lambdaNondirected, config.mu, true);
            } else {
                result.writeTimeSeriesCsv(resultsPath, resultItr, inDistribution, outDistribution,
                        rho0, lambdaDirected, lambdaNondirected, config.mu, true);
            }
        } catch (IOException e) {
            throw new RuntimeException("CSV output error (batch=" + batchIndex
                    + ", iteration=" + itr + ", in=" + inDistribution + ", out=" + outDistribution
                    + ", rho0=" + rho0 + ", lambdaDirected=" + lambdaDirected
                    + ", lambdaNondirected=" + lambdaNondirected + ", mu=" + config.mu + ")", e);
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
        }, "degree-distribution-progress-renderer");
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

    /** 実験設定。既存のSARと同様、必要に応じてソース上で変更する。 */
    private static final class SimulationConfig {
        final String networkType = "DegreeDistributionSAR";
        final String optionPath = "dis-pair";
        final int N = 500_000;
        final int kdMin = 5;
        final int kdMax = (int) Math.sqrt(N);
        final double gamma = 2.5;

        final boolean isFinal = true;
        final int batchSize = 10;
        final int itrs = 10;
        final double mu = 1.0;
        final double tMax = 200.0;
        final double[] lambdaDirectedList = { 2.0 };
        final double[] lambdaNondirectedList = { 0.0 };
        final double rho0Min = 0.0;
        final double rho0Max = 0.1;
        final double rho0Step = 0.001;
        final double[] rho0List = ArrayUtils.arange(rho0Min, rho0Max, rho0Step);
        final int threshold = 3;
        final boolean useGillespie = false;
    }
}
