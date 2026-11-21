package sim;

import sim.utils.ProgressBar;
import sim.network.DirectedGraph;
import sim.simulation.BondPercolationResult;
import sim.simulation.BondPercolationSimulator;
import sim.utils.ArrayUtils;
import sim.utils.SwitchUtils;
import sim.SAR.RandomizationMode;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Files;
import java.util.Random;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;

/**
 * ボンドパーコレーションシミュレーションのメインクラス。
 * 並列処理によるバッチシミュレーションを実行する。
 */
public class BondPercolation {
    private static final int PROGRESS_UPDATE_INTERVAL_MS = 100;

    private static final long SIM_BASE_SEED = 12345L;
    private static final long GRAPH_BASE_SEED = 42L;
    private static final long RANDOMIZATION_BASE_SEED = 4242L;

    public static void main(String[] args) throws Exception {
        SimulationConfig config = new SimulationConfig();

        final long totalTasks = (long) config.batchSize * config.itrs * config.pList.length;

        System.out.println("Total tasks: " + totalTasks);
        System.out.println(config.networkType + ": N=" + config.N + ", itrs=" + config.itrs
                + ", pList.length=" + config.pList.length);

        int parallelism = Runtime.getRuntime().availableProcessors();
        System.out.println("Parallelism: " + parallelism + " (available processors)");

        AtomicLong done = new AtomicLong(0);
        AtomicBoolean running = new AtomicBoolean(true);

        Thread renderer = createTotalProgressRenderer(done, totalTasks, running);
        renderer.start();

        try (ForkJoinPool pool = new ForkJoinPool(parallelism)) {
            DirectedGraph sharedGraph = config.randomizationMode == RandomizationMode.JOINT_DEGREE_CM
                    ? SARGraphPreparation.loadGraph(networkParameters(config), config.loadFromEdgeList,
                            config.randomizationMode, GRAPH_BASE_SEED, 0)
                    : null;
            Future<?> future = pool.submit(() -> IntStream.range(0, config.batchSize).parallel()
                    .forEach(batchIndex -> processBatch(batchIndex, config, sharedGraph, done)));
            future.get();
        } finally {
            running.set(false);
            renderer.join();
        }

        System.out.println("All tasks completed");
    }

    private static SwitchUtils.GraphGeneratorParams networkParameters(SimulationConfig config) {
        return SwitchUtils.GraphGeneratorParams.builder(config.networkType, config.N)
                .kdMin(config.kdMin)
                .kdMax(config.kdMax)
                .kInMin(config.kInMin)
                .kInMax(config.kInMax)
                .kOutMin(config.kOutMin)
                .kOutMax(config.kOutMax)
                .kuAve(config.kuAve)
                .gamma(config.gamma)
                .kuMin(config.kuMin)
                .kuMax(config.kuMax)
                .m0(config.m0)
                .m(config.m)
                .swapNum(config.swapNum)
                .gammaIn(config.gammaIn)
                .gammaOut(config.gammaOut)
                .corrA(config.corrA).k1(config.k1).k2(config.k2).p(config.p).r(config.r).build();
    }

    private static void processBatch(int batchIndex, SimulationConfig config,
            DirectedGraph sharedGraph, AtomicLong done) {
        DirectedGraph g = sharedGraph == null
                ? SARGraphPreparation.loadGraph(networkParameters(config), config.loadFromEdgeList,
                        config.randomizationMode, GRAPH_BASE_SEED, batchIndex)
                : sharedGraph;
        g = SARGraphPreparation.applyRandomization(g, config.randomizationMode,
                RANDOMIZATION_BASE_SEED + batchIndex, config.loadFromEdgeList, config.sizeMultiplier);

        OutputPaths paths = prepareOutputPaths(g, batchIndex, config);

        BondPercolationSimulator sim = new BondPercolationSimulator(g);

        for (int itr = 0; itr < config.itrs; itr++) {
            for (int pi = 0; pi < config.pList.length; pi++) {
                double p = config.pList[pi];
                long simSeed = SIM_BASE_SEED
                        + (long) batchIndex * config.itrs * config.pList.length
                        + (long) itr * config.pList.length
                        + pi;
                BondPercolationResult res = sim.run(p, new Random(simSeed));
                try {
                    res.writeFinalStateCsv(paths.results(), batchIndex * config.itrs + itr, p, true);
                    res.writeComponentSizeDistributionCsv(paths.componentSizes(),
                            batchIndex * config.itrs + itr, p, true);
                } catch (IOException e) {
                    throw new RuntimeException("CSV output error (batch=" + batchIndex
                            + ", itr=" + itr + ", p=" + p + "): " + e.getMessage(), e);
                }
                done.incrementAndGet();
            }
        }
    }

    private static OutputPaths prepareOutputPaths(DirectedGraph g, int batchIndex, SimulationConfig config) {
        Path outputDir = SwitchUtils.buildPercolationOutputDir(config.optionPath);
        Path networkPath = SARGraphPreparation.buildNetworkOutputPath(
                networkParameters(config), g.n, config.randomizationMode, config.sizeMultiplier);
        Path basePath = outputDir.resolve(networkPath);
        return resolveOutputPaths(basePath, batchIndex);
    }

    record OutputPaths(Path results, Path componentSizes) {
    }

    static OutputPaths resolveOutputPaths(Path directory, int batchIndex) {
        String idx = String.format(java.util.Locale.ROOT, "%02d", batchIndex);
        for (int index = 0;; index++) {
            String suffix = index == 0 ? "" : " (" + index + ")";
            Path results = directory.resolve("results_" + idx + suffix + ".csv");
            Path sizes = directory.resolve("component_sizes_" + idx + suffix + ".csv");
            if (!Files.exists(results) && !Files.exists(sizes)) {
                return new OutputPaths(results, sizes);
            }
        }
    }

    private static Thread createTotalProgressRenderer(AtomicLong done, long totalTasks, AtomicBoolean running) {
        return new Thread(() -> {
            long lastPrintedDone = 0;
            while (running.get()) {
                long d = done.get();
                if (d != lastPrintedDone) {
                    synchronized (System.out) {
                        renderTotalProgressBar(d, totalTasks);
                    }
                    lastPrintedDone = d;
                }
                if (d >= totalTasks)
                    break;
                try {
                    Thread.sleep(PROGRESS_UPDATE_INTERVAL_MS);
                } catch (InterruptedException e) {
                    running.set(false);
                    break;
                }
            }
            synchronized (System.out) {
                renderTotalProgressBar(totalTasks, totalTasks);
                System.out.println();
            }
        }, "progress-renderer");
    }

    private static void renderTotalProgressBar(long done, long total) {
        ProgressBar.render(done, total);
    }

    private static class SimulationConfig {
        final String networkType = "ER";
        final String optionPath = "play";
        final int N = (int) Math.pow(2, 16);
        final int kdMin = 2;
        final int kdMax = 1000;
        final int kInMin = 2;
        final int kInMax = (int) Math.sqrt(N);
        final int kOutMin = 2;
        final int kOutMax = (int) Math.sqrt(N);
        final double kuAve = 4.0;
        final int kuMin = 2;
        final int kuMax = (int) Math.sqrt(N);
        final int m0 = 6;
        final int m = 6;
        final double gamma = 3.3;
        final int swapNum = 0;

        // SchwartzDirectedSF 用パラメータ（他のネットワークタイプでは未使用）
        final Double gammaIn = 3.3; // λ_in
        final Double gammaOut = 3.3; // λ_out
        final Double corrA = null; // 相関確率 A ∈ [0, 1]

        // TwoPointDirectedCM: set all four when selecting this topology.
        final Integer k1 = null, k2 = null;
        final Double p = null, r = null;

        final boolean loadFromEdgeList = false;
        final RandomizationMode randomizationMode = RandomizationMode.NONE;
        final int sizeMultiplier = 1; // JOINT_DEGREE_CM 用の頂点数倍率

        final int batchSize = 10;
        final int itrs = 1000;
        final double pMin = 0.0;
        final double pMax = 1.0;
        final double pStep = 0.005;
        final double[] pList = ArrayUtils.arange(pMin, pMax, pStep);
    }
}
