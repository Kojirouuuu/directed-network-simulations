package sim;

import sim.network.DirectedGraph;
import sim.simulation.BondPercolationResult;
import sim.simulation.BondPercolationSimulator;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
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
    private static final int PROGRESS_BAR_LENGTH = 100;
    private static final int PROGRESS_UPDATE_INTERVAL_MS = 100;

    private static final long SIM_BASE_SEED = 12345L;
    private static final long GRAPH_BASE_SEED = 42L;

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
        DirectedGraph g;
        Path resultsPath;

        if (config.loadFromEdgeList) {
            Path networkPath = SwitchUtils.buildNetworkPath(
                    config.networkType, config.N,
                    null, config.kuAve,
                    config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                    config.kdMin, config.kdMax, config.kuMin, config.kuMax,
                    config.m0, config.m, config.gamma, config.swapNum);
            Path edgeListPath = Paths.get("out/edgelist").resolve(networkPath)
                    .resolve(String.format("%d.csv", batchIndex));
            try {
                g = DirectedGraph.loadFromEdgeList(config.networkType, edgeListPath);
            } catch (IOException e) {
                throw new RuntimeException("Failed to load edge list: " + edgeListPath, e);
            }
        } else {
            g = SwitchUtils.generateGraph(config.networkType, config.N,
                    null, config.kdMin, config.kdMax, config.kInMin, config.kInMax,
                    config.kOutMin, config.kOutMax,
                    config.kuMin, config.kuMax,
                    config.kuAve, config.gamma, config.m0, config.m, config.swapNum,
                    GRAPH_BASE_SEED + batchIndex);
        }

        resultsPath = prepareOutputPath(g, batchIndex, config);

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
                    res.writeFinalStateCsv(resultsPath, batchIndex * config.itrs + itr, p, true);
                } catch (IOException e) {
                    throw new RuntimeException("CSV output error (batch=" + batchIndex
                            + ", itr=" + itr + ", p=" + p + "): " + e.getMessage(), e);
                }
                done.incrementAndGet();
            }
        }
    }

    private static Path prepareOutputPath(DirectedGraph g, int batchIndex, SimulationConfig config) {
        String idx = String.format("%02d", batchIndex);
        Path outputDir = SwitchUtils.buildPercolationOutputDir(config.optionPath);
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, g.n,
                null, config.kuAve,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kdMin, config.kdMax, config.kuMin, config.kuMax,
                config.m0, config.m, config.gamma, config.swapNum);
        Path basePath = outputDir.resolve(networkPath);
        return PathsEx.resolveIndexed(basePath.resolve(String.format("results_%s.csv", idx)));
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
        int filled = (total == 0) ? PROGRESS_BAR_LENGTH
                : (int) Math.min(PROGRESS_BAR_LENGTH, (done * PROGRESS_BAR_LENGTH) / total);
        int percent = (total == 0) ? 100 : (int) Math.min(100, (done * 100) / total);
        String bar = "#".repeat(filled) + "-".repeat(PROGRESS_BAR_LENGTH - filled);
        System.out.print("\033[2K\rProgress [%s] %3d%% (%d/%d)".formatted(bar, percent, done, total));
        System.out.flush();
    }

    private static class SimulationConfig {
        final String networkType = "PowPow";
        final String optionPath = "test";
        final int N = 1_000_000;
        final int kdMin = 1;
        final int kdMax = 1000;
        final int kInMin = 1;
        final int kInMax = (int) Math.sqrt(N);
        final int kOutMin = 1;
        final int kOutMax = (int) Math.sqrt(N);
        final double kuAve = 10;
        final int kuMin = 1;
        final int kuMax = 1000;
        final int m0 = 6;
        final int m = 6;
        final double gamma = 2.5;
        final int swapNum = 0;
        final boolean loadFromEdgeList = false;

        final int batchSize = 10;
        final int itrs = 10;
        final double pMin = 0.0;
        final double pMax = 1.0;
        final double pStep = 0.01;
        final double[] pList = ArrayUtils.arange(pMin, pMax, pStep);
    }
}
