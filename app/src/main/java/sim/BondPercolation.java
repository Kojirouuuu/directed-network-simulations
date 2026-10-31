package sim;

import sim.utils.ProgressBar;
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
                .corrA(config.corrA).build();
    }

    private static void processBatch(int batchIndex, SimulationConfig config, AtomicLong done) {
        DirectedGraph g;
        Path resultsPath;

        if (config.loadFromEdgeList) {
            Path networkPath = SwitchUtils.buildNetworkPath(networkParameters(config));
            Path edgeListPath = Paths.get("out/edgelist").resolve(networkPath)
                    .resolve(String.format("%d.csv", batchIndex));
            try {
                g = DirectedGraph.loadFromEdgeList(config.networkType, edgeListPath);
            } catch (IOException e) {
                throw new RuntimeException("Failed to load edge list: " + edgeListPath, e);
            }
        } else {
            g = SwitchUtils.generateGraph(networkParameters(config), GRAPH_BASE_SEED + batchIndex);
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
        Path networkPath = SwitchUtils.buildNetworkPath(networkParameters(config).withN(g.n));
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
        ProgressBar.render(done, total);
    }

    private static class SimulationConfig {
        final String networkType = "DirectedCMInPow";
        final String optionPath = "2607";
        final int N = 100_000;
        final int kdMin = 2;
        final int kdMax = 1000;
        final int kInMin = 2;
        final int kInMax = (int) Math.sqrt(N);
        final int kOutMin = 2;
        final int kOutMax = (int) Math.sqrt(N);
        final double kuAve = 0;
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

        final boolean loadFromEdgeList = false;

        final int batchSize = 16;
        final int itrs = 8;
        final double pMin = 0.0;
        final double pMax = 1.0;
        final double pStep = 0.005;
        final double[] pList = ArrayUtils.arange(pMin, pMax, pStep);
    }
}
