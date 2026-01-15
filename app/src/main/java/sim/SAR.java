package sim;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedCM;
import sim.simulation.SARSimulator;
import sim.simulation.SARResult;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;

/**
 * SAR（Susceptible-Adopted-Recovered）シミュレーションのメインクラス。
 * 並列処理によるバッチシミュレーションを実行する。
 */
public class SAR {
    private static final int PROGRESS_BAR_LENGTH = 100; // 進捗バーの長さ
    private static final int PROGRESS_UPDATE_INTERVAL_MS = 100; // 進捗更新間隔（ミリ秒）

    private static final long RNG_BASE_SEED = 7L; // 乱数生成用のベースシード（閾値リスト・ノードシャッフル用）
    private static final long SIM_BASE_SEED = 12345L; // シミュレーション用のベースシード
    private static final long GRAPH_BASE_SEED = 42L; // グラフ生成用のベースシード

    // シードオフセット（異なる目的で異なるシードを生成するため）
    private static final long SEED_OFFSET_NODES = 2000L; // ノードシャッフル用オフセット
    
    /**
     * メインメソッド。
     *
     * @param args コマンドライン引数
     * @throws Exception 実行エラー
     */
    public static void main(String[] args) throws Exception {
        SimulationConfig config = new SimulationConfig();

        final int rho0Count = config.rho0List.length;
        final int lambdaDirectedCount = config.lambdaDirectedList.length;
        final long totalTasks = (long) config.batchSize * config.itrs * rho0Count * lambdaDirectedCount;
        System.out.println("Total tasks: " + totalTasks);
        System.out.println(config.networkType + ": N=" + config.N + ", itrs=" + config.itrs);

        int parallelism = Runtime.getRuntime().availableProcessors();
        System.out.println("Parallelism: " + parallelism + " (available processors)");

        int[] progressItr = new int[config.batchSize];
        AtomicLong done = new AtomicLong(0);
        AtomicBoolean running = new AtomicBoolean(true);

        Thread renderer = createTotalProgressRenderer(done, totalTasks, running);
        renderer.start();

        try (ForkJoinPool pool = new ForkJoinPool(parallelism)) {
            Future<?> future = pool.submit(() ->
                    IntStream.range(0, config.batchSize).parallel().forEach(batchIndex ->
                            processBatch(batchIndex, config, progressItr, done, totalTasks)
                    )
            );

            future.get();
        } finally {
            running.set(false);
            renderer.join();
        }

        System.out.println("All tasks completed");
    }
    
    /**
     * 1つのバッチを処理する。
     *
     * @param batchIndex  バッチインデックス
     * @param config      シミュレーション設定
     * @param progressItr 進捗記録用配列
     * @param done        完了タスク数のカウンタ
     * @param totalTasks  総タスク数
     */
    private static void processBatch(int batchIndex, SimulationConfig config,
            int[] progressItr, AtomicLong done, long totalTasks) {
        DirectedGraph g = switch (config.networkType) {
            case "DirectedCM" -> DirectedCM.generate("DirectedCM", config.N, config.kHat,
                    GRAPH_BASE_SEED + batchIndex);
            default -> throw new IllegalArgumentException("Unknown network type: " + config.networkType);
        };

        Path resultsPath = prepareOutputPath(g, batchIndex, config);

        for (int itr = 0; itr < config.itrs; itr++) {
            progressItr[batchIndex] = itr;

            for (int ri = 0; ri < config.rho0List.length; ri++) {
                double rho0 = config.rho0List[ri];
                for (int li = 0; li < config.lambdaDirectedList.length; li++) {
                    double lambdaDirected = config.lambdaDirectedList[li];

                    int[] thresholdList = new int[config.N];
                    Arrays.fill(thresholdList, config.threshold);

                    runSimulation(g, config, lambdaDirected, config.lambdaNondirected, config.mu, rho0,
                            thresholdList, batchIndex, itr, resultsPath);

                    done.incrementAndGet();
                }
            }
        }

        progressItr[batchIndex] = config.itrs;
    }
    
    /**
     * 出力パスを準備する。
     *
     * @param g          グラフ
     * @param batchIndex バッチインデックス
     * @param config     シミュレーション設定
     * @return 出力パス
     */
    private static Path prepareOutputPath(DirectedGraph g, int batchIndex, SimulationConfig config) {
        String idx = String.format("%02d", batchIndex);
        String networkPath = g.name;

        Path basePath = Paths.get(String.format("out/fastsar/%s/threshold=%d/N=%d", networkPath, config.threshold,
                config.N));
        return PathsEx.resolveIndexed(
                basePath.resolve(String.format("results_%s.csv", idx))
        );
    }
    
    /**
     * 1回のシミュレーションを実行する。
     *
     * @param g                   グラフ
     * @param config              シミュレーション設定
     * @param lambdaDirected      有向辺の感染率
     * @param lambdaNondirected   無向辺の感染率
     * @param mu                  回復率
     * @param rho0                初期感染率
     * @param thresholdList       各ノードの閾値リスト
     * @param batchIndex          バッチインデックス
     * @param itr                 イテレーション番号
     * @param resultsPath         結果出力パス
     */
    private static void runSimulation(DirectedGraph g, SimulationConfig config,
            double lambdaDirected, double lambdaNondirected, double mu, double rho0, int[] thresholdList,
            int batchIndex, int itr, Path resultsPath) {

        int initialInfectedNum = (int) (g.n * rho0);
        if (initialInfectedNum <= 0) {
            initialInfectedNum = 1;
        }

        int[] nodes = new int[g.n];
        for (int i = 0; i < g.n; i++) {
            nodes[i] = i;
        }
        long nodeSeed = RNG_BASE_SEED + (long) batchIndex * config.itrs + itr + SEED_OFFSET_NODES;
        ArrayUtils.shuffle(nodes, nodeSeed);
        int[] init = Arrays.copyOfRange(nodes, 0, initialInfectedNum);

        long simSeed = SIM_BASE_SEED + (long) batchIndex * config.itrs + itr;

        SARResult res = SARSimulator.simulate(
                g, lambdaDirected, config.lambdaNondirected, config.mu, config.tMax, thresholdList, init,
                simSeed
        );

        try {
            if (config.isFinal) {
                res.writeFinalStateCsv(resultsPath, itr, lambdaDirected, config.lambdaNondirected, config.mu,
                        true);
            } else {
                res.writeTimeSeriesCsv(resultsPath, itr, lambdaDirected, config.lambdaNondirected, config.mu,
                        true);
            }
        } catch (IOException e) {
            System.out.println("CSV output error (batch " + batchIndex + ", iteration " + itr + ", lambdaDirected "
                    + lambdaDirected + ", lambdaNondirected " + config.lambdaNondirected + ", mu " + config.mu
                    + "): " + e.getMessage());
            throw new RuntimeException(e);
        }
    }
    
    /**
     * 進捗表示スレッドを作成する。
     *
     * @param done       完了タスク数のカウンタ
     * @param totalTasks 総タスク数
     * @param running    実行中フラグ
     * @return 進捗表示スレッド
     */
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

                if (d >= totalTasks) {
                    break;
                }

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

    /**
     * 進捗バーを表示する。
     *
     * @param done  完了タスク数
     * @param total 総タスク数
     */
    private static void renderTotalProgressBar(long done, long total) {
        int filled = (total == 0) ? PROGRESS_BAR_LENGTH
                : (int) Math.min(PROGRESS_BAR_LENGTH, (done * PROGRESS_BAR_LENGTH) / total);

        int percent = (total == 0) ? 100 : (int) Math.min(100, (done * 100) / total);

        String bar = "#".repeat(filled) + "-".repeat(PROGRESS_BAR_LENGTH - filled);

        // 行をクリアしてから進捗を表示（\033[2K は行全体をクリア、\r は行頭に戻る）
        System.out.print("\033[2K\rProgress [%s] %3d%% (%d/%d)".formatted(bar, percent, done, total));
        System.out.flush();
    }
    
    /**
     * シミュレーション設定を保持する内部クラス。
     */
    private static class SimulationConfig {
        final String networkType = "DirectedCM"; // ネットワークタイプ
        final int N = 1_000_000; // 頂点数
        final int kHat = 4; // 平均次数
        final boolean isFinal = true; // 最終状態のみ出力するか
        final int batchSize = 16; // バッチサイズ
        final int itrs = 1000; // イテレーション数
        final double mu = 1.0; // 回復率
        final double tMax = 200.0; // シミュレーション終了時刻
        final double lambdaDirectedMin = 0.0;
        final double lambdaDirectedMax = 0.2;
        final double lambdaDirectedStep = 0.002;
        final double[] lambdaDirectedList = ArrayUtils.arange(lambdaDirectedMin, lambdaDirectedMax, lambdaDirectedStep); // 有向辺の感染率
        final double lambdaNondirected = 0.2; // 無向辺の感染率
        final double[] rho0List = { (double) 1 / N }; // 初期感染率のリスト
        final int threshold = 1; // 閾値
    }
}
