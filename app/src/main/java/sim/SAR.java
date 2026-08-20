package sim;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.DegreeSide;
import sim.simulation.SARGillespieSimulator;
import sim.simulation.SARSimulator;
import sim.simulation.SARResult;
import sim.utils.ArrayUtils;
import sim.utils.PathsEx;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
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
    private static final long GRAPH_RANDOMIZATION_BASE_SEED = 4242L; // グラフランダマイズ用のベースシード

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
        final int lambdaNondirectedCount = config.lambdaNondirectedList.length;
        final long totalTasks = config.runSarSimulations
                ? (long) config.batchSize * config.itrs * rho0Count * lambdaDirectedCount
                        * lambdaNondirectedCount
                : config.batchSize;

        System.out.println("Total tasks: " + totalTasks);
        System.out.println(config.networkType + ": N=" + config.N + ", itrs=" + config.itrs);
        System.out.println("Graph randomization: " + config.randomizationMode.pathLabel());
        System.out.println("Edge-list output: " + (config.writeEdgeList ? "enabled" : "disabled"));
        System.out.println("SAR simulations: " + (config.runSarSimulations ? "enabled" : "disabled"));

        int parallelism = Runtime.getRuntime().availableProcessors();
        System.out.println("Parallelism: " + parallelism + " (available processors)");

        int[] progressItr = new int[config.batchSize];
        AtomicLong done = new AtomicLong(0);
        AtomicBoolean running = new AtomicBoolean(true);

        Thread renderer = createTotalProgressRenderer(done, totalTasks, running);
        renderer.start();

        try (ForkJoinPool pool = new ForkJoinPool(parallelism)) {
            Future<?> future = pool.submit(() -> IntStream.range(0, config.batchSize).parallel()
                    .forEach(batchIndex -> processBatch(batchIndex, config, progressItr, done)));

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
     * @param batchIndex バッチインデックス
     * @param config シミュレーション設定
     * @param progressItr 進捗記録用配列
     * @param done 完了タスク数のカウンタ
     */
    private static void processBatch(int batchIndex, SimulationConfig config,
            int[] progressItr, AtomicLong done) {
        DirectedGraph g;
        if (config.loadFromEdgeList) {
            Path networkPath = SwitchUtils.buildNetworkPath(
                    config.networkType, config.N,
                    null, config.kuAve,
                    config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                    config.kdMin, config.kdMax, config.kuMin, config.kuMax, config.m0, config.m,
                    config.gamma, config.swapNum,
                    config.gammaIn, config.gammaOut, config.corrA);
            try {
                Path edgeListPath = resolveEdgeListPath(
                        Paths.get("out/edgelist"), networkPath,
                        config.randomizationMode.usesEdgeSwappedInput(), batchIndex);
                g = DirectedGraph.loadFromEdgeList(config.networkType, edgeListPath);
            } catch (IOException e) {
                throw new RuntimeException("Failed to load edge list for batch " + batchIndex, e);
            }
        } else {
            g = SwitchUtils.generateGraph(config.networkType, config.N,
                    null, config.kdMin, config.kdMax, config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                    config.kuMin, config.kuMax,
                    config.kuAve, config.gamma, config.m0, config.m, config.swapNum,
                    config.gammaIn, config.gammaOut, config.corrA,
                    GRAPH_BASE_SEED + batchIndex);
        }

        g = applyRandomization(
                g, config.randomizationMode,
                GRAPH_RANDOMIZATION_BASE_SEED + batchIndex,
                config.loadFromEdgeList);

        if (config.writeEdgeList && !config.loadFromEdgeList) {
            writeEdgeList(g, batchIndex, config);
        }

        if (!config.runSarSimulations) {
            progressItr[batchIndex] = config.itrs;
            done.incrementAndGet();
            return;
        }

        Path resultsPath = prepareOutputPath(g, batchIndex, config);

        for (int itr = 0; itr < config.itrs; itr++) {
            progressItr[batchIndex] = itr;

            for (int ri = 0; ri < config.rho0List.length; ri++) {
                double rho0 = config.rho0List[ri];
                for (int li = 0; li < config.lambdaDirectedList.length; li++) {
                    double lambdaDirected = config.lambdaDirectedList[li];
                    for (int lni = 0; lni < config.lambdaNondirectedList.length; lni++) {
                        double lambdaNondirected = config.lambdaNondirectedList[lni];

                        int[] thresholdList = new int[g.n];
                        Arrays.fill(thresholdList, config.threshold);

                        runSimulation(g, config, lambdaDirected, lambdaNondirected, config.mu, rho0,
                                thresholdList, batchIndex, itr, resultsPath);

                        done.incrementAndGet();
                    }
                }
            }
        }

        progressItr[batchIndex] = config.itrs;
    }

    /**
     * 読み込むエッジリストを解決する。通常の保存先を優先し、GraphGen が実頂点数や
     * 追加条件のサブディレクトリを使って保存したファイルも検索対象にする。
     */
    static Path resolveEdgeListPath(Path edgeListRoot, Path networkPath,
            boolean randomizedByEdgeSwaps, int batchIndex) throws IOException {
        Path randomizationPath = SwitchUtils.appendRandomizationPath(
                networkPath, randomizedByEdgeSwaps);
        String fileName = String.format("%d.csv", batchIndex);
        Path expectedPath = edgeListRoot.resolve(randomizationPath).resolve(fileName);
        if (Files.isRegularFile(expectedPath)) {
            return expectedPath;
        }

        // まず同じ N・ネットワーク条件の追加サブディレクトリ内を探す。
        List<Path> candidates = findEdgeLists(
                edgeListRoot.resolve(randomizationPath), fileName, null);
        if (candidates.isEmpty()) {
            // 実データ由来のグラフでは設定上の N と読み込み後の実頂点数が異なり得る。
            String randomizationDir = "randomization="
                    + (randomizedByEdgeSwaps ? "edge-swap" : "none");
            candidates = findEdgeLists(
                    edgeListRoot.resolve(networkPath.getName(0)), fileName, randomizationDir);
        }

        if (candidates.size() == 1) {
            return candidates.getFirst();
        }
        if (candidates.isEmpty()) {
            throw new IOException("Edge list not found: " + expectedPath);
        }
        throw new IOException("Multiple matching edge lists found for " + expectedPath + ": " + candidates);
    }

    private static List<Path> findEdgeLists(Path searchRoot, String fileName,
            String requiredDirectory) throws IOException {
        if (!Files.isDirectory(searchRoot)) {
            return List.of();
        }
        try (var paths = Files.walk(searchRoot)) {
            return paths
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().equals(fileName))
                    .filter(path -> requiredDirectory == null
                            || containsDirectory(path, requiredDirectory))
                    .sorted()
                    .toList();
        }
    }

    private static boolean containsDirectory(Path path, String directoryName) {
        for (Path part : path) {
            if (part.toString().equals(directoryName)) {
                return true;
            }
        }
        return false;
    }

    /** 選択したモードをグラフへ適用する。 */
    static DirectedGraph applyRandomization(DirectedGraph graph, RandomizationMode mode,
            long seed, boolean loadedFromEdgeList) {
        if (graph == null) {
            throw new IllegalArgumentException("graph must be non-null");
        }
        if (mode == null) {
            throw new IllegalArgumentException("mode must be non-null");
        }

        return switch (mode) {
            case NONE -> graph;
            case EDGE_SWAP -> loadedFromEdgeList ? graph : graph.randomizeByEdgeSwaps(seed);
            case SHUFFLE_IN_DEGREES -> graph.randomizeByShuffledDegreeSequence(DegreeSide.IN, seed);
            case SHUFFLE_OUT_DEGREES -> graph.randomizeByShuffledDegreeSequence(DegreeSide.OUT, seed);
        };
    }

    /** 最終的なランダマイズ方式を区別する出力パスを返す。 */
    static Path appendRandomizationPath(Path networkPath, RandomizationMode mode) {
        if (networkPath == null) {
            throw new IllegalArgumentException("networkPath must be non-null");
        }
        if (mode == null) {
            throw new IllegalArgumentException("mode must be non-null");
        }
        return networkPath.resolve("randomization=" + mode.pathLabel());
    }

    /**
     * 生成したグラフを、読み込み処理と共通のパス規約で書き出す。
     *
     * @param g グラフ
     * @param batchIndex バッチインデックス
     * @param config シミュレーション設定
     */
    private static void writeEdgeList(DirectedGraph g, int batchIndex, SimulationConfig config) {
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, g.n,
                null, config.kuAve,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kdMin, config.kdMax, config.kuMin, config.kuMax, config.m0, config.m,
                config.gamma, config.swapNum,
                config.gammaIn, config.gammaOut, config.corrA);
        Path edgeListPath = Paths.get("out/edgelist")
                .resolve(appendRandomizationPath(networkPath, config.randomizationMode))
                .resolve(String.format("%d.csv", batchIndex));
        try {
            g.writeEdgeList(edgeListPath);
        } catch (IOException e) {
            throw new RuntimeException("Failed to write edge list: " + edgeListPath, e);
        }
    }

    /**
     * 出力パスを準備する。
     *
     * @param g グラフ
     * @param batchIndex バッチインデックス
     * @param config シミュレーション設定
     * @return 出力パス
     */
    private static Path prepareOutputPath(DirectedGraph g, int batchIndex, SimulationConfig config) {
        String idx = String.format("%02d", batchIndex);
        Path outputDir = SwitchUtils.buildSimulationOutputDir(config.optionPath, config.threshold);
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, g.n,
                null, config.kuAve,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kdMin, config.kdMax, config.kuMin, config.kuMax, config.m0, config.m,
                config.gamma, config.swapNum,
                config.gammaIn, config.gammaOut, config.corrA);
        Path basePath = outputDir.resolve(
                appendRandomizationPath(networkPath, config.randomizationMode));
        return PathsEx.resolveIndexed(
                basePath.resolve(String.format("results_%s.csv", idx)));
    }

    /**
     * 1回のシミュレーションを実行する。
     *
     * @param g グラフ
     * @param config シミュレーション設定
     * @param lambdaDirected 有向辺の感染率
     * @param lambdaNondirected 無向辺の感染率
     * @param mu 回復率
     * @param rho0 初期感染率
     * @param thresholdList 各ノードの閾値リスト
     * @param batchIndex バッチインデックス
     * @param itr イテレーション番号
     * @param resultsPath 結果出力パス
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

        SARResult res;
        if (config.useGillespie) {
            res = config.isFinal
                    ? SARGillespieSimulator.simulate(
                            g, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                            thresholdList, init, simSeed)
                    : SARGillespieSimulator.simulate(
                            g, lambdaDirected, lambdaNondirected, config.mu, config.tMax, config.dt,
                            thresholdList, init, simSeed);
        } else {
            res = config.isFinal
                    ? SARSimulator.simulate(
                            g, lambdaDirected, lambdaNondirected, config.mu, config.tMax,
                            thresholdList, init, simSeed)
                    : SARSimulator.simulate(
                            g, lambdaDirected, lambdaNondirected, config.mu, config.tMax, config.dt,
                            thresholdList, init, simSeed);
        }

        try {
            if (config.isFinal) {
                res.writeFinalStateCsv(resultsPath, batchIndex * config.itrs + itr, rho0, lambdaDirected,
                        lambdaNondirected, config.mu,
                        true);
            } else {
                res.writeTimeSeriesCsv(resultsPath, batchIndex * config.itrs + itr, rho0, lambdaDirected,
                        lambdaNondirected, config.mu,
                        true);
            }
        } catch (IOException e) {
            System.out.println("CSV output error (batch " + batchIndex + ", iteration " + itr + ", lambdaDirected "
                    + rho0 + ", lambdaDirected " + lambdaDirected + ", lambdaNondirected " + lambdaNondirected + ", mu "
                    + config.mu
                    + "): " + e.getMessage());
            throw new RuntimeException(e);
        }
    }

    /**
     * 進捗表示スレッドを作成する。
     *
     * @param done 完了タスク数のカウンタ
     * @param totalTasks 総タスク数
     * @param running 実行中フラグ
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
     * @param done 完了タスク数
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

    /** SAR で使用するグラフランダマイズ方式。 */
    enum RandomizationMode {
        NONE("none", false),
        EDGE_SWAP("edge-swap", true),
        SHUFFLE_IN_DEGREES("in-degree-shuffle", true),
        SHUFFLE_OUT_DEGREES("out-degree-shuffle", true);

        private final String pathLabel;
        private final boolean usesEdgeSwappedInput;

        RandomizationMode(String pathLabel, boolean usesEdgeSwappedInput) {
            this.pathLabel = pathLabel;
            this.usesEdgeSwappedInput = usesEdgeSwappedInput;
        }

        String pathLabel() {
            return pathLabel;
        }

        boolean usesEdgeSwappedInput() {
            return usesEdgeSwappedInput;
        }
    }

    /**
     * シミュレーション設定を保持する内部クラス。
     */
    private static class SimulationConfig {
        // ネットワークの基本設定
        final String networkType = "rev-ego-Twitter"; // ネットワークタイプ
        final String optionPath = "lambda-ugokasu-real-2"; // オプションパス
        final int N = 500_000; // 頂点数

        // 次数パラメータ
        final int kdMin = 5; // 最小次数
        final int kdMax = (int) Math.sqrt(N); // 最大次数
        final int kInMin = 5; // 最小入次数
        final int kInMax = (int) Math.sqrt(N); // 最大入次数
        // final int kInMax = N; // 最大入次数
        final int kOutMin = 5; // 最小出次数
        final int kOutMax = (int) Math.sqrt(N); // 最大出次数
        // final int kOutMax = N; // 最大出次数
        final double kuAve = 10; // 平均次数
        final int kuMin = 5; // 最小次数
        final int kuMax = (int) Math.sqrt(N); // 最大次数

        // トポロジー生成パラメータ
        final int m0 = 6; // 初期完全グラフの頂点数
        final int m = 6; // 各新規ノードが接続する辺（弧）の数
        final double gamma = 2.5;
        final int swapNum = 0; // PowPow 用（null のとき 0 として扱う）

        // SchwartzDirectedSF 用パラメータ（他のネットワークタイプでは未使用）
        final Double gammaIn = gamma; // λ_in
        final Double gammaOut = gamma; // λ_out
        final Double corrA = null; // 相関確率 A ∈ [0, 1]

        // 入出力・グラフ処理
        /**
         * true のとき
         * out/edgelist/{networkPath}/randomization={mode}/{batchIndex}.csv
         * からグラフを読み込む（次数列シャッフル時は edge-swap 配下を使用）
         */
        final boolean loadFromEdgeList = true;
        final RandomizationMode randomizationMode = RandomizationMode.SHUFFLE_IN_DEGREES;
        final boolean writeEdgeList = false; // 生成したネットワークのエッジリストを書き出すか
        final boolean runSarSimulations = true; // SAR シミュレーションを実行するか

        // 実行回数
        final int batchSize = 10; // バッチサイズ
        final int itrs = 10; // イテレーション数

        // SAR シミュレーション設定
        final boolean isFinal = true; // 最終状態のみ出力するか
        final double dt = 0.1; // isFinal == false の時は dt 刻みで記録する。
        final boolean useGillespie = false; // true: Gillespie方式, false: イベント駆動方式
        final double mu = 1.0; // 回復率
        final double tMax = 200.0; // シミュレーション終了時刻

        // 伝播率
        final double lambdaDirectedMin = 0.0;
        final double lambdaDirectedMax = 2.0;
        final double lambdaDirectedStep = 0.02;
        final double[] lambdaDirectedList = ArrayUtils.arange(lambdaDirectedMin,
                lambdaDirectedMax, lambdaDirectedStep); // 有向辺の感染率
        // final double[] lambdaDirectedList = { 0.1, 0.2, 1.0, 2.0, 5.0, 10.0 };

        // final double lambdaNondirectedMin = 0.0;
        // final double lambdaNondirectedMax = 2.0;
        // final double lambdaNondirectedStep = 0.02;
        // final double[] lambdaNondirectedList =
        // ArrayUtils.arange(lambdaNondirectedMin, lambdaNondirectedMax,
        // lambdaNondirectedStep); // 無向辺の感染率
        final double[] lambdaNondirectedList = { 0.0 };

        // 初期採用率・閾値
        // final double rho0Min = 1.0e-5;
        // final double rho0Max = 2.0e-1;
        // final double rho0Step = 0.0005;
        // final int rho0Count = 100;
        // final double[] rho0List = ArrayUtils.arange(rho0Min, rho0Max, rho0Step);
        // final double[] rho0List = ArrayUtils.geomspace(rho0Min, rho0Max, rho0Count);
        final double[] rho0List = { 1e-2, 1e-3, 1e-4 }; // 初期感染率のリスト
        final int threshold = 3; // 閾値
    }
}
