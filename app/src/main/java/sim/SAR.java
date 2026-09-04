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
import java.time.Instant;
import java.util.List;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
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
        final long plannedSimulationCount = config.runSarSimulations ? totalTasks : 0L;
        final int parallelism = Runtime.getRuntime().availableProcessors();
        SARRunSummary summary = createRunSummary(config, plannedSimulationCount, parallelism);
        Thread shutdownHook = new Thread(summary::finishInterruptedBestEffort, "sar-summary-shutdown-hook");
        Runtime.getRuntime().addShutdownHook(shutdownHook);
        Throwable failure = null;

        try {
            DirectedGraph sharedJointDegreeSource = config.randomizationMode == RandomizationMode.JOINT_DEGREE_CM
                    ? loadGraph(config, 0)
                    : null;

            System.out.println("Total tasks: " + totalTasks);
            if (sharedJointDegreeSource == null) {
                System.out.println(config.networkType + ": N=" + config.N + ", itrs=" + config.itrs);
            } else {
                int targetN = checkedExpandedSize(sharedJointDegreeSource.n, config.sizeMultiplier, "vertex");
                System.out.println(config.networkType + ": sourceN=" + sharedJointDegreeSource.n
                        + ", sizeMultiplier=" + config.sizeMultiplier + ", targetN=" + targetN
                        + ", itrs=" + config.itrs);
            }
            System.out.println("Graph randomization: "
                    + effectiveRandomizationMode(config.randomizationMode, config.sizeMultiplier).pathLabel());
            System.out.println("Edge-list output: " + (config.writeEdgeList ? "enabled" : "disabled"));
            System.out.println("SAR simulations: " + (config.runSarSimulations ? "enabled" : "disabled"));

            System.out.println("Parallelism: " + parallelism + " (available processors)");

            int[] progressItr = new int[config.batchSize];
            AtomicLong done = new AtomicLong(0);
            AtomicBoolean running = new AtomicBoolean(true);

            Thread renderer = createTotalProgressRenderer(done, totalTasks, running);
            renderer.start();

            try (ForkJoinPool pool = new ForkJoinPool(parallelism)) {
                Future<?> future = pool.submit(() -> IntStream.range(0, config.batchSize).parallel()
                        .forEach(batchIndex -> processBatch(
                                batchIndex, config, sharedJointDegreeSource, progressItr, done, summary)));

                future.get();
            } finally {
                running.set(false);
                renderer.join();
            }

            System.out.println("All tasks completed");
        } catch (Exception e) {
            failure = e;
            summary.recordUnhandledErrorIfAbsent(e);
            throw e;
        } catch (Error e) {
            failure = e;
            summary.recordUnhandledErrorIfAbsent(e);
            throw e;
        } finally {
            try {
                if (failure == null) {
                    summary.finishSucceeded();
                } else {
                    summary.finishFailed();
                }
                if (summary.summaryPath() != null) {
                    System.out.println("Run summary: " + summary.summaryPath());
                }
            } catch (IOException e) {
                if (failure == null) {
                    throw e;
                }
                failure.addSuppressed(e);
                System.err.println("Failed to finalize SAR run summary: " + e.getMessage());
            } finally {
                try {
                    Runtime.getRuntime().removeShutdownHook(shutdownHook);
                } catch (IllegalStateException ignored) {
                    // JVM shutdown is already in progress; the hook will handle finalization.
                }
            }
        }
    }

    /**
     * 1つのバッチを処理する。
     *
     * @param batchIndex バッチインデックス
     * @param config シミュレーション設定
     * @param progressItr 進捗記録用配列
     * @param done 完了タスク数のカウンタ
     * @param summary 実行サマリー
     */
    private static void processBatch(int batchIndex, SimulationConfig config,
            DirectedGraph sharedJointDegreeSource,
            int[] progressItr, AtomicLong done, SARRunSummary summary) {
        Instant batchStartedAt = Instant.now();
        summary.recordBatchStarted(batchIndex, batchStartedAt);
        String stage = "graph-loading";
        SARRunSummary.ErrorContext errorContext = SARRunSummary.ErrorContext.forBatch(batchIndex);

        try {
            DirectedGraph g;
            if (sharedJointDegreeSource != null) {
                stage = "graph-randomization";
                g = applyRandomization(
                        sharedJointDegreeSource, config.randomizationMode,
                        GRAPH_RANDOMIZATION_BASE_SEED + batchIndex,
                        config.loadFromEdgeList, config.sizeMultiplier);
            } else {
                g = loadGraph(config, batchIndex);
                stage = "graph-randomization";
                g = applyRandomization(
                        g, config.randomizationMode,
                        GRAPH_RANDOMIZATION_BASE_SEED + batchIndex,
                        config.loadFromEdgeList, config.sizeMultiplier);
            }

            Path outputDirectory = buildOutputDirectory(g, config);
            stage = "summary-initialization";
            summary.registerBatchGraph(batchIndex, g.n, g.m, outputDirectory);

            if (config.writeEdgeList && !config.loadFromEdgeList) {
                stage = "edge-list-output";
                writeEdgeList(g, batchIndex, config);
            }

            if (!config.runSarSimulations) {
                progressItr[batchIndex] = config.itrs;
                done.incrementAndGet();
                summary.recordBatchCompleted(batchIndex, Instant.now());
                return;
            }

            stage = "result-path-preparation";
            Path resultsPath = prepareOutputPath(outputDirectory, batchIndex);
            summary.recordResultPath(batchIndex, resultsPath);

            for (int itr = 0; itr < config.itrs; itr++) {
                progressItr[batchIndex] = itr;

                for (int ri = 0; ri < config.rho0List.length; ri++) {
                    double rho0 = config.rho0List[ri];
                    for (int li = 0; li < config.lambdaDirectedList.length; li++) {
                        double lambdaDirected = config.lambdaDirectedList[li];
                        for (int lni = 0; lni < config.lambdaNondirectedList.length; lni++) {
                            double lambdaNondirected = config.lambdaNondirectedList[lni];
                            errorContext = SARRunSummary.ErrorContext.forSimulation(
                                    batchIndex, itr, rho0, lambdaDirected, lambdaNondirected, config.mu);
                            stage = "simulation";

                            int[] thresholdList = new int[g.n];
                            Arrays.fill(thresholdList, config.threshold);

                            runSimulation(g, config, lambdaDirected, lambdaNondirected, config.mu, rho0,
                                    thresholdList, batchIndex, itr, resultsPath);

                            done.incrementAndGet();
                            summary.recordSimulationCompleted();
                        }
                    }
                }
            }

            progressItr[batchIndex] = config.itrs;
            summary.recordBatchCompleted(batchIndex, Instant.now());
        } catch (RuntimeException | Error e) {
            summary.recordError(stage, errorContext, e);
            if (errorContext.iteration() != null) {
                summary.recordSimulationFailed();
            }
            summary.recordBatchFailed(batchIndex, Instant.now());
            throw e;
        }
    }

    /** 設定に従って、ランダム化前のグラフを読み込むか生成する。 */
    private static DirectedGraph loadGraph(SimulationConfig config, int batchIndex) {
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
                return DirectedGraph.loadFromEdgeList(config.networkType, edgeListPath);
            } catch (IOException e) {
                throw new RuntimeException("Failed to load edge list for batch " + batchIndex, e);
            }
        }
        return SwitchUtils.generateGraph(config.networkType, config.N,
                null, config.kdMin, config.kdMax, config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kuMin, config.kuMax,
                config.kuAve, config.gamma, config.m0, config.m, config.swapNum,
                config.gammaIn, config.gammaOut, config.corrA,
                GRAPH_BASE_SEED + batchIndex);
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
        return applyRandomization(graph, mode, seed, loadedFromEdgeList, 1);
    }

    /** 選択したモードとサイズ倍率をグラフへ適用する。 */
    static DirectedGraph applyRandomization(DirectedGraph graph, RandomizationMode mode,
            long seed, boolean loadedFromEdgeList, int sizeMultiplier) {
        if (graph == null) {
            throw new IllegalArgumentException("graph must be non-null");
        }
        if (mode == null) {
            throw new IllegalArgumentException("mode must be non-null");
        }
        if (sizeMultiplier < 1) {
            throw new IllegalArgumentException("sizeMultiplier must be at least 1");
        }
        if (mode != RandomizationMode.JOINT_DEGREE_CM && sizeMultiplier != 1) {
            throw new IllegalArgumentException(
                    "sizeMultiplier can exceed 1 only in JOINT_DEGREE_CM mode");
        }

        return switch (mode) {
            case NONE -> graph;
            case EDGE_SWAP -> loadedFromEdgeList ? graph : graph.randomizeByEdgeSwaps(seed);
            case SHUFFLE_IN_DEGREES -> graph.randomizeByShuffledDegreeSequence(DegreeSide.IN, seed);
            case SHUFFLE_OUT_DEGREES -> graph.randomizeByShuffledDegreeSequence(DegreeSide.OUT, seed);
            case JOINT_DEGREE_CM -> sizeMultiplier == 1
                    ? graph
                    : graph.expandByRepeatedJointDegreeSequence(sizeMultiplier, seed);
        };
    }

    /** 倍率1では元ネットワークを使うため、出力上はランダム化なしとして扱う。 */
    static RandomizationMode effectiveRandomizationMode(RandomizationMode mode, int sizeMultiplier) {
        if (mode == null) {
            throw new IllegalArgumentException("mode must be non-null");
        }
        if (sizeMultiplier < 1) {
            throw new IllegalArgumentException("sizeMultiplier must be at least 1");
        }
        return mode == RandomizationMode.JOINT_DEGREE_CM && sizeMultiplier == 1
                ? RandomizationMode.NONE
                : mode;
    }

    private static int checkedExpandedSize(int originalSize, int multiplier, String quantity) {
        if (multiplier < 1) {
            throw new IllegalArgumentException("sizeMultiplier must be at least 1");
        }
        try {
            return Math.multiplyExact(originalSize, multiplier);
        } catch (ArithmeticException e) {
            throw new IllegalArgumentException(
                    "expanded " + quantity + " count exceeds the supported int range", e);
        }
    }

    /** 最終的なランダマイズ方式を区別する出力パスを返す。 */
    static Path appendRandomizationPath(Path networkPath, RandomizationMode mode) {
        return appendRandomizationPath(networkPath, mode, 1);
    }

    /** サイズ倍率を考慮して、最終的なランダマイズ方式を区別する出力パスを返す。 */
    static Path appendRandomizationPath(Path networkPath, RandomizationMode mode, int sizeMultiplier) {
        if (networkPath == null) {
            throw new IllegalArgumentException("networkPath must be non-null");
        }
        RandomizationMode effectiveMode = effectiveRandomizationMode(mode, sizeMultiplier);
        return networkPath.resolve("randomization=" + effectiveMode.pathLabel());
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
                .resolve(appendRandomizationPath(
                        networkPath, config.randomizationMode, config.sizeMultiplier))
                .resolve(String.format("%d.csv", batchIndex));
        try {
            g.writeEdgeList(edgeListPath);
        } catch (IOException e) {
            throw new RuntimeException("Failed to write edge list: " + edgeListPath, e);
        }
    }

    /**
     * 実際のグラフサイズを使って出力ディレクトリを構築する。
     *
     * @param g グラフ
     * @param config シミュレーション設定
     * @return 出力ディレクトリ
     */
    private static Path buildOutputDirectory(DirectedGraph g, SimulationConfig config) {
        Path outputDir = SwitchUtils.buildSimulationOutputDir(config.optionPath, config.threshold);
        Path networkPath = SwitchUtils.buildNetworkPath(
                config.networkType, g.n,
                null, config.kuAve,
                config.kInMin, config.kInMax, config.kOutMin, config.kOutMax,
                config.kdMin, config.kdMax, config.kuMin, config.kuMax, config.m0, config.m,
                config.gamma, config.swapNum,
                config.gammaIn, config.gammaOut, config.corrA);
        return outputDir.resolve(
                appendRandomizationPath(
                        networkPath, config.randomizationMode, config.sizeMultiplier));
    }

    private static Path prepareOutputPath(Path outputDirectory, int batchIndex) {
        String idx = String.format("%02d", batchIndex);
        return PathsEx.resolveIndexed(
                outputDirectory.resolve(String.format("results_%s.csv", idx)));
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
            System.err.println("CSV output error (batch " + batchIndex + ", iteration " + itr + ", rho0 "
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

    private static SARRunSummary createRunSummary(
            SimulationConfig config, long plannedSimulationCount, int parallelism) {
        Map<String, Object> network = new LinkedHashMap<>();
        network.put("name", config.networkType);
        network.put("configuredN", config.N);
        network.put("sizeMultiplier", config.sizeMultiplier);
        network.put("randomizationMode",
                effectiveRandomizationMode(config.randomizationMode, config.sizeMultiplier).pathLabel());
        network.put("loadFromEdgeList", config.loadFromEdgeList);
        network.put("writeEdgeList", config.writeEdgeList);

        Map<String, Object> topology = new LinkedHashMap<>();
        topology.put("kdMin", config.kdMin);
        topology.put("kdMax", config.kdMax);
        topology.put("kInMin", config.kInMin);
        topology.put("kInMax", config.kInMax);
        topology.put("kOutMin", config.kOutMin);
        topology.put("kOutMax", config.kOutMax);
        topology.put("kuAve", config.kuAve);
        topology.put("kuMin", config.kuMin);
        topology.put("kuMax", config.kuMax);
        topology.put("m0", config.m0);
        topology.put("m", config.m);
        topology.put("gamma", config.gamma);
        topology.put("swapNum", config.swapNum);
        topology.put("gammaIn", config.gammaIn);
        topology.put("gammaOut", config.gammaOut);
        topology.put("corrA", config.corrA);
        network.put("topologyParameters", topology);

        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("rho0Values", config.rho0List);
        parameters.put("lambdaDirectedValues", config.lambdaDirectedList);
        parameters.put("lambdaNondirectedValues", config.lambdaNondirectedList);
        parameters.put("lambdaDirectedMin", config.lambdaDirectedMin);
        parameters.put("lambdaDirectedMax", config.lambdaDirectedMax);
        parameters.put("lambdaDirectedStep", config.lambdaDirectedStep);
        parameters.put("mu", config.mu);
        parameters.put("threshold", config.threshold);
        parameters.put("tMax", config.tMax);
        parameters.put("dt", config.dt);

        Map<String, Object> seeds = new LinkedHashMap<>();
        seeds.put("rngBaseSeed", RNG_BASE_SEED);
        seeds.put("simulationBaseSeed", SIM_BASE_SEED);
        seeds.put("graphBaseSeed", GRAPH_BASE_SEED);
        seeds.put("graphRandomizationBaseSeed", GRAPH_RANDOMIZATION_BASE_SEED);
        seeds.put("nodeSeedOffset", SEED_OFFSET_NODES);
        parameters.put("seeds", seeds);

        Map<String, Object> execution = new LinkedHashMap<>();
        execution.put("optionPath", config.optionPath);
        execution.put("batchSize", config.batchSize);
        execution.put("iterationsPerBatch", config.itrs);
        execution.put("parallelism", parallelism);
        execution.put("runSarSimulations", config.runSarSimulations);
        execution.put("finalStateOnly", config.isFinal);
        execution.put("useGillespie", config.useGillespie);

        return new SARRunSummary(plannedSimulationCount, network, parameters, execution);
    }

    /** SAR で使用するグラフランダマイズ方式。 */
    enum RandomizationMode {
        NONE("none", false), EDGE_SWAP("edge-swap", true), SHUFFLE_IN_DEGREES("in-degree-shuffle",
                true), SHUFFLE_OUT_DEGREES("out-degree-shuffle", true), JOINT_DEGREE_CM("joint-degree-cm", false);

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
        final int sizeMultiplier = 1; // 実ネットワークの同時次数分布を使う場合の頂点数倍率

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
