package sim;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.DegreeSide;
import sim.utils.SwitchUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import sim.SAR.SimulationConfig;
import sim.SAR.RandomizationMode;

/** Graph input, transformations and output paths for the SAR experiment. */
final class SARGraphPreparation {
    private SARGraphPreparation() {}

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


    /** 設定に従って、ランダム化前のグラフを読み込むか生成する。 */
    static DirectedGraph loadGraph(SimulationConfig config, int batchIndex) {
        if (config.loadFromEdgeList) {
            Path networkPath = SwitchUtils.buildNetworkPath(networkParameters(config));
            try {
                Path edgeListPath = resolveEdgeListPath(
                        Paths.get("out/edgelist"), networkPath,
                        config.randomizationMode.usesEdgeSwappedInput(), batchIndex);
                return DirectedGraph.loadFromEdgeList(config.networkType, edgeListPath);
            } catch (IOException e) {
                throw new RuntimeException("Failed to load edge list for batch " + batchIndex, e);
            }
        }
        return SwitchUtils.generateGraph(networkParameters(config), SAR.GRAPH_BASE_SEED + batchIndex);
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

    static List<Path> findEdgeLists(Path searchRoot, String fileName,
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

    static boolean containsDirectory(Path path, String directoryName) {
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

    static int checkedExpandedSize(int originalSize, int multiplier, String quantity) {
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
    static void writeEdgeList(DirectedGraph g, int batchIndex, SimulationConfig config) {
        Path networkPath = SwitchUtils.buildNetworkPath(networkParameters(config).withN(g.n));
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
    static Path buildOutputDirectory(DirectedGraph g, SimulationConfig config) {
        Path outputDir = SwitchUtils.buildSimulationOutputDir(config.optionPath, config.threshold);
        Path networkPath = SwitchUtils.buildNetworkPath(networkParameters(config).withN(g.n));
        return outputDir.resolve(
                appendRandomizationPath(
                        networkPath, config.randomizationMode, config.sizeMultiplier));
    }

}
