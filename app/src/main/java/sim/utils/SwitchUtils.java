package sim.utils;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedBA;
import sim.network.topology.DirectedCM;
import sim.network.topology.DirectedCMInPow;
import sim.network.topology.DirectedCMOutPow;
import sim.network.topology.undirected.ER;
import sim.network.topology.undirected.BA;

import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * ネットワークタイプに応じたパス構成・グラフ生成などの switch ロジックを一元化するユーティリティ。
 * 新規ネットワークタイプを追加する場合は、{@link #buildNetworkPath} および {@link #generateGraph} に
 * case
 * を追加する。
 */
public final class SwitchUtils {

    private static final String SAR_PREFIX = "sar";

    private SwitchUtils() {
    }

    /**
     * シミュレーションレベルの出力ディレクトリを構築する。
     *
     * @param optionPath オプション識別子（例: config.optionPath）
     * @param threshold 閾値
     * @return 出力ディレクトリ（例: sar/{optionPath}/threshold={threshold}）
     */
    public static Path buildSimulationOutputDir(String optionPath, int threshold) {
        return Paths.get("out", SAR_PREFIX, optionPath, String.format("threshold=%d", threshold));
    }

    /**
     * ネットワークパス（シミュレーションディレクトリ以下の相対パス）を構築する。
     * 使用されないパラメータは null でよい。ネットワークタイプに応じて必要なものだけ参照される。
     *
     * @param networkType ネットワークタイプ（DirectedCM, DirectedCMInPow, DirectedCMOutPow,
     *        ER, BA, DirectedBA）
     * @param N 頂点数
     * @param kuAve DirectedCM 用（null 可）
     * @param kInMin DirectedCMInPow 用（null 可）
     * @param kOutMin DirectedCMOutPow 用（null 可）
     * @param z ER 用の平均次数（null 可）
     * @param m0 DirectedBA 用（null 可）
     * @param m DirectedBA 用（null 可）
     * @return ネットワークパス（例: {networkType}/N={N}/{networkSpecific}/）
     */
    public static Path buildNetworkPath(String networkType, int N,
            Integer kuAve, Integer kInMin, Integer kOutMin, Double z, Integer m0, Integer m) {
        String networkSpecific = switch (networkType) {
            case "DirectedCM" -> String.format("kuAve=%d", requireNonNull(kuAve, "DirectedCM requires kuAve"));
            case "DirectedCMInPow" -> String.format("kInMin=%d",
                    requireNonNull(kInMin, "DirectedCMInPow requires kInMin"));
            case "DirectedCMOutPow" -> String.format("kOutMin=%d",
                    requireNonNull(kOutMin, "DirectedCMOutPow requires kOutMin"));
            case "ER" -> String.format("z=%.2f", requireNonNull(z, "ER requires z"));
            case "DirectedBA", "BA" -> String.format("m0=%d/m=%d", requireNonNull(m0, "BA requires m0"),
                    requireNonNull(m, "BA requires m"));
            default -> throw new IllegalArgumentException("Unknown network type: " + networkType);
        };

        return Paths.get(String.format("%s/N=%d/%s", networkType, N, networkSpecific));
    }

    private static <T> T requireNonNull(T value, String message) {
        if (value == null) {
            throw new IllegalArgumentException(message);
        }
        return value;
    }

    /**
     * ネットワークタイプと設定値から DirectedGraph を生成する。
     * 使用されないパラメータは null または 0 でよい。
     *
     * @param networkType ネットワークタイプ
     * @param N 頂点数
     * @param kHat DirectedCM 用（null 可）
     * @param kInMin DirectedCMInPow 用（null 可）
     * @param kInMax DirectedCMInPow 用（null 可）
     * @param kOutMin DirectedCMOutPow 用（null 可）
     * @param kOutMax DirectedCMOutPow 用（null 可）
     * @param kuAve 平均次数（DirectedCMInPow, DirectedCMOutPow, ER）
     * @param gamma DirectedCMInPow, DirectedCMOutPow 用
     * @param m0 DirectedBA 用（null 可）
     * @param m DirectedBA 用（null 可）
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generateGraph(String networkType, int N,
            Integer kHat, Integer kInMin, Integer kInMax, Integer kOutMin, Integer kOutMax,
            double kuAve, double gamma, Integer m0, Integer m, long seed) {
        GraphGeneratorParams params = new GraphGeneratorParams(
                networkType, N, kHat, kInMin, kInMax, kOutMin, kOutMax, kuAve, gamma, m0, m);
        return generateGraph(params, seed);
    }

    /**
     * ネットワークタイプとパラメータから DirectedGraph を生成する。
     *
     * @param params グラフ生成に必要なパラメータ（ネットワークタイプに応じて必要なフィールドのみ設定）
     * @param seed 乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generateGraph(GraphGeneratorParams params, long seed) {
        String name = params.networkType();
        int N = params.N();
        return switch (params.networkType()) {
            case "DirectedCM" -> DirectedCM.generate(name, N,
                    requireNonNull(params.kHat(), "DirectedCM requires kHat"), seed);
            case "DirectedCMInPow" -> DirectedCMInPow.generate(name, N,
                    requireNonNull(params.kInMin(), "DirectedCMInPow requires kInMin"),
                    requireNonNull(params.kInMax(), "DirectedCMInPow requires kInMax"),
                    params.kuAve(), params.gamma(), seed);
            case "DirectedCMOutPow" -> DirectedCMOutPow.generate(name, N,
                    requireNonNull(params.kOutMin(), "DirectedCMOutPow requires kOutMin"),
                    requireNonNull(params.kOutMax(), "DirectedCMOutPow requires kOutMax"),
                    params.kuAve(), params.gamma(), seed);
            case "ER" -> ER.generateERFromKAve(N, params.kuAve(), seed);
            case "DirectedBA" -> DirectedBA.generate(name, N,
                    requireNonNull(params.m0(), "DirectedBA requires m0"),
                    requireNonNull(params.m(), "DirectedBA requires m"), seed);
            case "BA" -> BA.generate(name, N,
                    requireNonNull(params.m0(), "BA requires m0"),
                    requireNonNull(params.m(), "BA requires m"), seed);
            default -> throw new IllegalArgumentException("Unknown network type: " + params.networkType());
        };
    }

    /**
     * グラフ生成に必要なパラメータを保持するレコード。
     * ネットワークタイプに応じて必要なフィールドのみ設定すればよい。
     */
    public record GraphGeneratorParams(
            String networkType,
            int N,
            Integer kHat,
            Integer kInMin,
            Integer kInMax,
            Integer kOutMin,
            Integer kOutMax,
            double kuAve,
            double gamma,
            Integer m0,
            Integer m) {
        public static GraphGeneratorParams forDirectedCM(int N, int kHat) {
            return new GraphGeneratorParams("DirectedCM", N, kHat, null, null, null, null, 0, 0, null, null);
        }

        public static GraphGeneratorParams forDirectedCMInPow(int N, int kInMin, int kInMax, double kuAve,
                double gamma) {
            return new GraphGeneratorParams("DirectedCMInPow", N, null, kInMin, kInMax, null, null, kuAve, gamma, null,
                    null);
        }

        public static GraphGeneratorParams forDirectedCMOutPow(int N, int kOutMin, int kOutMax, double kuAve,
                double gamma) {
            return new GraphGeneratorParams("DirectedCMOutPow", N, null, null, null, kOutMin, kOutMax, kuAve, gamma,
                    null, null);
        }

        public static GraphGeneratorParams forER(int N, double kuAve) {
            return new GraphGeneratorParams("ER", N, null, null, null, null, null, kuAve, 0, null, null);
        }

        public static GraphGeneratorParams forBA(int N, int m0, int m) {
            return new GraphGeneratorParams("BA", N, null, null, null, null, null, 0, 0, m0, m);
        }
    }

}
