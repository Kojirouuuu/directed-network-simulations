package sim.utils;

import sim.network.DirectedGraph;
import sim.network.topology.DirectedBA;
import sim.network.topology.DirectedCM;
import sim.network.topology.DirectedCMInPow;
import sim.network.topology.DirectedCMOutPow;
import sim.network.topology.PowPow;
import sim.network.topology.SameInOut;
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
         *        PowPow, SameInOut, ER, BA, DirectedBA）
         * @param N 頂点数
         * @param kdAve DirectedCM 用（null 可）
         * @param kuAve ER 用（null 可）
         * @param kInMin DirectedCMInPow 用（null 可）
         * @param kInMax DirectedCMInPow 用（null 可）
         * @param kOutMin DirectedCMOutPow 用（null 可）
         * @param kOutMax DirectedCMOutPow 用（null 可）
         * @param kdMin PowPow 用（null 可）
         * @param kdMax PowPow 用（null 可）
         * @param m0 DirectedBA 用（null 可）
         * @param m DirectedBA 用（null 可）
         * @param gamma DirectedCMInPow, DirectedCMOutPow, PowPow 用（null 可）
         * @param swapNum PowPow 用（null 可）
         * @return ネットワークパス（例: {networkType}/N={N}/{networkSpecific}/）
         */
        public static Path buildNetworkPath(String networkType, int N,
                        Integer kdAve, Double kuAve, Integer kInMin, Integer kInMax, Integer kOutMin, Integer kOutMax,
                        Integer kdMin, Integer kdMax, Integer m0, Integer m,
                        Double gamma, Integer swapNum) {
                String networkSpecific = switch (networkType) {
                        case "DirectedCM" -> String.format("kdAve=%d",
                                        requireNonNull(kdAve, "DirectedCM requires kdAve"));
                        case "DirectedCMInPow" -> String.format("gamma=%.2f/kInMin=%d/kInMax=%d",
                                        requireNonNull(gamma, "DirectedCMInPow requires gamma"),
                                        requireNonNull(kInMin, "DirectedCMInPow requires kInMin"),
                                        requireNonNull(kInMax, "DirectedCMInPow requires kInMax"));
                        case "DirectedCMOutPow" -> String.format("gamma=%.2f/kOutMin=%d/kOutMax=%d",
                                        requireNonNull(gamma, "DirectedCMOutPow requires gamma"),
                                        requireNonNull(kOutMin, "DirectedCMOutPow requires kOutMin"),
                                        requireNonNull(kOutMax, "DirectedCMOutPow requires kOutMax"));
                        case "PowPow" -> String.format("gamma=%.2f/kdMin=%d/kdMax=%d/swapNum=%d",
                                        requireNonNull(gamma, "PowPow requires gamma"),
                                        requireNonNull(kdMin, "PowPow requires kdMin"),
                                        requireNonNull(kdMax, "PowPow requires kdMax"),
                                        requireNonNull(swapNum, "PowPow requires swapNum"));
                        case "SameInOut" -> String.format("gamma=%.2f/kdMin=%d/kdMax=%d",
                                        requireNonNull(gamma, "SameInOut requires gamma"),
                                        requireNonNull(kdMin, "SameInOut requires kdMin"),
                                        requireNonNull(kdMax, "SameInOut requires kdMax"));
                        case "ER" -> String.format("kuAve=%.2f", requireNonNull(kuAve, "ER requires kuAve"));
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
         * @param kdAve DirectedCM 用（null 可）
         * @param kInMin DirectedCMInPow 用（null 可）
         * @param kInMax DirectedCMInPow 用（null 可）
         * @param kOutMin DirectedCMOutPow 用（null 可）
         * @param kOutMax DirectedCMOutPow 用（null 可）
         * @param kuAve 平均次数（DirectedCMInPow, DirectedCMOutPow, ER）
         * @param gamma DirectedCMInPow, DirectedCMOutPow, PowPow, SameInOut 用（null 可）
         * @param m0 DirectedBA 用（null 可）
         * @param m DirectedBA 用（null 可）
         * @param swapNum PowPow 用（null 可、0 として扱う）
         * @param seed 乱数シード
         * @return 生成された DirectedGraph
         */
        public static DirectedGraph generateGraph(String networkType, int N,
                        Integer kdAve, Integer kdMin, Integer kdMax, Integer kInMin, Integer kInMax, Integer kOutMin,
                        Integer kOutMax,
                        Double kuAve, Double gamma, Integer m0, Integer m, Integer swapNum, long seed) {
                GraphGeneratorParams params = new GraphGeneratorParams(
                                networkType, N, kdAve, kdMin, kdMax, kInMin, kInMax, kOutMin, kOutMax, kuAve, gamma, m0,
                                m, swapNum);
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
                                        requireNonNull(params.kdAve(), "DirectedCM requires kdAve"), seed);
                        case "DirectedCMInPow" -> DirectedCMInPow.generate(name, N,
                                        requireNonNull(params.kInMin(), "DirectedCMInPow requires kInMin"),
                                        requireNonNull(params.kInMax(), "DirectedCMInPow requires kInMax"),
                                        requireNonNull(params.kuAve(), "DirectedCMInPow requires kuAve"),
                                        requireNonNull(params.gamma(), "DirectedCMInPow requires gamma"), seed);
                        case "DirectedCMOutPow" -> DirectedCMOutPow.generate(name, N,
                                        requireNonNull(params.kOutMin(), "DirectedCMOutPow requires kOutMin"),
                                        requireNonNull(params.kOutMax(), "DirectedCMOutPow requires kOutMax"),
                                        requireNonNull(params.kuAve(), "DirectedCMOutPow requires kuAve"),
                                        requireNonNull(params.gamma(), "DirectedCMOutPow requires gamma"), seed);
                        case "ER" -> ER.generateERFromKAve(N, requireNonNull(params.kuAve(), "ER requires kuAve"),
                                        seed);
                        case "DirectedBA" -> DirectedBA.generate(name, N,
                                        requireNonNull(params.m0(), "DirectedBA requires m0"),
                                        requireNonNull(params.m(), "DirectedBA requires m"), seed);
                        case "BA" -> BA.generate(name, N,
                                        requireNonNull(params.m0(), "BA requires m0"),
                                        requireNonNull(params.m(), "BA requires m"), seed);
                        case "PowPow" -> PowPow.generate(name, N,
                                        requireNonNull(params.kdMin(), "PowPow requires kdMin"),
                                        requireNonNull(params.kdMax(), "PowPow requires kdMax"),
                                        requireNonNull(params.gamma(), "PowPow requires gamma"),
                                        requireNonNull(params.swapNum(), "PowPow requires swapNum"), seed);
                        case "SameInOut" -> SameInOut.generate(name, N,
                                        requireNonNull(params.kdMin(), "SameInOut requires kdMin"),
                                        requireNonNull(params.kdMax(), "SameInOut requires kdMax"),
                                        params.gamma(), seed);
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
                        Integer kdAve,
                        Integer kdMin,
                        Integer kdMax,
                        Integer kInMin,
                        Integer kInMax,
                        Integer kOutMin,
                        Integer kOutMax,
                        Double kuAve,
                        Double gamma,
                        Integer m0,
                        Integer m,
                        Integer swapNum) {
                public static GraphGeneratorParams forDirectedCM(int N, int kdAve) {
                        return new GraphGeneratorParams("DirectedCM", N, kdAve, null, null, null, null, null, null,
                                        null, null, null, null, null);
                }

                public static GraphGeneratorParams forDirectedCMInPow(int N, int kInMin, int kInMax, Double kuAve,
                                Double gamma) {
                        return new GraphGeneratorParams("DirectedCMInPow", N, null, null, null, kInMin, kInMax, null,
                                        null, kuAve, gamma, null, null, null);
                }

                public static GraphGeneratorParams forDirectedCMOutPow(int N, int kOutMin, int kOutMax, Double kuAve,
                                Double gamma) {
                        return new GraphGeneratorParams("DirectedCMOutPow", N, null, null, null, null, null, kOutMin,
                                        kOutMax, kuAve, gamma, null, null, null);
                }

                public static GraphGeneratorParams forPowPow(int N, int kMin, int kMax, Double gamma, int swapNum) {
                        return new GraphGeneratorParams("PowPow", N, null, kMin, kMax, null, null, null, null,
                                        null, gamma, null, null, swapNum);
                }

                public static GraphGeneratorParams forSameInOut(int N, int kMin, int kMax, Double gamma) {
                        return new GraphGeneratorParams("SameInOut", N, null, kMin, kMax, null, null, null, null,
                                        null, gamma, null, null, null);
                }

                public static GraphGeneratorParams forER(int N, Double kuAve) {
                        return new GraphGeneratorParams("ER", N, null, null, null, null, null, null, null, kuAve, null,
                                        null, null, null);
                }

                public static GraphGeneratorParams forBA(int N, int m0, int m) {
                        return new GraphGeneratorParams("BA", N, null, null, null, null, null, null, null, null, null,
                                        m0, m, null);
                }
        }

}
