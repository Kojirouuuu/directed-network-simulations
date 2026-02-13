package sim.network.topology;

import sim.network.DirectedGraph;

/**
 * 入次数をパワーロー分布で、出次数をランダムに配分する。
 * DirectedCMOutPow で生成したグラフの有向辺を全て逆向きにしたものとして定義する。
 * （OutPow: 出次数冪・入次数ランダム → 反転 → 入次数冪・出次数ランダム = InPow）
 */
public final class DirectedCMInPow {

    private DirectedCMInPow() {}

    /**
     * 指定されたパラメータから Directed + Nondirected CM を生成する。
     * 入次数はパワーロー分布、出次数は総数をランダムに配分する。
     * 実装は DirectedCMOutPow で出次数冪分布のグラフを生成し、有向辺を反転している。
     *
     * @param name    グラフ名
     * @param n       頂点数
     * @param kInMin  最小入次数
     * @param kInMax  最大入次数
     * @param kuAve   平均無向次数
     * @param gamma   パワーロー指数
     * @param seed    乱数シード
     * @return 生成された DirectedGraph
     */
    public static DirectedGraph generate(String name, int n, int kInMin, int kInMax, double kuAve, double gamma, long seed) {
        DirectedGraph gOut = DirectedCMOutPow.generate(name, n, kInMin, kInMax, kuAve, gamma, seed);
        return gOut.reverseDirectedEdges(name);
    }
}
