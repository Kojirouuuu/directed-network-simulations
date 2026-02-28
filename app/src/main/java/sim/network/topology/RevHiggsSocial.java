package sim.network.topology;

import sim.network.DirectedGraph;

import java.io.IOException;

/**
 * HiggsSocial の有向グラフを読み込み、辺を全て逆向きにしたグラフを提供する。
 *
 * @see HiggsSocial
 * @see DirectedGraph#reverseDirectedEdges(String)
 */
public final class RevHiggsSocial {

    private RevHiggsSocial() {
    }

    /**
     * デフォルトのエッジリストを HiggsSocial で読み込み、辺を反転したグラフを返す。
     *
     * @return 辺を逆向きにした DirectedGraph（名前は "rev-higgs-social"）
     * @throws IOException エッジリストが見つからない場合、または読み込みエラー
     */
    public static DirectedGraph loadReversedFromDefaultEdgeList() throws IOException {
        DirectedGraph g = HiggsSocial.loadFromDefaultEdgeList();
        return g.reverseDirectedEdges("rev-higgs-social");
    }
}
