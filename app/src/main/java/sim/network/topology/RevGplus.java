package sim.network.topology;

import sim.network.DirectedGraph;

import java.io.IOException;

/**
 * Gplus の有向グラフを読み込み、辺を全て逆向きにしたグラフを提供する。
 *
 * @see Gplus
 * @see DirectedGraph#reverseDirectedEdges(String)
 */
public final class RevGplus {

    private RevGplus() {
    }

    /**
     * デフォルトのエッジリストを Gplus で読み込み、辺を反転したグラフを返す。
     *
     * @return 辺を逆向きにした DirectedGraph（名前は "rev-gplus"）
     * @throws IOException エッジリストが見つからない場合、または読み込みエラー
     */
    public static DirectedGraph loadReversedFromDefaultEdgeList() throws IOException {
        DirectedGraph g = Gplus.loadFromDefaultEdgeList();
        return g.reverseDirectedEdges("rev-gplus");
    }
}
