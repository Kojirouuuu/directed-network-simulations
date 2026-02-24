package sim.network.topology;

import sim.network.DirectedGraph;

import java.io.IOException;
import java.nio.file.Path;

/**
 * EgoTwitter の生成メソッドを import して使用し、
 * 読み込んだ有向グラフの辺を全て逆向きにしたグラフを提供する。
 * 同一のエッジリスト（twitter_combined.txt 等）を読み、各辺 a→b を b→a にした DirectedGraph を返す。
 *
 * @see EgoTwitter
 * @see DirectedGraph#reverseDirectedEdges()
 */
public final class RevEgoTwitter {

    private RevEgoTwitter() {
    }

    /**
     * デフォルトのエッジリストを EgoTwitter で読み込み、辺を反転したグラフを返す。
     *
     * @return 辺を逆向きにした DirectedGraph（名前は "RevEgoTwitter"）
     * @throws IOException エッジリストが見つからない場合、または読み込みエラー
     */
    public static DirectedGraph loadReversedFromDefaultEdgeList() throws IOException {
        DirectedGraph g = EgoTwitter.loadFromDefaultEdgeList();
        return g.reverseDirectedEdges("RevEgoTwitter");
    }

    /**
     * 指定パスのエッジリストを EgoTwitter で読み込み、辺を反転したグラフを返す。
     *
     * @param name グラフ名（null の場合は "RevEgoTwitter"）
     * @param path エッジリストファイルのパス
     * @return 辺を逆向きにした DirectedGraph
     * @throws IOException ファイル読み込みエラー
     */
    public static DirectedGraph loadReversedFromEdgeList(String name, Path path) throws IOException {
        DirectedGraph g = EgoTwitter.loadFromEdgeList(null, path);
        String newName = name != null ? name : "RevEgoTwitter";
        return g.reverseDirectedEdges(newName);
    }

    /**
     * 指定パス文字列のエッジリストを EgoTwitter で読み込み、辺を反転したグラフを返す。
     *
     * @param name グラフ名（null の場合は "RevEgoTwitter"）
     * @param pathString ファイルパス文字列
     * @return 辺を逆向きにした DirectedGraph
     * @throws IOException ファイル読み込みエラー
     */
    public static DirectedGraph loadReversedFromEdgeList(String name, String pathString) throws IOException {
        return loadReversedFromEdgeList(name, Path.of(pathString));
    }
}
