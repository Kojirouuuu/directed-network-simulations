package sim.network.topology;

import sim.network.DirectedGraph;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * 空白区切りエッジリスト（1行が "a b" の形式）から有向グラフ a→b を読み込む。
 * 例: twitter_combined.txt のような "src dst" の並びから DirectedGraph を構築する。
 *
 * <p>データソース: SNAP (Stanford Network Analysis Project)
 * <a href="https://snap.stanford.edu/data/">https://snap.stanford.edu/data/</a>
 * には Twitter (ego-Twitter フォロー関係)、Google+、Wikipedia などの有向グラフが多数収録されている。
 * 本クラスでは Twitter ego-network（twitter_combined.txt 等）のフォロー関係を使用する。
 */
public final class EgoTwitter {

    private static final String DEFAULT_PATH1 = "app/data/real-network/twitter_combined.txt";
    private static final String DEFAULT_PATH2 = "data/real-network/twitter_combined.txt";

    private EgoTwitter() {}

    /**
     * デフォルトのエッジリストを読み込む。
     * まず {@code app/data/real-network/twitter_combined.txt} を探し、
     * 見つからなければ {@code data/real-network/twitter_combined.txt} を試す（いずれもカレントディレクトリ基準）。
     *
     * @return 構築された DirectedGraph
     * @throws IOException 両パスでファイルが見つからない場合、または読み込みエラー
     */
    public static DirectedGraph loadFromDefaultEdgeList() throws IOException {
        Path p = Path.of(DEFAULT_PATH1);
        if (!Files.exists(p)) {
            p = Path.of(DEFAULT_PATH2);
        }
        if (!Files.exists(p)) {
            throw new IOException("Edge list not found. Tried: " + DEFAULT_PATH1 + " and " + DEFAULT_PATH2);
        }
        return loadFromEdgeList("ego-Twitter", p);
    }

    /**
     * 指定パスのファイルを読み、各行を "a b" と解釈して有向辺 a→b の DirectedGraph を返す。
     * 空行および '#' で始まる行は無視する。頂点IDは出現順に 0..n-1 にマッピングする。
     *
     * @param name グラフ名（null の場合はファイル名を使用）
     * @param path エッジリストファイルのパス
     * @return 構築された DirectedGraph
     * @throws IOException ファイル読み込みエラー
     */
    public static DirectedGraph loadFromEdgeList(String name, Path path) throws IOException {
        return DirectedGraph.loadFromEdgeList(name, path);
    }

    /**
     * 指定パス文字列からエッジリストを読み、有向グラフを構築する。
     *
     * @param name グラフ名（null の場合はファイル名を使用）
     * @param pathString ファイルパス文字列（例: "/Users/.../twitter_combined.txt"）
     * @return 構築された DirectedGraph
     * @throws IOException ファイル読み込みエラー
     */
    public static DirectedGraph loadFromEdgeList(String name, String pathString) throws IOException {
        return loadFromEdgeList(name, Path.of(pathString));
    }
}
