package sim.network.topology;

import sim.network.DirectedGraph;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Higgs Twitter ソーシャルネットワークのエッジリストから有向グラフを読み込む。
 * 空白区切り "src dst" 形式。EgoTwitter と同じ読み込みロジックを使用する。
 *
 * <p>
 * データソース: SNAP (Stanford Network Analysis Project)
 * <a href="https://snap.stanford.edu/data/">https://snap.stanford.edu/data/</a>
 * の higgs-social_network.edgelist を使用する。
 */
public final class HiggsSocial {

    private static final String DEFAULT_PATH1 = "app/data/real-network/higgs-social_network.edgelist";
    private static final String DEFAULT_PATH2 = "data/real-network/higgs-social_network.edgelist";

    private HiggsSocial() {
    }

    /**
     * デフォルトのエッジリストを読み込む。
     * まず {@code app/data/real-network/higgs-social_network.edgelist} を探し、
     * 見つからなければ {@code data/real-network/higgs-social_network.edgelist} を試す。
     *
     * @return 構築された DirectedGraph（名前は "higgs-social"）
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
        return loadFromEdgeList("higgs-social", p);
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
