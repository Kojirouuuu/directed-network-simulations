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
 *
 * @see EgoTwitter#loadFromEdgeList(String, Path)
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
        return EgoTwitter.loadFromEdgeList("higgs-social", p);
    }
}
