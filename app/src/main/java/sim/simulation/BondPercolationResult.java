package sim.simulation;

import sim.utils.PathsEx;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

/**
 * ボンドパーコレーションのシミュレーション結果。
 * GWCC（最大弱連結成分）と GSCC（最大強連結成分）の統計を保持する。
 */
public final class BondPercolationResult {
    public final int n;

    // GWCC（弱連結成分）
    public final int gwccSize;
    public final int gwccSecondSize;
    public final double gwccMeanSize;

    // GSCC（強連結成分）
    public final int gsccSize;
    public final int gsccSecondSize;
    public final double gsccMeanSize;

    BondPercolationResult(int n,
            int gwccSize, int gwccSecondSize, double gwccMeanSize,
            int gsccSize, int gsccSecondSize, double gsccMeanSize) {
        this.n = n;
        this.gwccSize = gwccSize;
        this.gwccSecondSize = gwccSecondSize;
        this.gwccMeanSize = gwccMeanSize;
        this.gsccSize = gsccSize;
        this.gsccSecondSize = gsccSecondSize;
        this.gsccMeanSize = gsccMeanSize;
    }

    /**
     * 最終状態 CSV を出力する。
     *
     * @param path   出力先パス
     * @param itr    イテレーション番号
     * @param p      パーコレーション確率
     * @param append true のとき追記、false のとき新規ファイル（番号付き）
     * @throws IOException ファイル書き込みエラー
     */
    public void writeFinalStateCsv(Path path, int itr, double p, boolean append) throws IOException {
        if (!append) {
            path = PathsEx.resolveIndexed(path);
        }
        Files.createDirectories(path.getParent());
        boolean writeHeader = true;
        if (Files.exists(path)) {
            try {
                writeHeader = Files.size(path) == 0L;
            } catch (IOException ignored) {
                // fallback to writing header
            }
        }
        try (BufferedWriter bw = Files.newBufferedWriter(path,
                StandardOpenOption.CREATE,
                append ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING);
                PrintWriter out = new PrintWriter(bw)) {
            if (writeHeader) {
                out.println("itr,p,gwcc_size,gwcc_fraction,gwcc_second_size,gwcc_mean_size,gscc_size,gscc_fraction,gscc_second_size,gscc_mean_size");
            }
            out.printf(Locale.ROOT, "%d,%.9f,%d,%.9f,%d,%.9f,%d,%.9f,%d,%.9f%n",
                    itr, p,
                    gwccSize, (double) gwccSize / n,
                    gwccSecondSize, gwccMeanSize,
                    gsccSize, (double) gsccSize / n,
                    gsccSecondSize, gsccMeanSize);
        }
    }
}
