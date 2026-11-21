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
 * GWCC（最大弱連結成分）と GSCC（最大強連結成分）の統計、および全成分のサイズ分布を保持する。
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

    private final int[] wccSizeCounts;
    private final int[] sccSizeCounts;

    BondPercolationResult(int n,
            int gwccSize, int gwccSecondSize, double gwccMeanSize,
            int gsccSize, int gsccSecondSize, double gsccMeanSize,
            int[] wccSizeCounts, int[] sccSizeCounts) {
        this.n = n;
        this.gwccSize = gwccSize;
        this.gwccSecondSize = gwccSecondSize;
        this.gwccMeanSize = gwccMeanSize;
        this.gsccSize = gsccSize;
        this.gsccSecondSize = gsccSecondSize;
        this.gsccMeanSize = gsccMeanSize;
        this.wccSizeCounts = wccSizeCounts;
        this.sccSizeCounts = sccSizeCounts;
    }

    /**
     * Writes size frequencies for all WCCs and SCCs, including largest components
     * and isolated vertices. Only nonzero counts are emitted, in ascending size order.
     * When append is false, an existing file is preserved using an indexed name.
     */
    public void writeComponentSizeDistributionCsv(Path path, int itr, double p, boolean append) throws IOException {
        if (!append) {
            path = PathsEx.resolveIndexed(path);
        }
        if (path.getParent() != null) {
            Files.createDirectories(path.getParent());
        }
        boolean writeHeader = !Files.exists(path) || Files.size(path) == 0;
        try (BufferedWriter out = Files.newBufferedWriter(path,
                StandardOpenOption.CREATE,
                append ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING)) {
            if (writeHeader) {
                out.write("itr,p,component_type,size,count");
                out.newLine();
            }
            writeSizeCounts(out, itr, p, "wcc", wccSizeCounts);
            writeSizeCounts(out, itr, p, "scc", sccSizeCounts);
        }
    }

    private static void writeSizeCounts(BufferedWriter out, int itr, double p, String type,
            int[] counts) throws IOException {
        for (int size = 1; size < counts.length; size++) {
            if (counts[size] != 0) {
                out.write(String.format(Locale.ROOT, "%d,%.9f,%s,%d,%d", itr, p, type, size, counts[size]));
                out.newLine();
            }
        }
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
