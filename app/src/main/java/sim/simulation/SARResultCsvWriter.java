package sim.simulation;

import sim.utils.PathsEx;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

/** CSV serialization for SAR results. */
final class SARResultCsvWriter {
    private final SARResult result;

    SARResultCsvWriter(SARResult result) { this.result = result; }

    /**
     * 集計時系列CSVを出力する。
     * 解析時にパラメータも横に展開したいケース向け。
     *
     * @param path   出力先のパス
     * @param itr    イテレーション番号
     * @param lambdaDirected  有向辺の感染率
     * @param lambdaNondirected 無向辺の感染率
     * @param mu 回復率
     * @param append 追記モードの場合 true
     * @throws IOException ファイル書き込みエラー
     */
    public void writeTimeSeriesCsv(Path path, int itr, double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        writeTimeSeriesCsvInternal(path, itr, null, null, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    /**
     * 入次数・出次数分布ラベルを含む集計時系列CSVを出力する。
     */
    public void writeTimeSeriesCsv(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        writeTimeSeriesCsvInternal(path, itr, inDegreeDistribution, outDegreeDistribution,
                rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    private void writeTimeSeriesCsvInternal(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
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
                if (inDegreeDistribution == null) {
                    out.println("itr,rho_0,lambda_d,lambda_u,mu,time,A,R,Phi");
                } else {
                    out.println("itr,in_degree_distribution,out_degree_distribution,rho_0,lambda_d,lambda_u,mu,time,A,R,Phi");
                }
            }
            for (int i = 0; i < result.times.size(); i++) {
                if (inDegreeDistribution == null) {
                    out.printf(Locale.ROOT, "%d,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                            itr, rho0, lambdaDirected, lambdaNondirected, mu, result.times.get(i), result.A.get(i), result.R.get(i),
                            result.Phi.get(i));
                } else {
                    out.printf(Locale.ROOT, "%d,%s,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                            itr, inDegreeDistribution, outDegreeDistribution, rho0, lambdaDirected,
                            lambdaNondirected, mu, result.times.get(i), result.A.get(i), result.R.get(i), result.Phi.get(i));
                }
            }
        }
    }

    /**
     * 最終状態CSVを出力する。
     * 解析時にパラメータも横に展開したいケース向け。
     *
     * @param path   出力先のパス
     * @param itr    イテレーション番号
     * @param lambdaDirected  有向辺の感染率
     * @param lambdaNondirected 無向辺の感染率
     * @param mu 回復率
     * @param append 追記モードの場合 true
     * @throws IOException ファイル書き込みエラー
     */
    public void writeFinalStateCsv(Path path, int itr, double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        writeFinalStateCsvInternal(path, itr, null, null, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    /**
     * 入次数・出次数分布ラベルを含む最終状態CSVを出力する。
     */
    public void writeFinalStateCsv(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        writeFinalStateCsvInternal(path, itr, inDegreeDistribution, outDegreeDistribution,
                rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    private void writeFinalStateCsvInternal(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
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
                if (inDegreeDistribution == null) {
                    out.println("itr,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi");
                } else {
                    out.println("itr,in_degree_distribution,out_degree_distribution,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi");
                }
            }
            if (inDegreeDistribution == null) {
                out.printf(Locale.ROOT, "%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                        itr, rho0, lambdaDirected, lambdaNondirected, mu, result.times.get(result.times.size() - 1),
                        result.initialAdoptedTime, result.finalAdoptedTime, result.A.get(result.A.size() - 1), result.R.get(result.R.size() - 1),
                        result.Phi.get(result.Phi.size() - 1));
            } else {
                out.printf(Locale.ROOT, "%d,%s,%s,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d,%d%n",
                        itr, inDegreeDistribution, outDegreeDistribution, rho0, lambdaDirected,
                        lambdaNondirected, mu, result.times.get(result.times.size() - 1), result.initialAdoptedTime,
                        result.finalAdoptedTime, result.A.get(result.A.size() - 1), result.R.get(result.R.size() - 1), result.Phi.get(result.Phi.size() - 1));
            }
        }
    }
}
