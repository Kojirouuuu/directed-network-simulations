package sim.simulation;

import sim.utils.PathsEx;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.List;
import java.util.Locale;

/**
 * SAR（Susceptible-Adopted-Recovered）シミュレーション結果。
 * 可視化や再利用しやすいCSV出力機能を提供する。
 */
public final class SARResult {
    public final int n; // 頂点数
    public final List<Double> times; // 時刻のリスト

    public final List<Integer> S; // 各時刻の Susceptible 数
    public final List<Integer> A; // 各時刻の Adopted 数
    public final List<Integer> R; // 各時刻の Recovered 数

    public final double initialAdoptedTime; // 初期採用時刻
    public final double finalAdoptedTime; // 最終採用時刻

    public final double[] tInfect; // 各ノードの感染成立時刻（未感染は NaN）
    public final double[] tRecover; // 各ノードの回復成立時刻（未回復は NaN）

    SARResult(int n, List<Double> times, List<Integer> S, List<Integer> A, List<Integer> R,
            double initialAdoptedTime, double finalAdoptedTime, double[] tInfect, double[] tRecover) {
        this.n = n;
        this.times = times;
        this.S = S;
        this.A = A;
        this.R = R;
        this.initialAdoptedTime = initialAdoptedTime;
        this.finalAdoptedTime = finalAdoptedTime;
        this.tInfect = tInfect;
        this.tRecover = tRecover;
    }

    /**
     * 集計時系列CSV（itr,alpha,beta,lambda,rho0,time,A,R）を追記モードで出力する。
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
                out.println("itr,rho0,lambdaDirected,lambdaNondirected,mu,time,A,R");
            }
            for (int i = 0; i < times.size(); i++) {
                out.printf(Locale.ROOT, "%d,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d%n",
                        itr, rho0, lambdaDirected, lambdaNondirected, mu, times.get(i), A.get(i), R.get(i));
            }
        }
    }

    /**
     * 最終状態CSV（itr,alpha,beta,lambda,rho0,time,initialAdoptedTime,finalAdoptedTime,A,R）を追記モードで出力する。
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
                out.println("itr,rho0,lambdaDirected,lambdaNondirected,mu,time,initialAdoptedTime,finalAdoptedTime,A,R");
            }
            out.printf(Locale.ROOT, "%d,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d,%d%n",
                    itr, rho0, lambdaDirected, lambdaNondirected, mu, times.get(times.size() - 1), initialAdoptedTime,
                    finalAdoptedTime, A.get(A.size() - 1), R.get(R.size() - 1));
        }
    }
}
