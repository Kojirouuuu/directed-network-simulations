package sim.simulation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

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
    public final List<Integer> Phi; // あと一回の伝播で採用されるノード数

    public final double initialAdoptedTime; // 初期採用時刻
    public final double finalAdoptedTime; // 最終採用時刻

    public final double[] tInfect; // 各ノードの感染成立時刻（未感染は NaN）
    public final double[] tRecover; // 各ノードの回復成立時刻（未回復は NaN）

    SARResult(int n, List<Double> times, List<Integer> S, List<Integer> A, List<Integer> R,
            List<Integer> Phi, double initialAdoptedTime, double finalAdoptedTime, double[] tInfect, double[] tRecover) {
        this.n = n;
        this.times = times;
        this.S = S;
        this.A = A;
        this.R = R;
        this.Phi = Phi;
        this.initialAdoptedTime = initialAdoptedTime;
        this.finalAdoptedTime = finalAdoptedTime;
        this.tInfect = tInfect;
        this.tRecover = tRecover;
    }

    public void writeTimeSeriesCsv(Path path, int itr, double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        new SARResultCsvWriter(this).writeTimeSeriesCsv(path, itr, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    public void writeTimeSeriesCsv(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        new SARResultCsvWriter(this).writeTimeSeriesCsv(path, itr, inDegreeDistribution, outDegreeDistribution, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    public void writeFinalStateCsv(Path path, int itr, double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        new SARResultCsvWriter(this).writeFinalStateCsv(path, itr, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

    public void writeFinalStateCsv(Path path, int itr,
            String inDegreeDistribution, String outDegreeDistribution,
            double rho0, double lambdaDirected, double lambdaNondirected, double mu,
            boolean append) throws IOException {
        new SARResultCsvWriter(this).writeFinalStateCsv(path, itr, inDegreeDistribution, outDegreeDistribution, rho0, lambdaDirected, lambdaNondirected, mu, append);
    }

}
