package sim.simulation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import sim.simulation.SPIRGillespieSimulator.ModelMode;

/** Population trajectories and per-node transition times. Unreached times are NaN. */
public final class SPIRResult {
    public final int n;
    public final List<Double> times = new ArrayList<>();
    public final List<Integer> S = new ArrayList<>(), P = new ArrayList<>(), I = new ArrayList<>(),
                               R = new ArrayList<>();
    public final double[] tInfect, tOnset, tRecover;
    /** Observer membership and isolation are independent of the S/P/I/R health state. */
    public final boolean[] observers, quarantined;
    public final double[] tQuarantine;
    public final List<Integer> Q = new ArrayList<>();
    int quarantinedCount, observerCount;

    int recovered, finalS, finalP, finalI;
    double finalTime;
    public int observerCount() {
        return observerCount;
    }

    SPIRResult(int n) {
        this.n = n;
        tInfect = new double[n];
        tOnset = new double[n];
        tRecover = new double[n];
        observers = new boolean[n];
        quarantined = new boolean[n];
        tQuarantine = new double[n];
        Arrays.fill(tInfect, Double.NaN);
        Arrays.fill(tOnset, Double.NaN);
        Arrays.fill(tRecover, Double.NaN);
        Arrays.fill(tQuarantine, Double.NaN);
    }

    public void writeCsv(Path path, int itr, double rho0, ModelMode mode, double lambda1, double lambda2,
            double mu1, double mu2, boolean finalOnly, boolean append) throws IOException {
        SPIRResultCsvWriter.write(this, path, itr, rho0, mode, lambda1, lambda2, mu1, mu2, finalOnly, append);
    }

    /** Observer-aware CSV appends phi, observer_count and cumulative quarantined population Q. */
    public void writeCsv(Path path, int itr, double rho0, double phi, ModelMode mode, double lambda1,
            double lambda2, double mu1, double mu2, boolean finalOnly, boolean append) throws IOException {
        SPIRResultCsvWriter.writeWithObservers(
                this, path, itr, rho0, phi, mode, lambda1, lambda2, mu1, mu2, finalOnly, append);
    }
    public static CsvBatch openCsv(Path path) throws IOException {
        return new CsvBatch(path);
    }

    public static final class CsvBatch implements AutoCloseable {
        private final SPIRResultCsvWriter writer;
        private CsvBatch(Path path) throws IOException {
            writer = new SPIRResultCsvWriter(path, true, true);
        }
        public void write(SPIRResult result, int itr, double rho, double phi, ModelMode mode, double lambda1,
                double lambda2, double mu1, double mu2, boolean last) throws IOException {
            if (!Double.isFinite(phi) || phi < 0 || phi > 1)
                throw new IllegalArgumentException("invalid phi");
            writer.writeResult(result, itr, rho, phi, mode, lambda1, lambda2, mu1, mu2, last);
        }
        @Override
        public void close() throws IOException {
            writer.close();
        }
    }
}
