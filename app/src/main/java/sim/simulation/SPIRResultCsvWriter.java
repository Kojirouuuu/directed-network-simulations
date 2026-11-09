package sim.simulation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Locale;
import sim.simulation.SPIRGillespieSimulator.ModelMode;
import sim.utils.PathsEx;

final class SPIRResultCsvWriter implements AutoCloseable {
    private final BufferedWriter out;
    SPIRResultCsvWriter(Path path, boolean append, boolean observers) throws IOException {
        if (!append)
            path = PathsEx.resolveIndexed(path);
        if (path.getParent() != null)
            Files.createDirectories(path.getParent());
        String expected = "model,itr,rho_0,lambda1,lambda2,mu1,mu2,time,S,P,I,R"
                + (observers ? ",phi,observer_count,Q" : "");
        boolean header = !Files.exists(path) || Files.size(path) == 0;
        if (!header && append) {
            try (var reader = Files.newBufferedReader(path)) {
                if (!expected.equals(reader.readLine()))
                    throw new IOException("CSV schema mismatch: " + path);
            }
        }
        out = Files.newBufferedWriter(path, StandardOpenOption.CREATE,
                append ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING);
        if (header)
            out.write(expected + "\n");
    }
    @Override
    public void close() throws IOException {
        out.close();
    }

    static void write(SPIRResult result, Path path, int itr, double rho, ModelMode mode, double lambda1,
            double lambda2, double mu1, double mu2, boolean last, boolean append) throws IOException {
        writeInternal(result, path, itr, rho, null, mode, lambda1, lambda2, mu1, mu2, last, append);
    }

    static void writeWithObservers(SPIRResult result, Path path, int itr, double rho, double phi,
            ModelMode mode, double lambda1, double lambda2, double mu1, double mu2, boolean last,
            boolean append) throws IOException {
        if (!Double.isFinite(phi) || phi < 0 || phi > 1)
            throw new IllegalArgumentException("phi must be in [0,1]");
        writeInternal(result, path, itr, rho, phi, mode, lambda1, lambda2, mu1, mu2, last, append);
    }

    private static void writeInternal(SPIRResult result, Path path, int itr, double rho, Double phi,
            ModelMode mode, double lambda1, double lambda2, double mu1, double mu2, boolean last,
            boolean append) throws IOException {
        try (var writer = new SPIRResultCsvWriter(path, append, phi != null)) {
            writer.writeResult(result, itr, rho, phi, mode, lambda1, lambda2, mu1, mu2, last);
        }
    }

    void writeResult(SPIRResult result, int itr, double rho, Double phi, ModelMode mode, double lambda1,
            double lambda2, double mu1, double mu2, boolean last) throws IOException {
        String prefix = String.format(Locale.ROOT, "%s,%d,%.9f,%s,%.9f,%s,%.9f,",
                mode.name().toLowerCase(Locale.ROOT), itr, rho,
                mode == ModelMode.SIR ? "" : String.format(Locale.ROOT, "%.9f", lambda1), lambda2,
                mode == ModelMode.SIR ? "" : String.format(Locale.ROOT, "%.9f", mu1), mu2);
        if (last) {
            row(out, prefix, result.finalTime, result.finalS, result.finalP, result.finalI, result.recovered,
                    phi, result.observerCount(), result.quarantinedCount);
        } else {
            for (int k = 0; k < result.times.size(); k++) {
                row(out, prefix, result.times.get(k), result.S.get(k), result.P.get(k), result.I.get(k),
                        result.R.get(k), phi, result.observerCount(), result.Q.get(k));
            }
        }
    }

    private static void row(BufferedWriter out, String prefix, double time, int s, int p, int i, int r,
            Double phi, int observers, int quarantined) throws IOException {
        out.write(prefix + String.format(Locale.ROOT, "%.9f,%d,%d,%d,%d", time, s, p, i, r));
        if (phi != null)
            out.write(String.format(Locale.ROOT, ",%.9f,%d,%d", phi, observers, quarantined));
        out.newLine();
    }
}
