package sim.simulation;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class SARResultCsvTest {

    @TempDir
    Path tempDir;

    @Test
    void existingFinalStateFormatRemainsUnchanged() throws Exception {
        Path path = tempDir.resolve("legacy.csv");
        result().writeFinalStateCsv(path, 3, 0.1, 2.0, 0.0, 1.0, true);

        List<String> lines = Files.readAllLines(path);
        assertEquals("itr,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi",
                lines.get(0));
        assertEquals(11, lines.get(1).split(",").length);
    }

    @Test
    void distributionLabelsAreWrittenToFinalStateCsv() throws Exception {
        Path path = tempDir.resolve("comparison.csv");
        result().writeFinalStateCsv(path, 3, "Pow", "Poi", 0.1, 2.0, 0.0, 1.0, true);

        List<String> lines = Files.readAllLines(path);
        assertEquals("itr,in_degree_distribution,out_degree_distribution,rho_0,lambda_d,lambda_u,mu,time,initial_adopted_time,final_adopted_time,A,R,Phi",
                lines.get(0));
        assertEquals("Pow", lines.get(1).split(",")[1]);
        assertEquals("Poi", lines.get(1).split(",")[2]);
        assertEquals(13, lines.get(1).split(",").length);
    }

    @Test
    void distributionLabelsAreWrittenToEveryTimeSeriesRow() throws Exception {
        Path path = tempDir.resolve("timeseries.csv");
        result().writeTimeSeriesCsv(path, 3, "Poi", "Pow", 0.1, 2.0, 0.0, 1.0, true);

        List<String> lines = Files.readAllLines(path);
        assertEquals("itr,in_degree_distribution,out_degree_distribution,rho_0,lambda_d,lambda_u,mu,time,A,R,Phi",
                lines.get(0));
        assertEquals(3, lines.size());
        for (String row : lines.subList(1, lines.size())) {
            assertEquals("Poi", row.split(",")[1]);
            assertEquals("Pow", row.split(",")[2]);
            assertEquals(11, row.split(",").length);
        }
    }

    private static SARResult result() {
        return new SARResult(
                10,
                List.of(0.0, 1.5),
                List.of(9, 4),
                List.of(1, 2),
                List.of(0, 4),
                List.of(2, 1),
                0.25,
                1.25,
                new double[10],
                new double[10]);
    }
}
