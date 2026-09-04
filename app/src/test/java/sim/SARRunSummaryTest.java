package sim;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SARRunSummaryTest {
    @TempDir
    Path tempDir;

    @Test
    void writesSuccessfulRunWithParametersCountsAndBatchMetadata() throws IOException {
        SARRunSummary summary = newSummary(2);
        Instant batchStartedAt = Instant.parse("2026-09-03T01:00:00Z");
        Path resultPath = tempDir.resolve("results_00.csv");

        summary.recordBatchStarted(0, batchStartedAt);
        summary.registerBatchGraph(0, 100, 250, tempDir);
        JsonObject running = read(summary.summaryPath());
        assertEquals("running", running.get("status").getAsString());
        assertTrue(running.getAsJsonObject("timing").get("endedAt").isJsonNull());
        summary.recordResultPath(0, resultPath);
        summary.recordSimulationCompleted();
        summary.recordSimulationCompleted();
        summary.recordBatchCompleted(0, batchStartedAt.plusSeconds(3));
        summary.finishSucceeded();

        JsonObject root = read(summary.summaryPath());
        assertEquals(1, root.get("schemaVersion").getAsInt());
        assertEquals("succeeded", root.get("status").getAsString());
        assertFalse(root.getAsJsonObject("timing").get("endedAt").isJsonNull());

        JsonObject execution = root.getAsJsonObject("execution");
        assertEquals(2, execution.get("plannedSimulationCount").getAsLong());
        assertEquals(2, execution.get("completedSimulationCount").getAsLong());
        assertEquals(0, execution.get("failedSimulationCount").getAsLong());
        assertEquals(1, execution.get("completedBatchCount").getAsLong());

        JsonObject parameters = root.getAsJsonObject("parameters");
        assertEquals(2, parameters.getAsJsonArray("lambdaDirectedValues").size());

        JsonObject batch = root.getAsJsonArray("batches").get(0).getAsJsonObject();
        assertEquals("succeeded", batch.get("status").getAsString());
        assertEquals(100, batch.get("actualN").getAsInt());
        assertEquals(250, batch.get("edgeCount").getAsInt());
        assertEquals(3_000, batch.get("durationMillis").getAsLong());
        assertEquals(resultPath.toString(), batch.get("resultCsv").getAsString());
        assertTrue(root.getAsJsonArray("errors").isEmpty());

        JsonObject system = root.getAsJsonObject("system");
        assertTrue(system.has("osName"));
        assertTrue(system.has("osArchitecture"));
        assertTrue(system.has("javaVersion"));
        assertTrue(system.get("availableProcessors").getAsInt() > 0);

        try (var files = Files.list(tempDir)) {
            assertTrue(files.noneMatch(path -> path.getFileName().toString().endsWith(".tmp")));
        }
    }

    @Test
    void writesFailedRunWithSimulationContextAndStackTrace() throws IOException {
        SARRunSummary summary = newSummary(4);
        Instant batchStartedAt = Instant.parse("2026-09-03T01:00:00Z");
        IllegalStateException failure = new IllegalStateException("simulation failed");

        summary.recordBatchStarted(2, batchStartedAt);
        summary.registerBatchGraph(2, 50, 75, tempDir);
        summary.recordSimulationCompleted();
        summary.recordSimulationFailed();
        summary.recordError("simulation",
                SARRunSummary.ErrorContext.forSimulation(2, 3, 0.01, 0.2, 0.0, 1.0),
                failure);
        summary.recordBatchFailed(2, batchStartedAt.plusSeconds(1));
        summary.finishFailed();

        JsonObject root = read(summary.summaryPath());
        assertEquals("failed", root.get("status").getAsString());
        assertEquals(1, root.getAsJsonObject("execution").get("completedSimulationCount").getAsLong());
        assertEquals(1, root.getAsJsonObject("execution").get("failedSimulationCount").getAsLong());
        assertEquals(1, root.getAsJsonObject("execution").get("failedBatchCount").getAsLong());

        JsonObject error = root.getAsJsonArray("errors").get(0).getAsJsonObject();
        assertEquals("simulation", error.get("stage").getAsString());
        assertEquals(2, error.get("batchIndex").getAsInt());
        assertEquals(3, error.get("iteration").getAsInt());
        assertEquals(0.2, error.get("lambdaDirected").getAsDouble());
        assertEquals(IllegalStateException.class.getName(), error.get("exceptionClass").getAsString());
        assertTrue(error.get("stackTrace").getAsString().contains("simulation failed"));
    }

    @Test
    void usesUniqueFilesForRunsInTheSameOutputDirectory() throws IOException {
        SARRunSummary first = newSummary(0);
        SARRunSummary second = newSummary(0);

        first.recordBatchStarted(0, Instant.now());
        first.registerBatchGraph(0, 10, 20, tempDir);
        first.finishSucceeded();
        second.recordBatchStarted(0, Instant.now());
        second.registerBatchGraph(0, 10, 20, tempDir);
        second.finishInterruptedBestEffort();

        assertNotEquals(first.summaryPath(), second.summaryPath());
        assertEquals("interrupted", read(second.summaryPath()).get("status").getAsString());
    }

    @Test
    void doesNotWriteSummaryBeforeAnOutputDirectoryIsKnown() throws IOException {
        SARRunSummary summary = newSummary(1);
        summary.recordUnhandledErrorIfAbsent(new IllegalArgumentException("early failure"));
        summary.finishFailed();

        assertNull(summary.summaryPath());
        try (var files = Files.list(tempDir)) {
            assertEquals(0, files.count());
        }
    }

    private static SARRunSummary newSummary(long plannedSimulationCount) {
        Map<String, Object> network = new LinkedHashMap<>();
        network.put("name", "test-network");
        network.put("configuredN", 100);

        Map<String, Object> parameters = new LinkedHashMap<>();
        parameters.put("lambdaDirectedValues", new double[] { 0.0, 0.2 });

        Map<String, Object> execution = new LinkedHashMap<>();
        execution.put("batchSize", 1);
        execution.put("iterationsPerBatch", 2);
        return new SARRunSummary(plannedSimulationCount, network, parameters, execution);
    }

    private static JsonObject read(Path path) throws IOException {
        return JsonParser.parseString(Files.readString(path)).getAsJsonObject();
    }
}
