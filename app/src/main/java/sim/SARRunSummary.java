package sim;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.AtomicMoveNotSupportedException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.Duration;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.UUID;
import java.util.concurrent.atomic.AtomicLong;

/** Records one SAR invocation in a machine-readable JSON summary. */
final class SARRunSummary {
    private static final int SCHEMA_VERSION = 1;
    private static final DateTimeFormatter RUN_ID_TIME = DateTimeFormatter
            .ofPattern("yyyyMMdd'T'HHmmssSSS'Z'")
            .withZone(ZoneOffset.UTC);
    private static final Gson GSON = new GsonBuilder()
            .serializeNulls()
            .setPrettyPrinting()
            .create();

    record ErrorContext(Integer batchIndex, Integer iteration,
            Double rho0, Double lambdaDirected, Double lambdaNondirected, Double mu) {
        static ErrorContext forBatch(int batchIndex) {
            return new ErrorContext(batchIndex, null, null, null, null, null);
        }

        static ErrorContext forSimulation(int batchIndex, int iteration,
                double rho0, double lambdaDirected, double lambdaNondirected, double mu) {
            return new ErrorContext(batchIndex, iteration, rho0, lambdaDirected, lambdaNondirected, mu);
        }
    }

    private static final class BatchRecord {
        private final int batchIndex;
        private final String startedAt;
        private String endedAt;
        private Long durationMillis;
        private String status = "running";
        private Integer actualN;
        private Integer edgeCount;
        private String resultCsv;

        private BatchRecord(int batchIndex, Instant startedAt) {
            this.batchIndex = batchIndex;
            this.startedAt = startedAt.toString();
        }

        private Map<String, Object> snapshot() {
            Map<String, Object> value = new LinkedHashMap<>();
            value.put("batchIndex", batchIndex);
            value.put("status", status);
            value.put("startedAt", startedAt);
            value.put("endedAt", endedAt);
            value.put("durationMillis", durationMillis);
            value.put("actualN", actualN);
            value.put("edgeCount", edgeCount);
            value.put("resultCsv", resultCsv);
            return value;
        }
    }

    private final String runId;
    private final Instant startedAt;
    private final long plannedSimulationCount;
    private final Map<String, Object> network;
    private final Map<String, Object> parameters;
    private final Map<String, Object> executionConfiguration;
    private final Map<String, Object> system;
    private final Map<Integer, BatchRecord> batches = new TreeMap<>();
    private final List<Map<String, Object>> errors = new ArrayList<>();
    private final AtomicLong completedSimulationCount = new AtomicLong();
    private final AtomicLong failedSimulationCount = new AtomicLong();
    private final AtomicLong completedBatchCount = new AtomicLong();
    private final AtomicLong failedBatchCount = new AtomicLong();

    private String status = "running";
    private Instant endedAt;
    private Path summaryPath;
    private boolean finalized;

    SARRunSummary(long plannedSimulationCount,
            Map<String, Object> network,
            Map<String, Object> parameters,
            Map<String, Object> executionConfiguration) {
        this.startedAt = Instant.now();
        this.runId = RUN_ID_TIME.format(startedAt) + "_" + UUID.randomUUID().toString().substring(0, 8);
        this.plannedSimulationCount = plannedSimulationCount;
        this.network = new LinkedHashMap<>(network);
        this.parameters = new LinkedHashMap<>(parameters);
        this.executionConfiguration = new LinkedHashMap<>(executionConfiguration);
        this.system = minimalSystemInformation();
    }

    synchronized void recordBatchStarted(int batchIndex, Instant batchStartedAt) {
        batches.putIfAbsent(batchIndex, new BatchRecord(batchIndex, batchStartedAt));
    }

    synchronized void registerBatchGraph(int batchIndex, int actualN, int edgeCount, Path outputDirectory) {
        BatchRecord batch = requireBatch(batchIndex);
        batch.actualN = actualN;
        batch.edgeCount = edgeCount;

        if (summaryPath == null) {
            try {
                Files.createDirectories(outputDirectory);
                summaryPath = outputDirectory.resolve("summary_" + runId + ".json");
                persist();
            } catch (IOException e) {
                throw new UncheckedIOException("Failed to create SAR run summary in " + outputDirectory, e);
            }
        }
    }

    synchronized void recordResultPath(int batchIndex, Path resultPath) {
        requireBatch(batchIndex).resultCsv = resultPath.toString();
    }

    void recordSimulationCompleted() {
        completedSimulationCount.incrementAndGet();
    }

    void recordSimulationFailed() {
        failedSimulationCount.incrementAndGet();
    }

    synchronized void recordBatchCompleted(int batchIndex, Instant batchEndedAt) {
        BatchRecord batch = requireBatch(batchIndex);
        if (!"succeeded".equals(batch.status)) {
            completedBatchCount.incrementAndGet();
        }
        finishBatch(batch, "succeeded", batchEndedAt);
    }

    synchronized void recordBatchFailed(int batchIndex, Instant batchEndedAt) {
        BatchRecord batch = requireBatch(batchIndex);
        if (!"failed".equals(batch.status)) {
            failedBatchCount.incrementAndGet();
        }
        finishBatch(batch, "failed", batchEndedAt);
    }

    synchronized void recordError(String stage, ErrorContext context, Throwable error) {
        Map<String, Object> value = new LinkedHashMap<>();
        value.put("occurredAt", Instant.now().toString());
        value.put("stage", stage);
        value.put("batchIndex", context == null ? null : context.batchIndex());
        value.put("iteration", context == null ? null : context.iteration());
        value.put("rho0", context == null ? null : context.rho0());
        value.put("lambdaDirected", context == null ? null : context.lambdaDirected());
        value.put("lambdaNondirected", context == null ? null : context.lambdaNondirected());
        value.put("mu", context == null ? null : context.mu());
        value.put("exceptionClass", error.getClass().getName());
        value.put("message", error.getMessage());
        value.put("stackTrace", stackTrace(error));
        errors.add(value);
    }

    synchronized void recordUnhandledErrorIfAbsent(Throwable error) {
        if (errors.isEmpty()) {
            recordError("run", null, error);
        }
    }

    synchronized void finishSucceeded() throws IOException {
        finish("succeeded");
    }

    synchronized void finishFailed() throws IOException {
        finish("failed");
    }

    synchronized void finishInterruptedBestEffort() {
        if (finalized || summaryPath == null) {
            return;
        }
        try {
            finish("interrupted");
        } catch (IOException e) {
            System.err.println("Failed to finalize interrupted SAR run summary: " + e.getMessage());
        }
    }

    synchronized Path summaryPath() {
        return summaryPath;
    }

    private void finish(String finalStatus) throws IOException {
        if (finalized) {
            return;
        }
        status = finalStatus;
        endedAt = Instant.now();
        finalized = true;
        if (summaryPath != null) {
            persist();
        }
    }

    private BatchRecord requireBatch(int batchIndex) {
        BatchRecord batch = batches.get(batchIndex);
        if (batch == null) {
            throw new IllegalStateException("Batch was not started: " + batchIndex);
        }
        return batch;
    }

    private static void finishBatch(BatchRecord batch, String status, Instant batchEndedAt) {
        batch.status = status;
        batch.endedAt = batchEndedAt.toString();
        batch.durationMillis = Duration.between(Instant.parse(batch.startedAt), batchEndedAt).toMillis();
    }

    private void persist() throws IOException {
        Path parent = summaryPath.getParent();
        Path temporary = Files.createTempFile(parent, "." + summaryPath.getFileName(), ".tmp");
        try {
            Files.writeString(temporary, GSON.toJson(snapshot()), StandardCharsets.UTF_8);
            try {
                Files.move(temporary, summaryPath,
                        StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING);
            } catch (AtomicMoveNotSupportedException e) {
                Files.move(temporary, summaryPath, StandardCopyOption.REPLACE_EXISTING);
            }
        } finally {
            Files.deleteIfExists(temporary);
        }
    }

    private Map<String, Object> snapshot() {
        Map<String, Object> value = new LinkedHashMap<>();
        value.put("schemaVersion", SCHEMA_VERSION);
        value.put("runId", runId);
        value.put("status", status);

        Map<String, Object> timing = new LinkedHashMap<>();
        timing.put("startedAt", startedAt.toString());
        timing.put("endedAt", endedAt == null ? null : endedAt.toString());
        timing.put("durationMillis", endedAt == null ? null : Duration.between(startedAt, endedAt).toMillis());
        value.put("timing", timing);

        value.put("network", new LinkedHashMap<>(network));
        value.put("parameters", new LinkedHashMap<>(parameters));

        Map<String, Object> execution = new LinkedHashMap<>(executionConfiguration);
        execution.put("plannedSimulationCount", plannedSimulationCount);
        execution.put("completedSimulationCount", completedSimulationCount.get());
        execution.put("failedSimulationCount", failedSimulationCount.get());
        execution.put("completedBatchCount", completedBatchCount.get());
        execution.put("failedBatchCount", failedBatchCount.get());
        value.put("execution", execution);

        value.put("system", new LinkedHashMap<>(system));
        value.put("batches", batches.values().stream().map(BatchRecord::snapshot).toList());
        value.put("errors", new ArrayList<>(errors));
        return value;
    }

    private static Map<String, Object> minimalSystemInformation() {
        Map<String, Object> value = new LinkedHashMap<>();
        value.put("osName", System.getProperty("os.name"));
        value.put("osVersion", System.getProperty("os.version"));
        value.put("osArchitecture", System.getProperty("os.arch"));
        value.put("javaVersion", System.getProperty("java.version"));
        value.put("availableProcessors", Runtime.getRuntime().availableProcessors());
        return value;
    }

    private static String stackTrace(Throwable error) {
        StringWriter buffer = new StringWriter();
        error.printStackTrace(new PrintWriter(buffer));
        return buffer.toString();
    }
}
