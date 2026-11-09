package sim;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import sim.simulation.SPIRGillespieSimulator.ModelMode;

class SPIRTest {
    @org.junit.jupiter.api.io.TempDir java.nio.file.Path temp;

    /** Independent of the editable experiment defaults in SPIR.java. */
    private SPIR.SimulationConfig testConfig() {
        SPIR.SimulationConfig config = new SPIR.SimulationConfig();
        config.outputRoot = temp;
        config.mode = ModelMode.SPIR;
        config.engine = sim.simulation.SPIRSimulator.Engine.GILLESPIE;
        config.optionPath = "test-fixture";
        config.network = sim.utils.SwitchUtils.GraphGeneratorParams.builder("ER", 20).kuAve(4.0).build();
        config.loadFromEdgeList = false;
        config.writeEdgeList = false;
        config.isFinal = true;
        config.randomizationMode = SAR.RandomizationMode.NONE;
        config.sizeMultiplier = 1;
        config.batchSize = 1;
        config.itrs = 1;
        config.tMax = 2;
        config.dt = 0.1;
        config.phi = 0;
        config.rho0List = new double[] {0.1};
        config.equalLambda = false;
        config.lambda1List = new double[] {0.2};
        config.lambda2List = new double[] {0.2};
        config.mu1List = new double[] {1};
        config.mu2List = new double[] {1};
        return config;
    }


    @org.junit.jupiter.params.ParameterizedTest
    @org.junit.jupiter.params.provider.ValueSource(booleans = {false, true})
    void batchRunnerWritesCsvAndSummaryForBothModes(boolean equalLambda) throws Exception {
        for (ModelMode mode : ModelMode.values()) {
            SPIR.SimulationConfig config = testConfig();
            config.mode = mode;
            config.equalLambda = equalLambda;
            int conditions = mode == ModelMode.SIR ? 4 : equalLambda ? 8 : 16;
            config.outputRoot = temp;
            config.optionPath = "spir-example";
            config.network = sim.utils.SwitchUtils.GraphGeneratorParams.builder("BA", 10).m0(2).m(2).build();
            config.lambda1List = new double[] {0, 1};
            config.lambda2List = new double[] {0, 1};
            config.mu1List = new double[] {0, 1};
            config.mu2List = new double[] {0, 1};
            config.rho0List = new double[] {0};
            config.batchSize = 2;
            config.itrs = 2;
            config.tMax = 2;
            SPIR.run(config);
            java.nio.file.Path root =
                    temp.resolve("spir/spir-example/model=" + mode.name().toLowerCase(java.util.Locale.ROOT));
            try (var files = java.nio.file.Files.walk(root)) {
                var csvs = files.filter(p -> p.toString().endsWith(".csv")).toList();
                assertEquals(2, csvs.size());
                for (var csv : csvs) {
                    var lines = java.nio.file.Files.readAllLines(csv);
                    assertEquals(conditions * 2 + 1, lines.size());
                    if (equalLambda && mode == ModelMode.SPIR) {
                        for (String line : lines.subList(1, lines.size())) {
                            String[] cells = line.split(",", -1);
                            assertEquals(cells[3], cells[4]);
                        }
                    }
                }
            }
            try (var files = java.nio.file.Files.walk(root)) {
                var json = files.filter(p -> p.toString().endsWith(".json")).findFirst().orElseThrow();
                var summary = com.google.gson.JsonParser.parseString(java.nio.file.Files.readString(json))
                                      .getAsJsonObject();
                assertEquals("succeeded", summary.get("status").getAsString());
                assertEquals(conditions * 4,
                        summary.getAsJsonObject("execution").get("completedSimulationCount").getAsLong());
                assertEquals(2, summary.getAsJsonArray("batches").size());
                assertEquals(equalLambda, summary.getAsJsonObject("parameters").get("equalLambda").getAsBoolean());
            }
        }
    }

    @Test
    void sweepsAllRatesButOnlyEffectiveSirRates() {
        double[] values = {0, 1};
        assertEquals(16, SPIR.rates(ModelMode.SPIR, values, values, values, values).size());
        var sir = SPIR.rates(ModelMode.SIR, values, values, values, values);
        assertEquals(4, sir.size());
        assertEquals(4, sir.stream().distinct().count());
    }

    @Test
    void equalLambdaUsesLambda2ValuesAndStillSweepsBothTransitionRates() {
        double[] lambda1 = {9, 10};
        double[] lambda2 = {0, 0.2, 0.4};
        double[] mu = {0.5, 1};
        var rates = SPIR.rates(ModelMode.SPIR, lambda1, lambda2, mu, mu, true);
        assertEquals(12, rates.size());
        assertEquals(12, rates.stream().distinct().count());
        assertTrue(rates.stream().allMatch(r -> r.lambda1() == r.lambda2()));
        assertEquals(java.util.Set.of(0.0, 0.2, 0.4), rates.stream()
                .map(SPIR.Rates::lambda1).collect(java.util.stream.Collectors.toSet()));
        assertEquals(SPIR.rates(ModelMode.SIR, lambda1, lambda2, mu, mu),
                SPIR.rates(ModelMode.SIR, lambda1, lambda2, mu, mu, true));
    }
    @Test
    void modesHaveSeparateOutputPathsAndUseActualN() {
        SPIR.SimulationConfig config = testConfig();
        var graph = sim.network.DirectedGraph.fromEdgeListWithUndirectedFlag(
                "test", 3, new int[0], new int[0], new boolean[0]);
        String path = SPIR.outputDirectory(graph, config).toString();
        assertTrue(path.contains("model=spir"));
        assertTrue(path.contains("N=3"));
        assertFalse(path.contains("threshold="));
    }
    @Test
    void samplesExactObserverFractionAndInitialInfectionsOutsideObservers() {
        var sample = SPIR.samplePopulation(100_000, 0.001, 0.25, 2007, 424242);
        assertEquals(25_000, sample.observers().length);
        assertEquals(100, sample.initial().length);
        boolean[] observers = new boolean[100_000];
        for (int node : sample.observers()) {
            assertFalse(observers[node]);
            observers[node] = true;
        }
        for (int node : sample.initial()) assertFalse(observers[node]);
        var repeated = SPIR.samplePopulation(100_000, 0.001, 0.25, 2007, 424242);
        assertArrayEquals(sample.observers(), repeated.observers());
        assertArrayEquals(sample.initial(), repeated.initial());
        var none = SPIR.samplePopulation(100, 0.1, 0, 2007, 424242);
        assertArrayEquals(sim.simulation.InitialAdopters.sample(100, 0.1, 2007), none.initial());
        assertEquals(0, none.observers().length);
    }

    @Test
    void validatesPhiAndAvailabilityOfInitialNonObservers() {
        for (double phi : new double[] {-0.1, 1.1, Double.NaN, Double.POSITIVE_INFINITY}) {
            assertThrows(IllegalArgumentException.class, () -> SPIR.samplePopulation(100, 0.1, phi, 1, 2));
        }
        assertThrows(IllegalArgumentException.class, () -> SPIR.samplePopulation(100, 0.01, 1, 1, 2));
        assertThrows(IllegalArgumentException.class, () -> SPIR.samplePopulation(100, 0.1, 0.95, 1, 2));
    }

    @Test
    void observerBatchWritesPhiAndObserverCountsToCsvAndSummary() throws Exception {
        SPIR.SimulationConfig config = testConfig();
        config.outputRoot = temp;
        config.optionPath = "observer-batch";
        config.network = sim.utils.SwitchUtils.GraphGeneratorParams.builder("ER", 20).kuAve(4.0).build();
        config.phi = 0.25;
        config.rho0List = new double[] {0.1};
        config.lambda1List = new double[] {0.2};
        config.lambda2List = new double[] {0.2};
        config.mu1List = new double[] {1};
        config.mu2List = new double[] {1};
        config.batchSize = 1;
        config.itrs = 1;
        config.tMax = 100;
        SPIR.run(config);
        try (var files = java.nio.file.Files.walk(temp.resolve("spir/observer-batch"))) {
            var paths = files.filter(java.nio.file.Files::isRegularFile).toList();
            var csv = paths.stream().filter(p -> p.toString().endsWith(".csv")).findFirst().orElseThrow();
            assertTrue(csv.toString().contains("phi=0.25"));
            var lines = java.nio.file.Files.readAllLines(csv);
            assertTrue(lines.get(0).endsWith(",phi,observer_count,Q"));
            assertEquals("0.250000000", lines.get(1).split(",")[12]);
            assertEquals("5", lines.get(1).split(",")[13]);
            var json = paths.stream().filter(p -> p.toString().endsWith(".json")).findFirst().orElseThrow();
            var summary = com.google.gson.JsonParser.parseString(java.nio.file.Files.readString(json)).getAsJsonObject();
            assertEquals(0.25, summary.getAsJsonObject("parameters").get("phi").getAsDouble());
            assertTrue(summary.getAsJsonObject("parameters").getAsJsonObject("seeds").has("observer"));
        }
    }

}
