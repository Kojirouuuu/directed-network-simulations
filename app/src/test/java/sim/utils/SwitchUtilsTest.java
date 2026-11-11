package sim.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Path;

import org.junit.jupiter.api.Test;

class SwitchUtilsTest {

    @Test
    void twoPointParametersGenerateGraphAndPreservePathsAcrossLocales() {
        var originalLocale = java.util.Locale.getDefault();
        try {
            java.util.Locale.setDefault(java.util.Locale.GERMANY);
            for (double p : new double[] {0.5, 0.125, 0.12500000000000003, 1e-7}) {
                var params = SwitchUtils.GraphGeneratorParams.builder("TwoPointDirectedCM", 10_000_000)
                        .k1(2).k2(8).p(p).r(-0.0).build();
                String decimal = p == 1e-7 ? "0.0000001" : Double.toString(p);
                var expected = Path.of("TwoPointDirectedCM/N=10000000/k1=2/k2=8/p=" + decimal + "/r=0");
                assertEquals(expected, SwitchUtils.buildNetworkPath(params));
                assertEquals(expected.resolve("randomization=none"),
                        SwitchUtils.appendRandomizationPath(expected, false));
                assertEquals(expected.resolve("randomization=edge-swap"),
                        SwitchUtils.appendRandomizationPath(expected, true));
            }
            var params = SwitchUtils.GraphGeneratorParams.builder("TwoPointDirectedCM", 100)
                    .k1(2).k2(8).p(0.5).r(0.5).build();
            assertEquals(500, SwitchUtils.generateGraph(params, 42).m);
            assertEquals(Path.of("TwoPointDirectedCM/N=200/k1=2/k2=8/p=0.5/r=0.5"),
                    SwitchUtils.buildNetworkPath(params.withN(200)));
        } finally {
            java.util.Locale.setDefault(originalLocale);
        }
    }

    @Test
    void twoPointDispatchRequiresAllParameters() {
        var missing = SwitchUtils.GraphGeneratorParams.builder("TwoPointDirectedCM", 100).build();
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.generateGraph(missing, 42));
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.buildNetworkPath(missing));
    }

    @Test
    void baFactoryAssignsAttachmentParametersAndMatchesDirectGenerator() {
        var params = SwitchUtils.GraphGeneratorParams.forBA(20, 3, 2);
        org.junit.jupiter.api.Assertions.assertEquals(3, params.m0());
        org.junit.jupiter.api.Assertions.assertEquals(2, params.m());
        org.junit.jupiter.api.Assertions.assertNull(params.kuMin());
        var graph = SwitchUtils.generateGraph(params, 42L);
        var direct = sim.network.topology.undirected.BA.generate("BA", 20, 3, 2, 42L);
        org.junit.jupiter.api.Assertions.assertEquals(direct.m, graph.m);
        for (int i = 0; i < graph.m; i++) {
            org.junit.jupiter.api.Assertions.assertEquals(direct.getOutNeighbor(i), graph.getOutNeighbor(i));
        }
    }

    @Test
    void namedParametersPreservePathsAndAllowActualGraphSize() {
        var params = SwitchUtils.GraphGeneratorParams.builder("BA", 20).m0(3).m(2).build();
        org.junit.jupiter.api.Assertions.assertEquals(java.nio.file.Path.of("BA/N=20/m0=3/m=2"),
                SwitchUtils.buildNetworkPath(params));
        org.junit.jupiter.api.Assertions.assertEquals(java.nio.file.Path.of("BA/N=40/m0=3/m=2"),
                SwitchUtils.buildNetworkPath(params.withN(40)));
        org.junit.jupiter.api.Assertions.assertEquals(20, params.N());
    }

    @Test
    void buildsAndGeneratesDirectedEr() {
        Path path = SwitchUtils.buildNetworkPath(
                "DirectedER", 1_000,
                null, 6.0, null, null, null, null,
                null, null, null, null, null, null,
                null, null, null, null, null);

        assertEquals(Path.of("DirectedER/N=1000/kAve=6.00"), path);

        var graph = SwitchUtils.generateGraph(
                SwitchUtils.GraphGeneratorParams.forDirectedER(1_000, 6.0), 42L);
        assertEquals("DirectedER", graph.name);
        assertEquals(1_000, graph.n);
        for (int edge = 0; edge < graph.m; edge++) {
            assertFalse(graph.isOutUndirected(edge));
        }
    }

    @Test
    void directedErRequiresMeanDegree() {
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.buildNetworkPath(
                "DirectedER", 100,
                null, null, null, null, null, null,
                null, null, null, null, null, null,
                null, null, null, null, null));

        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.generateGraph(
                SwitchUtils.GraphGeneratorParams.forDirectedER(100, null), 42L));
    }

    @Test
    void buildsDegreeDistributionSarPath() {
        Path path = SwitchUtils.buildNetworkPath(
                "DegreeDistributionSAR", 500_000,
                null, null, null, null, null, null,
                5, 707, null, null, null, null,
                2.5, null, null, null, null);

        assertEquals(Path.of("DegreeDistributionSAR/N=500000/gamma=2.50/kdMin=5/kdMax=707"), path);
    }

    @Test
    void appendsRandomizationPath() {
        Path networkPath = Path.of("DirectedER/N=1000/kAve=6.00");

        assertEquals(
                Path.of("DirectedER/N=1000/kAve=6.00/randomization=edge-swap"),
                SwitchUtils.appendRandomizationPath(networkPath, true));
        assertEquals(
                Path.of("DirectedER/N=1000/kAve=6.00/randomization=none"),
                SwitchUtils.appendRandomizationPath(networkPath, false));
    }

    @Test
    void buildsStandardEdgeListPath() {
        Path networkPath = Path.of("DirectedER/N=1000/kAve=6.00");

        assertEquals(
                Path.of("out/edgelist/DirectedER/N=1000/kAve=6.00/randomization=edge-swap/2.csv"),
                SwitchUtils.buildEdgeListPath(networkPath, true, 2));
        assertEquals(
                Path.of("out/edgelist/DirectedER/N=1000/kAve=6.00/randomization=none/2.csv"),
                SwitchUtils.buildEdgeListPath(networkPath, false, 2));
    }

    @Test
    void degreeDistributionSarPathRequiresItsParameters() {
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.buildNetworkPath(
                "DegreeDistributionSAR", 100,
                null, null, null, null, null, null,
                null, 10, null, null, null, null,
                2.5, null, null, null, null));
    }
}
