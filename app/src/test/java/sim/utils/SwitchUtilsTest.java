package sim.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Path;

import org.junit.jupiter.api.Test;

class SwitchUtilsTest {

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
    void degreeDistributionSarPathRequiresItsParameters() {
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.buildNetworkPath(
                "DegreeDistributionSAR", 100,
                null, null, null, null, null, null,
                null, 10, null, null, null, null,
                2.5, null, null, null, null));
    }
}
