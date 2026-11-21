package sim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class SAREdgeListPathTest {
    @TempDir
    Path tempDir;

    @Test
    void sharedOutputPathsPreserveTopologyAndUseActualSizeAndEffectiveRandomization() {
        var params = sim.utils.SwitchUtils.GraphGeneratorParams.builder("BA", 100).m0(3).m(2).build();
        assertEquals(Path.of("BA/N=80/m0=3/m=2/randomization=none"),
                SARGraphPreparation.buildNetworkOutputPath(params, 80, SAR.RandomizationMode.NONE, 1));
        for (SAR.RandomizationMode mode : SAR.RandomizationMode.values()) {
            assertEquals(SAR.appendRandomizationPath(Path.of("BA/N=80/m0=3/m=2"), mode, 1),
                    SARGraphPreparation.buildNetworkOutputPath(params, 80, mode, 1));
        }
        assertEquals(Path.of("BA/N=200/m0=3/m=2/randomization=joint-degree-cm"),
                SARGraphPreparation.buildNetworkOutputPath(params, 200, SAR.RandomizationMode.JOINT_DEGREE_CM, 2));
        var twoPoint = sim.utils.SwitchUtils.GraphGeneratorParams.builder("TwoPointDirectedCM", 100)
                .k1(2).k2(8).p(0.5).r(0.5).build();
        assertEquals(Path.of("TwoPointDirectedCM/N=200/k1=2/k2=8/p=0.5/r=0.5/randomization=none"),
                SARGraphPreparation.buildNetworkOutputPath(twoPoint, 200, SAR.RandomizationMode.NONE, 1));
    }

    @Test
    void resolvesSwapEdgeListStoredUnderActualNodeCountAndExtraDirectory() throws IOException {
        Path edgeList = tempDir.resolve(
                "rev-higgs-social/N=456626/randomization=edge-swap/reciprocity=0.010000/0.csv");
        Files.createDirectories(edgeList.getParent());
        Files.createFile(edgeList);

        Path resolved = SAR.resolveEdgeListPath(
                tempDir, Path.of("rev-higgs-social/N=500000"), true, 0);

        assertEquals(edgeList, resolved);
    }

    @Test
    void prefersTheStandardPath() throws IOException {
        Path standard = tempDir.resolve("DirectedER/N=1000/kAve=6.00/randomization=edge-swap/2.csv");
        Files.createDirectories(standard.getParent());
        Files.createFile(standard);

        assertEquals(standard, SAR.resolveEdgeListPath(
                tempDir, Path.of("DirectedER/N=1000/kAve=6.00"), true, 2));
    }

    @Test
    void rejectsAmbiguousFallbackMatches() throws IOException {
        for (int n : new int[] { 10, 20 }) {
            Path edgeList = tempDir.resolve(
                    "rev-higgs-social/N=" + n + "/randomization=edge-swap/0.csv");
            Files.createDirectories(edgeList.getParent());
            Files.createFile(edgeList);
        }

        assertThrows(IOException.class, () -> SAR.resolveEdgeListPath(
                tempDir, Path.of("rev-higgs-social/N=500000"), true, 0));
    }
}
