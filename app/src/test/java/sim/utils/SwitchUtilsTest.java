package sim.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.nio.file.Path;

import org.junit.jupiter.api.Test;

class SwitchUtilsTest {

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
    void degreeDistributionSarPathRequiresItsParameters() {
        assertThrows(IllegalArgumentException.class, () -> SwitchUtils.buildNetworkPath(
                "DegreeDistributionSAR", 100,
                null, null, null, null, null, null,
                null, 10, null, null, null, null,
                2.5, null, null, null, null));
    }
}
