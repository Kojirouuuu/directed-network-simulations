package sim.simulation;

import org.junit.jupiter.api.Test;
import sim.utils.ArrayUtils;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class InitialAdoptersTest {
    @Test
    void samplingPreservesShuffleOrderAndMinimumOneAdopter() {
        int[] original = {0, 1, 2, 3, 4, 5, 6, 7};
        ArrayUtils.shuffle(original, 2007L);
        assertArrayEquals(Arrays.copyOfRange(original, 0, 4), InitialAdopters.sample(8, 0.5, 2007L));
        assertArrayEquals(Arrays.copyOfRange(original, 0, 1), InitialAdopters.sample(8, 0, 2007L));
    }

}
