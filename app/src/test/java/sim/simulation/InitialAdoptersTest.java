package sim.simulation;

import org.junit.jupiter.api.Test;
import sim.utils.ArrayUtils;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class InitialAdoptersTest {
    @Test
    void samplingPreservesShuffleOrderAndMinimumOneAdopter() {
        int[] original = {0, 1, 2, 3, 4, 5, 6, 7};
        ArrayUtils.shuffle(original, 2007L);
        assertArrayEquals(Arrays.copyOfRange(original, 0, 4), InitialAdopters.sample(8, 0.5, 2007L));
        assertArrayEquals(Arrays.copyOfRange(original, 0, 1), InitialAdopters.sample(8, 0, 2007L));
    }

    @Test
    void excludesObserversWithoutChangingRequestedInitialCount() {
        boolean[] observers = {true, false, true, false, false, false, false, false};
        int[] sample = InitialAdopters.sample(8, 0.5, observers, 2007L);
        assertEquals(4, sample.length);
        assertEquals(4, Arrays.stream(sample).distinct().count());
        for (int node : sample) assertFalse(observers[node]);
        assertArrayEquals(sample, InitialAdopters.sample(8, 0.5, observers, 2007L));
        assertArrayEquals(InitialAdopters.sample(8, 0.5, 2007L),
                InitialAdopters.sample(8, 0.5, new boolean[8], 2007L));
    }

    @Test
    void rejectsInsufficientNonObserversAndInvalidInputs() {
        assertThrows(IllegalArgumentException.class,
                () -> InitialAdopters.sample(2, 0, new boolean[] {true, true}, 1));
        assertThrows(IllegalArgumentException.class,
                () -> InitialAdopters.sample(2, 1, new boolean[] {true, false}, 1));
        assertThrows(IllegalArgumentException.class,
                () -> InitialAdopters.sample(2, Double.NaN, new boolean[2], 1));
    }

}
