package sim.network.topology;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.Random;

import org.junit.jupiter.api.Test;

class PowPowSwapTest {
    @Test
    void optimizedSwapMatchesOriginalStubSampling() throws Exception {
        int[] ko = { 2, 5, 1, 8, 4, 7, 3 };
        int[] ki = { 7, 2, 8, 1, 5, 3, 4 };
        Method method = PowPow.class.getDeclaredMethod(
                "correlateBySwap", int[].class, int[].class, int.class, long.class);
        method.setAccessible(true);

        for (int swapNum : new int[] { -500, -30, 30, 500 }) {
            int[] actual = ki.clone();
            int[] expected = ki.clone();
            originalSwap(ko, expected, swapNum, 1234L);
            method.invoke(null, ko, actual, swapNum, 1234L);
            assertArrayEquals(expected, actual, "swapNum=" + swapNum);
        }
        InvocationTargetException error = assertThrows(InvocationTargetException.class,
                () -> method.invoke(null, ko, ki.clone(), Integer.MIN_VALUE, 1234L));
        assertInstanceOf(IllegalArgumentException.class, error.getCause());
    }

    private static void originalSwap(int[] ko, int[] ki, int swapNum, long seed) {
        int n = ko.length;
        int m = Arrays.stream(ko).sum();
        Random rng = new Random(seed);
        for (int trial = 0; trial < Math.abs(swapNum); trial++) {
            int[] outStubs = new int[m];
            int[] inStubs = new int[m];
            int pOut = 0;
            int pIn = 0;
            for (int u = 0; u < n; u++) {
                for (int t = 0; t < ko[u]; t++) outStubs[pOut++] = u;
                for (int t = 0; t < ki[u]; t++) inStubs[pIn++] = u;
            }
            if (pOut != m || pIn != m) throw new IllegalStateException();
            int i = outStubs[rng.nextInt(m)];
            int j = outStubs[rng.nextInt(m)];
            if (i == j) continue;
            boolean koGreater = ko[i] > ko[j];
            boolean kiGreater = ki[i] > ki[j];
            boolean swap = swapNum > 0
                    ? koGreater != kiGreater
                    : koGreater == kiGreater && ki[i] != ki[j];
            if (swap) {
                int value = ki[i];
                ki[i] = ki[j];
                ki[j] = value;
            }
        }
    }
}
