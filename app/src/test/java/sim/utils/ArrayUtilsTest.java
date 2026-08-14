package sim.utils;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

class ArrayUtilsTest {

    @Test
    void generatesGeometricallySpacedValuesIncludingEndpoints() {
        assertArrayEquals(
                new double[] { 0.0001, 0.001, 0.01, 0.1 },
                ArrayUtils.geomspace(0.0001, 0.1, 4),
                1.0e-15);
    }

    @Test
    void generatesDescendingGeometricallySpacedValues() {
        assertArrayEquals(
                new double[] { 100.0, 10.0, 1.0 },
                ArrayUtils.geomspace(100.0, 1.0, 3),
                1.0e-12);
    }

    @Test
    void oneElementGeomspaceContainsStart() {
        assertArrayEquals(new double[] { 0.01 }, ArrayUtils.geomspace(0.01, 1.0, 1));
    }

    @Test
    void geomspaceRejectsNonPositiveOrNonFiniteEndpointsAndInvalidCount() {
        assertThrows(IllegalArgumentException.class, () -> ArrayUtils.geomspace(0.0, 1.0, 3));
        assertThrows(IllegalArgumentException.class, () -> ArrayUtils.geomspace(1.0, -1.0, 3));
        assertThrows(IllegalArgumentException.class,
                () -> ArrayUtils.geomspace(Double.NaN, 1.0, 3));
        assertThrows(IllegalArgumentException.class,
                () -> ArrayUtils.geomspace(1.0, Double.POSITIVE_INFINITY, 3));
        assertThrows(IllegalArgumentException.class, () -> ArrayUtils.geomspace(1.0, 10.0, 0));
    }
}
