package sim;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

class BondPercolationOutputPathTest {
    @TempDir Path temp;

    @Test
    void outputPairSkipsSuffixOccupiedByEitherFile() throws IOException {
        var first = BondPercolation.resolveOutputPaths(temp, 0);
        assertEquals(temp.resolve("results_00.csv"), first.results());
        assertEquals(temp.resolve("component_sizes_00.csv"), first.componentSizes());
        Files.createFile(first.componentSizes());
        var second = BondPercolation.resolveOutputPaths(temp, 0);
        assertEquals(temp.resolve("results_00 (1).csv"), second.results());
        assertEquals(temp.resolve("component_sizes_00 (1).csv"), second.componentSizes());
        Files.createFile(second.results());
        var third = BondPercolation.resolveOutputPaths(temp, 0);
        assertEquals(temp.resolve("results_00 (2).csv"), third.results());
        assertEquals(temp.resolve("component_sizes_00 (2).csv"), third.componentSizes());
    }
}
