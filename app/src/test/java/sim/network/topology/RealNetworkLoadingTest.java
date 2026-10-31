package sim.network.topology;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import sim.network.DirectedGraph;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.assertEquals;

class RealNetworkLoadingTest {
    @TempDir Path directory;

    @Test
    void realNetworkWrappersRetainEncounterOrderAndParsingRules() throws Exception {
        Path path = directory.resolve("edges.txt");
        Files.writeString(path, "# comment\n\nnode-z node-a\nignored\nnode-a node-q extra\nnode-z node-a\n");
        for (DirectedGraph graph : new DirectedGraph[]{
                DirectedGraph.loadFromEdgeList(null, path), EgoTwitter.loadFromEdgeList(null, path),
                Gplus.loadFromEdgeList(null, path), HiggsSocial.loadFromEdgeList(null, path)}) {
            assertEquals("edges.txt", graph.name);
            assertEquals(3, graph.n);
            assertEquals(3, graph.m);
            assertEquals(1, graph.getOutNeighbor(0));
            assertEquals(1, graph.getOutNeighbor(1));
            assertEquals(2, graph.getOutNeighbor(2));
        }
    }
}
