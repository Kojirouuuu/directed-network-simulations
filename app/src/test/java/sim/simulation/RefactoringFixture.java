package sim.simulation;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.HexFormat;
import sim.network.DirectedGraph;

/** Snapshot of traversal order and trajectories, captured from the pre-refactoring implementation. */
final class RefactoringFixture {
    static String fingerprint() throws Exception {
        DirectedGraph graph = DirectedGraph.fromEdgeListWithUndirectedFlag("fixture", 8,
                new int[]{0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7},
                new int[]{1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,1}, new boolean[16]);
        StringBuilder snapshot = new StringBuilder();
        for (DirectedGraph variant : new DirectedGraph[]{graph, graph.reverseDirectedEdges(),
                graph.randomizeByEdgeSwaps(42), graph.increaseReciprocity(0, 0, 20, 42),
                graph.randomizeByShuffledDegreeSequence(DirectedGraph.DegreeSide.IN, 42),
                graph.randomizeByShuffledDegreeSequence(DirectedGraph.DegreeSide.OUT, 42),
                graph.expandByRepeatedJointDegreeSequence(2, 42)}) {
            snapshot.append(variant.name).append(':').append(variant.n).append(':').append(variant.m);
            for (int u = 0; u < variant.n; u++) {
                var range = variant.outNeighborRange(u);
                for (int i = range.start; i < range.end; i++) {
                    snapshot.append('/').append(u).append('>').append(variant.getOutNeighbor(i))
                            .append(':').append(variant.isOutUndirected(i));
                }
            }
            snapshot.append('/').append(variant.countReciprocalArcs()).append('/')
                    .append(variant.countFeedForwardLoops());
        }
        for (boolean gillespie : new boolean[]{false, true}) {
            for (boolean fixed : new boolean[]{false, true}) {
                int[] thresholds = new int[8];
                Arrays.fill(thresholds, 2);
                SARResult result = gillespie
                        ? (fixed ? SARGillespieSimulator.simulate(graph, 2, 0, 1, 5, 0.25, thresholds, new int[]{0,1}, 12345)
                                 : SARGillespieSimulator.simulate(graph, 2, 0, 1, 5, thresholds, new int[]{0,1}, 12345))
                        : (fixed ? SARSimulator.simulate(graph, 2, 0, 1, 5, 0.25, thresholds, new int[]{0,1}, 12345)
                                 : SARSimulator.simulate(graph, 2, 0, 1, 5, thresholds, new int[]{0,1}, 12345));
                snapshot.append(result.times).append(result.S).append(result.A).append(result.R).append(result.Phi)
                        .append(result.initialAdoptedTime).append(result.finalAdoptedTime)
                        .append(Arrays.toString(result.tInfect)).append(Arrays.toString(result.tRecover));
            }
        }
        return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256")
                .digest(snapshot.toString().getBytes(StandardCharsets.UTF_8)));
    }

    public static void main(String[] args) throws Exception {
        System.out.println(fingerprint());
    }
}
