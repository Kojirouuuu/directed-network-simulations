package sim.network;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import static sim.network.GraphOperations.edgeKey;
import static sim.network.GraphOperations.edgeCountAfterRemoval;

/** Mutable edge-index view used while searching for FFL-increasing swaps. */
final class FflRewiringState {
    final int vertexCount;
    final int edgeCount;
    final int[] outEdgePtr;
    final int[] inEdgePtr;
    final int[] sources;
    final int[] destinations;
    final int[] inEdgeIndices;
    final int[] inPositionByEdge;
    final boolean[] representativeEdge;
    final Map<Long, Integer> edgeCounts;
    final long[] cumulativeWedgeWeights;
    final long totalWedgeWeight;
    long feedForwardLoopCount;

    FflRewiringState(DirectedGraph graph) {
        vertexCount = graph.n;
        edgeCount = graph.m;
        var view = graph.csrView();
        outEdgePtr = view.outPtr();
        inEdgePtr = view.inPtr();
        sources = new int[edgeCount];
        destinations = Arrays.copyOf(view.outIdx(), edgeCount);
        inEdgeIndices = new int[edgeCount];
        inPositionByEdge = new int[edgeCount];
        representativeEdge = new boolean[edgeCount];
        edgeCounts = new HashMap<>();

        int[] nextInPosition = Arrays.copyOf(inEdgePtr, vertexCount);
        for (int source = 0; source < vertexCount; source++) {
            for (int edge = outEdgePtr[source]; edge < outEdgePtr[source + 1]; edge++) {
                sources[edge] = source;
                int destination = destinations[edge];
                int position = nextInPosition[destination]++;
                inEdgeIndices[position] = edge;
                inPositionByEdge[edge] = position;

                long key = edgeKey(source, destination);
                int previousCount = edgeCounts.getOrDefault(key, 0);
                edgeCounts.put(key, previousCount + 1);
                if (previousCount == 0) {
                    representativeEdge[edge] = true;
                }
            }
        }

        cumulativeWedgeWeights = new long[vertexCount];
        long cumulative = 0L;
        for (int vertex = 0; vertex < vertexCount; vertex++) {
            long inDegree = inEdgePtr[vertex + 1] - inEdgePtr[vertex];
            long outDegree = outEdgePtr[vertex + 1] - outEdgePtr[vertex];
            cumulative = Math.addExact(cumulative, Math.multiplyExact(inDegree, outDegree));
            cumulativeWedgeWeights[vertex] = cumulative;
        }
        totalWedgeWeight = cumulative;
    }

    long countAllFeedForwardLoops() {
        long count = 0L;
        for (int edge = 0; edge < edgeCount; edge++) {
            if (!representativeEdge[edge]) {
                continue;
            }
            int source = sources[edge];
            int destination = destinations[edge];
            if (source != destination) {
                count = Math.addExact(count, countOutgoingIncomingIntersection(source, destination));
            }
        }
        return count;
    }

    boolean tryIncreasingSwap(Random random) {
        if (totalWedgeWeight == 0L) {
            return false;
        }

        int middle = sampleMiddleVertex(random);
        int inDegree = inEdgePtr[middle + 1] - inEdgePtr[middle];
        int outDegree = outEdgePtr[middle + 1] - outEdgePtr[middle];
        int firstPathEdge = inEdgeIndices[inEdgePtr[middle] + random.nextInt(inDegree)];
        int secondPathEdge = outEdgePtr[middle] + random.nextInt(outDegree);
        int source = sources[firstPathEdge];
        int destination = destinations[secondPathEdge];

        if (source == middle || middle == destination || source == destination
                || edgeExists(source, destination)) {
            return false;
        }

        int firstSwapEdge = outEdgePtr[source]
                + random.nextInt(outEdgePtr[source + 1] - outEdgePtr[source]);
        int secondSwapEdge = inEdgeIndices[inEdgePtr[destination]
                + random.nextInt(inEdgePtr[destination + 1] - inEdgePtr[destination])];
        int firstOldDestination = destinations[firstSwapEdge];
        int secondSource = sources[secondSwapEdge];

        // The two path edges must remain in place so that a→b→c is still present after the swap.
        if (firstOldDestination == middle || secondSource == middle
                || firstSwapEdge == secondSwapEdge || secondSource == firstOldDestination) {
            return false;
        }

        long oldFirst = edgeKey(source, firstOldDestination);
        long oldSecond = edgeKey(secondSource, destination);
        long newFirst = edgeKey(source, destination);
        long newSecond = edgeKey(secondSource, firstOldDestination);
        if (source == firstOldDestination || secondSource == destination
                || edgeCounts.getOrDefault(oldFirst, 0) != 1
                || edgeCounts.getOrDefault(oldSecond, 0) != 1
                || newFirst == newSecond
                || edgeCountAfterRemoval(edgeCounts, newFirst, oldFirst, oldSecond) > 0
                || edgeCountAfterRemoval(edgeCounts, newSecond, oldFirst, oldSecond) > 0) {
            return false;
        }

        long before = feedForwardLoopCount;
        swapDestinations(firstSwapEdge, secondSwapEdge);
        if (feedForwardLoopCount > before) {
            return true;
        }

        swapDestinations(firstSwapEdge, secondSwapEdge);
        if (feedForwardLoopCount != before) {
            throw new IllegalStateException("FFL count was not restored after rejecting a swap");
        }
        return false;
    }

    int sampleMiddleVertex(Random random) {
        long ticket = random.nextLong(totalWedgeWeight) + 1L;
        int low = 0;
        int high = cumulativeWedgeWeights.length - 1;
        while (low < high) {
            int middle = (low + high) >>> 1;
            if (cumulativeWedgeWeights[middle] >= ticket) {
                high = middle;
            } else {
                low = middle + 1;
            }
        }
        return low;
    }

    void swapDestinations(int firstEdge, int secondEdge) {
        int firstDestination = destinations[firstEdge];
        int secondDestination = destinations[secondEdge];

        removeEdge(firstEdge);
        removeEdge(secondEdge);

        destinations[firstEdge] = secondDestination;
        destinations[secondEdge] = firstDestination;

        int firstInPosition = inPositionByEdge[firstEdge];
        int secondInPosition = inPositionByEdge[secondEdge];
        inEdgeIndices[firstInPosition] = secondEdge;
        inEdgeIndices[secondInPosition] = firstEdge;
        inPositionByEdge[firstEdge] = secondInPosition;
        inPositionByEdge[secondEdge] = firstInPosition;

        addEdge(firstEdge);
        addEdge(secondEdge);
    }

    void removeEdge(int edge) {
        int source = sources[edge];
        int destination = destinations[edge];
        long key = edgeKey(source, destination);
        int count = edgeCounts.getOrDefault(key, 0);
        if (count <= 0) {
            throw new IllegalStateException("Cannot remove an edge that does not exist");
        }

        if (count == 1) {
            feedForwardLoopCount = Math.subtractExact(
                    feedForwardLoopCount, countFeedForwardLoopsContaining(source, destination));
            edgeCounts.remove(key);
            representativeEdge[edge] = false;
            return;
        }

        edgeCounts.put(key, count - 1);
        if (representativeEdge[edge]) {
            representativeEdge[edge] = false;
            int replacement = findParallelEdge(source, destination, edge);
            if (replacement < 0) {
                throw new IllegalStateException("Could not find a representative parallel edge");
            }
            representativeEdge[replacement] = true;
        }
    }

    void addEdge(int edge) {
        int source = sources[edge];
        int destination = destinations[edge];
        long key = edgeKey(source, destination);
        int count = edgeCounts.getOrDefault(key, 0);
        if (count == 0) {
            edgeCounts.put(key, 1);
            representativeEdge[edge] = true;
            feedForwardLoopCount = Math.addExact(
                    feedForwardLoopCount, countFeedForwardLoopsContaining(source, destination));
        } else {
            edgeCounts.put(key, count + 1);
            representativeEdge[edge] = false;
        }
    }

    int findParallelEdge(int source, int destination, int excludedEdge) {
        for (int edge = outEdgePtr[source]; edge < outEdgePtr[source + 1]; edge++) {
            if (edge != excludedEdge && destinations[edge] == destination) {
                return edge;
            }
        }
        return -1;
    }

    long countFeedForwardLoopsContaining(int source, int destination) {
        if (source == destination) {
            return 0L;
        }
        long count = countCommonOutgoing(source, destination);
        count = Math.addExact(count, countCommonIncoming(source, destination));
        return Math.addExact(count, countOutgoingIncomingIntersection(source, destination));
    }

    long countCommonOutgoing(int first, int second) {
        int firstDegree = outEdgePtr[first + 1] - outEdgePtr[first];
        int secondDegree = outEdgePtr[second + 1] - outEdgePtr[second];
        int iterated = firstDegree <= secondDegree ? first : second;
        int checked = iterated == first ? second : first;
        long count = 0L;
        for (int edge = outEdgePtr[iterated]; edge < outEdgePtr[iterated + 1]; edge++) {
            if (!representativeEdge[edge]) {
                continue;
            }
            int candidate = destinations[edge];
            if (candidate != first && candidate != second && edgeExists(checked, candidate)) {
                count++;
            }
        }
        return count;
    }

    long countCommonIncoming(int first, int second) {
        int firstDegree = inEdgePtr[first + 1] - inEdgePtr[first];
        int secondDegree = inEdgePtr[second + 1] - inEdgePtr[second];
        int iterated = firstDegree <= secondDegree ? first : second;
        int checked = iterated == first ? second : first;
        long count = 0L;
        for (int position = inEdgePtr[iterated]; position < inEdgePtr[iterated + 1]; position++) {
            int edge = inEdgeIndices[position];
            if (!representativeEdge[edge]) {
                continue;
            }
            int candidate = sources[edge];
            if (candidate != first && candidate != second && edgeExists(candidate, checked)) {
                count++;
            }
        }
        return count;
    }

    long countOutgoingIncomingIntersection(int source, int destination) {
        int outDegree = outEdgePtr[source + 1] - outEdgePtr[source];
        int inDegree = inEdgePtr[destination + 1] - inEdgePtr[destination];
        long count = 0L;
        if (outDegree <= inDegree) {
            for (int edge = outEdgePtr[source]; edge < outEdgePtr[source + 1]; edge++) {
                if (!representativeEdge[edge]) {
                    continue;
                }
                int middle = destinations[edge];
                if (middle != source && middle != destination && edgeExists(middle, destination)) {
                    count++;
                }
            }
        } else {
            for (int position = inEdgePtr[destination]; position < inEdgePtr[destination + 1]; position++) {
                int edge = inEdgeIndices[position];
                if (!representativeEdge[edge]) {
                    continue;
                }
                int middle = sources[edge];
                if (middle != source && middle != destination && edgeExists(source, middle)) {
                    count++;
                }
            }
        }
        return count;
    }

    boolean edgeExists(int source, int destination) {
        return edgeCounts.getOrDefault(edgeKey(source, destination), 0) > 0;
    }
}
