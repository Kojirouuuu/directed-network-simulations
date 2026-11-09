package sim.simulation;

import java.util.*;
import sim.network.DirectedGraph;

/** Continuous-time SPIR, with an exact SIR mode that bypasses P. */
public final class SPIRSimulator {
    public enum ModelMode { SPIR, SIR }
    public enum Engine { GILLESPIE, EVENT_DRIVEN }
    public enum RecordingMode { FINAL_ONLY, FIXED_INTERVAL, EVENT_TIMES }
    private Engine engine = Engine.GILLESPIE;
    private RecordingMode recording = RecordingMode.EVENT_TIMES;
    private int[] watchingObservers = new int[16];
    private int watchingCount;
    private final EventHeap queue = new EventHeap();
    private long sequence;

    /** Immutable graph indices reusable across runs and threads. */
    public static final class PreparedGraph {
        final DirectedGraph graph;
        final int[] incomingPtr, incomingEdges;
        public PreparedGraph(DirectedGraph graph) {
            if (graph == null) throw new IllegalArgumentException("graph is required");
            this.graph = graph;
            incomingPtr = new int[graph.n + 1];
            incomingEdges = new int[graph.m];
            for (int e = 0; e < graph.m; e++) incomingPtr[graph.getOutNeighbor(e) + 1]++;
            for (int n = 1; n <= graph.n; n++) incomingPtr[n] += incomingPtr[n - 1];
            int[] cursor = incomingPtr.clone();
            for (int e = 0; e < graph.m; e++) incomingEdges[cursor[graph.getOutNeighbor(e)]++] = e;
        }
    }

    private enum Status { S, P, I, R }
    private final DirectedGraph graph;
    private final ModelMode mode;
    private final double lambda1, lambda2, mu1, mu2, end;
    private final Random random;
    private final FixedIntervalSampler sampler;
    private final Status[] state;
    private final int[] incomingPtr, incomingEdges;
    private final Pool presymptomaticNodes, infectiousNodes, presymptomaticEdges, infectiousEdges;
    private final SPIRResult result;
    private final boolean[] pendingObservers;
    private boolean initializing;
    private static final class Pool {
        final int[] values, positions;
        int size;
        Pool(int capacity) {
            values = new int[capacity];
            positions = new int[capacity];
            Arrays.fill(positions, -1);
        }
        void add(int value) {
            if (positions[value] >= 0)
                return;
            positions[value] = size;
            values[size++] = value;
        }
        void remove(int value) {
            if (value >= positions.length)
                return;
            int pos = positions[value];
            if (pos < 0)
                return;
            int last = values[--size];
            values[pos] = last;
            positions[last] = pos;
            positions[value] = -1;
        }
        int pick(Random rng) {
            return values[rng.nextInt(size)];
        }
    }

    private SPIRSimulator(DirectedGraph graph, ModelMode mode, double lambda1, double lambda2, double mu1,
            double mu2, double end, long seed, FixedIntervalSampler sampler, int[] observers,
            PreparedGraph prepared, Engine engine) {
        if (graph == null || mode == null)
            throw new IllegalArgumentException("graph and mode are required");
        for (double rate : new double[] {lambda1, lambda2, mu1, mu2})
            if (!Double.isFinite(rate) || rate < 0)
                throw new IllegalArgumentException("rates must be finite and non-negative");
        if (!Double.isFinite(end) || end <= 0)
            throw new IllegalArgumentException("tMax must be finite and positive");
        this.graph = graph;
        this.engine = engine;
        this.mode = mode;
        this.lambda1 = lambda1;
        this.lambda2 = lambda2;
        this.mu1 = mu1;
        this.mu2 = mu2;
        this.end = end;
        this.random = new Random(seed);
        this.sampler = sampler;
        state = new Status[graph.n];
        Arrays.fill(state, Status.S);
        incomingPtr = prepared.incomingPtr;
        incomingEdges = prepared.incomingEdges;
        presymptomaticNodes = new Pool(graph.n);
        infectiousNodes = new Pool(graph.n);
        presymptomaticEdges = new Pool(engine == Engine.GILLESPIE ? graph.m : 0);
        infectiousEdges = new Pool(engine == Engine.GILLESPIE ? graph.m : 0);
        result = new SPIRResult(graph.n);
        pendingObservers = new boolean[graph.n];
        if (observers == null)
            throw new IllegalArgumentException("observer nodes are required");
        for (int node : observers) {
            if (node < 0 || node >= graph.n)
                throw new IllegalArgumentException("invalid observer: " + node);
            if (!result.observers[node]) {
                result.observers[node] = true;
                result.observerCount++;
            }
        }
    }

    private void activate(int node) {
        if (engine == Engine.EVENT_DRIVEN) {
            schedule(node);
            return;
        }
        if (result.quarantined[node])
            return;
        var range = graph.outNeighborRange(node);
        for (int edge = range.start; edge < range.end; edge++)
            if (state[graph.getOutNeighbor(edge)] == Status.S
                    && !result.quarantined[graph.getOutNeighbor(edge)])
                (state[node] == Status.P ? presymptomaticEdges : infectiousEdges).add(edge);
    }

    private void quarantine(int node, double time) {
        if (result.quarantined[node])
            return;
        result.quarantined[node] = true;
        result.tQuarantine[node] = time;
        result.quarantinedCount++;
        var outgoing = graph.outNeighborRange(node);
        for (int edge = outgoing.start; edge < outgoing.end; edge++) {
            presymptomaticEdges.remove(edge);
            infectiousEdges.remove(edge);
        }
        for (int k = incomingPtr[node]; k < incomingPtr[node + 1]; k++) {
            int edge = incomingEdges[k];
            presymptomaticEdges.remove(edge);
            infectiousEdges.remove(edge);
        }
    }

    /** Monitor adjacency ignores arc direction; transmission still follows the original arcs. */
    private void isolateNeighborhood(int observer, double time) {
        quarantine(observer, time);
        var outgoing = graph.outNeighborRange(observer);
        for (int edge = outgoing.start; edge < outgoing.end; edge++) {
            int neighbor = graph.getOutNeighbor(edge);
            if (state[neighbor] != Status.I)
                quarantine(neighbor, time);
        }
        var incoming = graph.inNeighborRange(observer);
        for (int edge = incoming.start; edge < incoming.end; edge++) {
            int neighbor = graph.getInNeighbor(edge);
            if (state[neighbor] != Status.I)
                quarantine(neighbor, time);
        }
    }

    private void addWatchingObserver(int node) {
        if (result.observers[node] && !result.quarantined[node] && !pendingObservers[node]) {
            pendingObservers[node] = true;
            if (watchingCount == watchingObservers.length)
                watchingObservers = Arrays.copyOf(watchingObservers, watchingCount * 2);
            watchingObservers[watchingCount++] = node;
        }
    }

    private void detectOnset(int node, double time) {
        if (result.observerCount == 0 || result.quarantined[node])
            return;
        // All observers active at onset react together, independently of processing order.
        watchingCount = 0;
        addWatchingObserver(node);
        var outgoing = graph.outNeighborRange(node);
        for (int edge = outgoing.start; edge < outgoing.end; edge++)
            addWatchingObserver(graph.getOutNeighbor(edge));
        var incoming = graph.inNeighborRange(node);
        for (int edge = incoming.start; edge < incoming.end; edge++)
            addWatchingObserver(graph.getInNeighbor(edge));
        for (int index = 0; index < watchingCount; index++) {
            int observer = watchingObservers[index];
            isolateNeighborhood(observer, time);
            pendingObservers[observer] = false;
        }
    }

    private void infect(int node, double time) {
        state[node] = mode == ModelMode.SPIR ? Status.P : Status.I;
        result.tInfect[node] = time;
        if (state[node] == Status.P)
            presymptomaticNodes.add(node);
        else {
            infectiousNodes.add(node);
            result.tOnset[node] = time;
        }
        for (int k = incomingPtr[node]; k < incomingPtr[node + 1]; k++) {
            int edge = incomingEdges[k];
            presymptomaticEdges.remove(edge);
            infectiousEdges.remove(edge);
        }
        activate(node);
        if (state[node] == Status.I && !initializing)
            detectOnset(node, time);
    }

    private void recordState(double time) {
        if (recording == RecordingMode.FINAL_ONLY)
            return;
        result.times.add(time);
        result.S.add(graph.n - presymptomaticNodes.size - infectiousNodes.size - result.recovered);
        result.P.add(presymptomaticNodes.size);
        result.I.add(infectiousNodes.size);
        result.R.add(result.recovered);
        result.Q.add(result.quarantinedCount);
    }

    private SPIRResult run(int[] initial) {
        if (initial == null)
            throw new IllegalArgumentException("initial nodes are required");
        for (int node : initial) {
            if (node < 0 || node >= graph.n)
                throw new IllegalArgumentException("invalid initial node: " + node);
            if (result.observers[node])
                throw new IllegalArgumentException("initial infection is an observer: " + node);
        }
        initializing = true;
        for (int node : initial) {
            if (state[node] == Status.S)
                infect(node, 0);
        }
        initializing = false;
        if (mode == ModelMode.SIR)
            for (int node : initial) detectOnset(node, 0);
        recordState(0);
        double time = 0;
        while (engine == Engine.GILLESPIE) {
            double presymptomaticTransmissionRate = presymptomaticEdges.size * lambda1,
                   infectiousTransmissionRate = infectiousEdges.size * lambda2,
                   onsetRate = presymptomaticNodes.size * mu1, recoveryRate = infectiousNodes.size * mu2;
            double total =
                    presymptomaticTransmissionRate + infectiousTransmissionRate + onsetRate + recoveryRate;
            if (!Double.isFinite(total))
                throw new IllegalArgumentException("total rate overflow");
            if (total == 0)
                break;
            double next = time - Math.log1p(-random.nextDouble()) / total;
            if (next >= end)
                break;
            sampler.recordBefore(next, this::recordState);
            time = next;
            double choice = random.nextDouble() * total;
            if (choice < presymptomaticTransmissionRate)
                infect(graph.getOutNeighbor(presymptomaticEdges.pick(random)), time);
            else if (choice < presymptomaticTransmissionRate + infectiousTransmissionRate)
                infect(graph.getOutNeighbor(infectiousEdges.pick(random)), time);
            else if (choice < presymptomaticTransmissionRate + infectiousTransmissionRate + onsetRate) {
                onset(presymptomaticNodes.pick(random), time);
            } else {
                recover(infectiousNodes.pick(random), time);
            }
            if (sampler.isEnabled())
                sampler.recordThrough(time, this::recordState);
            else
                recordState(time);
        }
        if (engine == Engine.EVENT_DRIVEN)
            runQueue();
        sampler.recordRemaining(this::recordState);
        result.finalS = graph.n - presymptomaticNodes.size - infectiousNodes.size - result.recovered;
        result.finalP = presymptomaticNodes.size;
        result.finalI = infectiousNodes.size;
        result.finalTime = end;
        if (!sampler.isEnabled())
            recordState(end);
        if (recording == RecordingMode.FINAL_ONLY) {
            recording = RecordingMode.EVENT_TIMES;
            recordState(end);
        }
        return result;
    }

    public static SPIRResult simulate(PreparedGraph prepared, ModelMode mode, Engine engine,
            RecordingMode recording, double lambda1, double lambda2, double mu1, double mu2, double tMax,
            double dt, int[] initial, int[] observers, long seed) {
        Objects.requireNonNull(engine);
        Objects.requireNonNull(recording);
        var simulator = new SPIRSimulator(prepared.graph, mode, lambda1, lambda2, mu1, mu2, tMax, seed,
                recording == RecordingMode.FIXED_INTERVAL ? FixedIntervalSampler.create(dt, tMax)
                                                          : FixedIntervalSampler.disabled(),
                observers, prepared, engine);
        simulator.engine = engine;
        simulator.recording = recording;
        return simulator.run(initial);
    }

    private void onset(int node, double time) {
        presymptomaticNodes.remove(node);
        infectiousNodes.add(node);
        state[node] = Status.I;
        result.tOnset[node] = time;
        var range = graph.outNeighborRange(node);
        for (int e = range.start; e < range.end; e++) presymptomaticEdges.remove(e);
        activate(node);
        detectOnset(node, time);
    }

    private void recover(int node, double time) {
        infectiousNodes.remove(node);
        state[node] = Status.R;
        result.recovered++;
        result.tRecover[node] = time;
        var range = graph.outNeighborRange(node);
        for (int e = range.start; e < range.end; e++) infectiousEdges.remove(e);
    }

    private double exponential(double rate) {
        return rate == 0 ? Double.POSITIVE_INFINITY : -Math.log1p(-random.nextDouble()) / rate;
    }

    private double clock;
    private void schedule(int node) {
        boolean pre = state[node] == Status.P;
        double transition = clock + exponential(pre ? mu1 : mu2);
        if (transition < end)
            queue.add(transition, pre ? 0 : 1, node, node, sequence++);
        if (result.quarantined[node])
            return;
        var range = graph.outNeighborRange(node);
        for (int e = range.start; e < range.end; e++) {
            int target = graph.getOutNeighbor(e);
            if (state[target] != Status.S || result.quarantined[target])
                continue;
            double time = clock + exponential(pre ? lambda1 : lambda2);
            if (time < transition && time < end)
                queue.add(time, pre ? 2 : 3, node, target, sequence++);
        }
    }

    private void runQueue() {
        while (queue.size > 0) {
            double time = queue.times[0];
            int kind = queue.kinds[0], source = queue.sources[0], target = queue.targets[0];
            queue.remove();
            sampler.recordBefore(time, this::recordState);
            clock = time;
            if (kind == 0 && state[source] == Status.P)
                onset(source, time);
            else if (kind == 1 && state[source] == Status.I)
                recover(source, time);
            else if (kind >= 2 && !result.quarantined[source] && !result.quarantined[target]
                    && state[target] == Status.S && state[source] == (kind == 2 ? Status.P : Status.I))
                infect(target, time);
            if (!sampler.isEnabled())
                recordState(time);
            if (queue.size == 0 || queue.times[0] != time)
                sampler.recordThrough(time, this::recordState);
        }
    }

    /** Primitive heap: stable ordering, health transitions before transmission at equal times. */
    private static final class EventHeap {
        int size;
        double[] times = new double[16];
        int[] kinds = new int[16], sources = new int[16], targets = new int[16];
        long[] sequences = new long[16];
        boolean less(int a, int b) {
            int c = Double.compare(times[a], times[b]);
            if (c != 0)
                return c < 0;
            c = Integer.compare(kinds[a], kinds[b]);
            return c == 0 ? sequences[a] < sequences[b] : c < 0;
        }
        void swap(int a, int b) {
            double t = times[a];
            times[a] = times[b];
            times[b] = t;
            int k = kinds[a];
            kinds[a] = kinds[b];
            kinds[b] = k;
            k = sources[a];
            sources[a] = sources[b];
            sources[b] = k;
            k = targets[a];
            targets[a] = targets[b];
            targets[b] = k;
            long q = sequences[a];
            sequences[a] = sequences[b];
            sequences[b] = q;
        }
        void add(double time, int kind, int source, int target, long sequence) {
            if (size == times.length) {
                int capacity = Math.multiplyExact(size, 2);
                times = Arrays.copyOf(times, capacity);
                kinds = Arrays.copyOf(kinds, capacity);
                sources = Arrays.copyOf(sources, capacity);
                targets = Arrays.copyOf(targets, capacity);
                sequences = Arrays.copyOf(sequences, capacity);
            }
            int pos = size++;
            times[pos] = time;
            kinds[pos] = kind;
            sources[pos] = source;
            targets[pos] = target;
            sequences[pos] = sequence;
            while (pos > 0) {
                int parent = (pos - 1) / 2;
                if (!less(pos, parent))
                    break;
                swap(pos, parent);
                pos = parent;
            }
        }
        void remove() {
            swap(0, --size);
            int pos = 0;
            while (2 * pos + 1 < size) {
                int child = 2 * pos + 1;
                if (child + 1 < size && less(child + 1, child))
                    child++;
                if (!less(child, pos))
                    break;
                swap(pos, child);
                pos = child;
            }
        }
    }
}
