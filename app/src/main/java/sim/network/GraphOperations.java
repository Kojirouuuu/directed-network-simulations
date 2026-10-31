package sim.network;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import sim.network.DirectedGraph.DegreeSide;
import static sim.network.DirectedGraph.fromEdgeListWithUndirectedFlag;

/** Internal graph metrics and transformations; preserves CSR traversal and random draw order. */
final class GraphOperations {
    private final DirectedGraph graph;
    private final String name;
    private final int n, m;
    private final int[] outPtr, outIdx, inPtr;
    private final boolean[] outIsUndirected;

    GraphOperations(DirectedGraph graph) {
        this.graph = graph;
        name = graph.name;
        n = graph.n;
        m = graph.m;
        var view = graph.csrView();
        outPtr = view.outPtr();
        outIdx = view.outIdx();
        inPtr = view.inPtr();
        outIsUndirected = view.outIsUndirected();
    }

    /**
     * 逆向きの辺と対応する有向辺（相互アーク）の本数を返す。
     * 多重辺については、頂点対ごとに両方向の本数の小さい方を対応付ける。
     * 例えば u→v が3本、v→u が1本なら相互アークは2本と数える。
     * 自己ループは相互アークに含めない。
     *
     * @return 相互アーク数
     */
    public int countReciprocalArcs() {
        return countReciprocalArcs(buildEdgeCounts());
    }

    /**
     * 全アークに占める相互アークの割合を返す。
     *
     * @return 相互アーク数 / m。空グラフでは 0.0
     */
    public double reciprocity() {
        return m == 0 ? 0.0 : countReciprocalArcs() / (double) m;
    }

    /**
     * フィードフォワードループ（FFL）の個数を返す。
     *
     * <p>相異なる順序付き3頂点 (a, b, c) に対し、a→b, b→c, a→c がすべて
     * 存在するとき1個と数える。追加の逆辺などは許容する。自己ループは無視し、
     * 平行アークは存在辺1本として扱う。</p>
     *
     * @return FFL 数
     */
    public long countFeedForwardLoops() {
        long cached = graph.cachedFeedForwardLoopCount;
        if (cached >= 0L) {
            return cached;
        }

        FflRewiringState state = new FflRewiringState(graph);
        long count = state.countAllFeedForwardLoops();
        graph.cachedFeedForwardLoopCount = count;
        return count;
    }

    /**
     * 入次数・出次数を保存する辺スワップにより、FFL 数を目標値まで増加させる。
     * 元のグラフは変更しない。
     *
     * <p>a→b→c があり a→c がないとき、別の辺 a→x と y→c を選び、
     * a→c と y→x に交換する。新たな自己ループ・多重辺になる提案は棄却し、
     * FFL 数が厳密に増える交換だけを採択する。</p>
     *
     * @param targetNffl 目標 FFL 数
     * @param maxIncreaseAttempts FFL を増加させる候補の最大試行回数
     * @param seed 乱数シード
     * @return リワイヤリング後の新しい DirectedGraph
     * @throws IllegalArgumentException 引数が不正、または無向由来辺を含む場合
     * @throws IllegalStateException 最大試行回数内に目標へ到達できなかった場合
     */
    public DirectedGraph increaseFeedForwardLoops(long targetNffl, long maxIncreaseAttempts, long seed) {
        if (targetNffl < 0L) {
            throw new IllegalArgumentException("targetNffl must be non-negative");
        }
        if (maxIncreaseAttempts < 0L) {
            throw new IllegalArgumentException("maxIncreaseAttempts must be non-negative");
        }
        for (boolean isUndirected : outIsUndirected) {
            if (isUndirected) {
                throw new IllegalArgumentException("increaseFeedForwardLoops requires a purely directed graph");
            }
        }

        FflRewiringState state = new FflRewiringState(graph);
        long initialNffl = graph.cachedFeedForwardLoopCount;
        if (initialNffl < 0L) {
            initialNffl = state.countAllFeedForwardLoops();
            graph.cachedFeedForwardLoopCount = initialNffl;
        }
        state.feedForwardLoopCount = initialNffl;

        Random random = new Random(seed);
        long attempts = 0L;
        long acceptedSwaps = 0L;
        while (state.feedForwardLoopCount < targetNffl && attempts < maxIncreaseAttempts) {
            if (state.tryIncreasingSwap(random)) {
                acceptedSwaps++;
            }
            attempts++;
        }

        if (state.feedForwardLoopCount < targetNffl) {
            throw new IllegalStateException(
                    "Could not reach target NFFL " + targetNffl
                            + " in " + attempts + " attempts; accepted " + acceptedSwaps
                            + " increasing swaps; achieved " + state.feedForwardLoopCount);
        }

        DirectedGraph rewired = fromEdgeListWithUndirectedFlag(
                name + "_ffl", n, state.sources, state.destinations, new boolean[m]);
        rewired.cachedFeedForwardLoopCount = state.feedForwardLoopCount;
        return rewired;
    }

    /**
     * 入次数・出次数を保存する辺スワップにより、相互辺割合を目標値まで増加させる。
     * 目標到達後は、相互アーク数が変化しないスワップを繰り返して偏りを緩和する。
     * 元のグラフは変更しない。
     *
     * <p>辺 a→b と c→d に対し、a→d と c→b への交換を提案する。新たな自己ループ・
     * 多重辺になる提案は棄却する。既存の自己ループ・多重辺は入力として許容する。</p>
     *
     * @param targetReciprocity 目標相互辺割合（0以上1以下）
     * @param maxIncreaseAttempts 相互辺を増加させる候補の最大試行回数
     * @param neutralSwapAttempts 相互辺数を保存する候補の試行回数
     * @param seed 乱数シード
     * @return リワイヤリング後の新しい DirectedGraph
     * @throws IllegalArgumentException 引数が不正、または無向由来辺を含む場合
     * @throws IllegalStateException 最大試行回数内に目標へ到達できなかった場合
     */
    public DirectedGraph increaseReciprocity(double targetReciprocity, long maxIncreaseAttempts,
            long neutralSwapAttempts, long seed) {
        if (!Double.isFinite(targetReciprocity) || targetReciprocity < 0.0 || targetReciprocity > 1.0) {
            throw new IllegalArgumentException("targetReciprocity must be finite and in [0, 1]");
        }
        if (maxIncreaseAttempts < 0) {
            throw new IllegalArgumentException("maxIncreaseAttempts must be non-negative");
        }
        if (neutralSwapAttempts < 0) {
            throw new IllegalArgumentException("neutralSwapAttempts must be non-negative");
        }
        for (boolean isUndirected : outIsUndirected) {
            if (isUndirected) {
                throw new IllegalArgumentException("increaseReciprocity requires a purely directed graph");
            }
        }

        int[] srcs = new int[m];
        int[] dsts = new int[m];
        int edgeIndex = 0;
        for (int u = 0; u < n; u++) {
            for (int i = outPtr[u]; i < outPtr[u + 1]; i++) {
                srcs[edgeIndex] = u;
                dsts[edgeIndex] = outIdx[i];
                edgeIndex++;
            }
        }

        Map<Long, Integer> edgeCounts = buildEdgeCounts(srcs, dsts);
        int reciprocalArcs = countReciprocalArcs(edgeCounts);
        Random random = new Random(seed);

        long requiredReciprocalArcs = (long) Math.ceil(targetReciprocity * m);
        long missingReciprocalArcs = Math.max(0L, requiredReciprocalArcs - reciprocalArcs);
        long theoreticalMinimumAttempts = (missingReciprocalArcs + 3L) / 4L;
        if (maxIncreaseAttempts < theoreticalMinimumAttempts) {
            throw new IllegalStateException(
                    "At least " + theoreticalMinimumAttempts
                            + " accepted swaps are theoretically required to reach target reciprocity "
                            + targetReciprocity + ", but maxIncreaseAttempts is " + maxIncreaseAttempts);
        }

        long attempts = 0L;
        long acceptedIncreaseSwaps = 0L;
        while (!hasReachedTarget(reciprocalArcs, targetReciprocity) && attempts < maxIncreaseAttempts) {
            int updatedReciprocalArcs = tryDegreePreservingSwap(
                    srcs, dsts, edgeCounts, reciprocalArcs, random, true);
            if (updatedReciprocalArcs > reciprocalArcs) {
                acceptedIncreaseSwaps++;
            }
            reciprocalArcs = updatedReciprocalArcs;
            attempts++;
        }

        if (!hasReachedTarget(reciprocalArcs, targetReciprocity)) {
            double achieved = m == 0 ? 0.0 : reciprocalArcs / (double) m;
            throw new IllegalStateException(
                    "Could not reach target reciprocity " + targetReciprocity
                            + " in " + attempts + " attempts; accepted " + acceptedIncreaseSwaps
                            + " increasing swaps; achieved " + achieved);
        }

        for (long attempt = 0L; attempt < neutralSwapAttempts; attempt++) {
            reciprocalArcs = tryDegreePreservingSwap(
                    srcs, dsts, edgeCounts, reciprocalArcs, random, false);
        }

        return fromEdgeListWithUndirectedFlag(
                name + "_reciprocity", n, srcs, dsts, new boolean[m]);
    }

    /**
     * 指定した側の次数列を頂点間でシャッフルし、全ての有向辺をランダムに繋ぎ直す。
     * 元のグラフは変更しない。
     *
     * <p>{@code IN} の場合は各頂点の出次数を固定したまま入次数を再配置し、
     * {@code OUT} の場合は各頂点の入次数を固定したまま出次数を再配置する。
     * どちらの場合も入次数分布・出次数分布と辺数は保存される。</p>
     *
     * <p>再接続は有向 Configuration Model として行うため、自己ループと多重辺を
     * 許容する。</p>
     *
     * @param side 頂点間でシャッフルする次数列
     * @param seed 乱数シード
     * @return 次数列の再配置とランダム再接続を行った新しい DirectedGraph
     * @throws IllegalArgumentException side が null、または無向由来辺を含む場合
     */
    public DirectedGraph randomizeByShuffledDegreeSequence(DegreeSide side, long seed) {
        if (side == null) {
            throw new IllegalArgumentException("side must be non-null");
        }
        for (boolean isUndirected : outIsUndirected) {
            if (isUndirected) {
                throw new IllegalArgumentException(
                        "randomizeByShuffledDegreeSequence requires a purely directed graph");
            }
        }

        int[] inDegrees = new int[n];
        int[] outDegrees = new int[n];
        for (int vertex = 0; vertex < n; vertex++) {
            inDegrees[vertex] = inPtr[vertex + 1] - inPtr[vertex];
            outDegrees[vertex] = outPtr[vertex + 1] - outPtr[vertex];
        }

        Random random = new Random(seed);
        shuffle(side == DegreeSide.IN ? inDegrees : outDegrees, random);

        int[] outStubs = new int[m];
        int[] inStubs = new int[m];
        int outPosition = 0;
        int inPosition = 0;
        for (int vertex = 0; vertex < n; vertex++) {
            Arrays.fill(outStubs, outPosition, outPosition + outDegrees[vertex], vertex);
            outPosition += outDegrees[vertex];
            Arrays.fill(inStubs, inPosition, inPosition + inDegrees[vertex], vertex);
            inPosition += inDegrees[vertex];
        }

        if (outPosition != m || inPosition != m) {
            throw new IllegalStateException(
                    "Stub count mismatch: out=" + outPosition + ", in=" + inPosition + ", m=" + m);
        }

        shuffle(outStubs, random);
        shuffle(inStubs, random);

        String suffix = side == DegreeSide.IN
                ? "_in_degree_shuffled"
                : "_out_degree_shuffled";
        return fromEdgeListWithUndirectedFlag(
                name + suffix, n, outStubs, inStubs, new boolean[m]);
    }

    /**
     * 各頂点の入次数・出次数の組を指定倍率だけ複製し、有向 Configuration Model として
     * 全ての辺をランダムに繋ぎ直した拡張グラフを返す。元のグラフは変更しない。
     *
     * <p>頂点ごとの {@code (inDegree, outDegree)} を組のまま複製するため、拡張後の
     * 経験的同時次数分布は元グラフと厳密に一致する。自己ループと多重辺を許容する。</p>
     *
     * @param multiplier 頂点数と辺数の倍率（2以上）
     * @param seed スタブ再接続用の乱数シード
     * @return 同時次数分布を保存した拡張グラフ
     * @throws IllegalArgumentException 倍率が2未満、無向由来辺を含む、または拡張後の
     *         頂点数・辺数が {@code int} の範囲を超える場合
     */
    public DirectedGraph expandByRepeatedJointDegreeSequence(int multiplier, long seed) {
        if (multiplier < 2) {
            throw new IllegalArgumentException("multiplier must be at least 2");
        }
        for (boolean isUndirected : outIsUndirected) {
            if (isUndirected) {
                throw new IllegalArgumentException(
                        "expandByRepeatedJointDegreeSequence requires a purely directed graph");
            }
        }

        final int expandedN;
        final int expandedM;
        try {
            expandedN = Math.multiplyExact(n, multiplier);
            expandedM = Math.multiplyExact(m, multiplier);
        } catch (ArithmeticException e) {
            throw new IllegalArgumentException(
                    "expanded vertex or edge count exceeds the supported int range", e);
        }

        int[] sources = new int[expandedM];
        int[] destinations = new int[expandedM];
        int sourcePosition = 0;
        int destinationPosition = 0;
        for (int copy = 0; copy < multiplier; copy++) {
            int vertexOffset = copy * n;
            for (int vertex = 0; vertex < n; vertex++) {
                int expandedVertex = vertexOffset + vertex;
                int outDegree = outPtr[vertex + 1] - outPtr[vertex];
                int inDegree = inPtr[vertex + 1] - inPtr[vertex];
                Arrays.fill(sources, sourcePosition, sourcePosition + outDegree, expandedVertex);
                Arrays.fill(destinations, destinationPosition, destinationPosition + inDegree, expandedVertex);
                sourcePosition += outDegree;
                destinationPosition += inDegree;
            }
        }

        if (sourcePosition != expandedM || destinationPosition != expandedM) {
            throw new IllegalStateException(
                    "Stub count mismatch: out=" + sourcePosition + ", in=" + destinationPosition
                            + ", m=" + expandedM);
        }

        shuffle(destinations, new Random(seed));
        return fromEdgeListWithUndirectedFlag(
                name + "_joint_degree_cm_x" + multiplier,
                expandedN, sources, destinations, new boolean[expandedM]);
    }

    /**
     * 入次数・出次数を保存するランダムな辺スワップにより、グラフをランダム化する。
     * 元のグラフは変更しない。
     *
     * <p>辺 a→b と c→d に対し、a→d と c→b への交換を提案し、指標による
     * 採択条件を設けず、構造的に有効な交換をすべて採択する。新たな自己ループ・
     * 多重辺になる提案は棄却する。既存の自己ループ・多重辺は入力として許容する。</p>
     *
     * <p>成功した交換が辺数の10倍に達するまで試行する。辺数の100倍の試行でも
     * 目標に達しない場合は例外を送出する。</p>
     *
     * @param seed 乱数シード
     * @return ランダム化後の新しい DirectedGraph
     * @throws IllegalArgumentException 無向由来辺を含む場合
     * @throws IllegalStateException 試行上限内に必要な交換回数へ到達できなかった場合
     */
    public DirectedGraph randomizeByEdgeSwaps(long seed) {
        for (boolean isUndirected : outIsUndirected) {
            if (isUndirected) {
                throw new IllegalArgumentException("randomizeByEdgeSwaps requires a purely directed graph");
            }
        }

        int[] srcs = new int[m];
        int[] dsts = new int[m];
        int edgeIndex = 0;
        for (int source = 0; source < n; source++) {
            for (int edge = outPtr[source]; edge < outPtr[source + 1]; edge++) {
                srcs[edgeIndex] = source;
                dsts[edgeIndex] = outIdx[edge];
                edgeIndex++;
            }
        }

        Map<Long, Integer> edgeCounts = buildEdgeCounts(srcs, dsts);
        Random random = new Random(seed);
        long targetAcceptedSwaps = 10L * m;
        long maxAttempts = 100L * m;
        long attempts = 0L;
        long acceptedSwaps = 0L;

        while (acceptedSwaps < targetAcceptedSwaps && attempts < maxAttempts) {
            if (tryRandomDegreePreservingSwap(srcs, dsts, edgeCounts, random)) {
                acceptedSwaps++;
            }
            attempts++;
        }

        if (acceptedSwaps < targetAcceptedSwaps) {
            throw new IllegalStateException(
                    "Could not complete " + targetAcceptedSwaps
                            + " random edge swaps in " + attempts + " attempts; accepted "
                            + acceptedSwaps);
        }

        return fromEdgeListWithUndirectedFlag(
                name + "_randomized", n, srcs, dsts, new boolean[m]);
    }

    private boolean hasReachedTarget(int reciprocalArcs, double targetReciprocity) {
        if (m == 0) {
            return targetReciprocity == 0.0;
        }
        return reciprocalArcs / (double) m >= targetReciprocity;
    }

    static int tryDegreePreservingSwap(int[] srcs, int[] dsts, Map<Long, Integer> edgeCounts,
            int reciprocalArcs, Random random, boolean requireIncrease) {
        int edgeCount = srcs.length;
        if (edgeCount < 2) {
            return reciprocalArcs;
        }

        int first = random.nextInt(edgeCount);
        int second = random.nextInt(edgeCount - 1);
        if (second >= first) {
            second++;
        }

        int a = srcs[first];
        int b = dsts[first];
        int c = srcs[second];
        int d = dsts[second];

        if (a == d || c == b || a == c || b == d) {
            return reciprocalArcs;
        }

        long oldFirst = edgeKey(a, b);
        long oldSecond = edgeKey(c, d);
        long newFirst = edgeKey(a, d);
        long newSecond = edgeKey(c, b);
        if (newFirst == newSecond) {
            return reciprocalArcs;
        }

        if (edgeCountAfterRemoval(edgeCounts, newFirst, oldFirst, oldSecond) > 0
                || edgeCountAfterRemoval(edgeCounts, newSecond, oldFirst, oldSecond) > 0) {
            return reciprocalArcs;
        }

        // 新しい2辺のどちらにも交換後の逆辺がなければ、相互アーク数は増えない。
        // 増加フェーズでは大多数の候補をここで棄却し、以降の差分計算を避ける。
        if (requireIncrease
                && projectedEdgeCount(edgeCounts, edgeKey(d, a), oldFirst, oldSecond, newFirst, newSecond) == 0
                && projectedEdgeCount(edgeCounts, edgeKey(b, c), oldFirst, oldSecond, newFirst, newSecond) == 0) {
            return reciprocalArcs;
        }

        long firstOldDyad = dyadKey(a, b);
        long secondOldDyad = dyadKey(c, d);
        long firstNewDyad = dyadKey(a, d);
        long secondNewDyad = dyadKey(c, b);
        int before = reciprocalContribution(
                edgeCounts, firstOldDyad, secondOldDyad, firstNewDyad, secondNewDyad);
        int after = projectedReciprocalContribution(
                edgeCounts, oldFirst, oldSecond, newFirst, newSecond,
                firstOldDyad, secondOldDyad, firstNewDyad, secondNewDyad);
        int updatedReciprocalArcs = reciprocalArcs + after - before;
        boolean accept = requireIncrease
                ? updatedReciprocalArcs > reciprocalArcs
                : updatedReciprocalArcs == reciprocalArcs;

        if (!accept) {
            return reciprocalArcs;
        }

        decrementEdgeCount(edgeCounts, oldFirst);
        decrementEdgeCount(edgeCounts, oldSecond);
        incrementEdgeCount(edgeCounts, newFirst);
        incrementEdgeCount(edgeCounts, newSecond);
        dsts[first] = d;
        dsts[second] = b;
        return updatedReciprocalArcs;
    }

    static boolean tryRandomDegreePreservingSwap(
            int[] srcs, int[] dsts, Map<Long, Integer> edgeCounts, Random random) {
        int edgeCount = srcs.length;
        if (edgeCount < 2) {
            return false;
        }

        int first = random.nextInt(edgeCount);
        int second = random.nextInt(edgeCount - 1);
        if (second >= first) {
            second++;
        }

        int a = srcs[first];
        int b = dsts[first];
        int c = srcs[second];
        int d = dsts[second];
        if (a == d || c == b || a == c || b == d) {
            return false;
        }

        long oldFirst = edgeKey(a, b);
        long oldSecond = edgeKey(c, d);
        long newFirst = edgeKey(a, d);
        long newSecond = edgeKey(c, b);
        if (newFirst == newSecond
                || edgeCountAfterRemoval(edgeCounts, newFirst, oldFirst, oldSecond) > 0
                || edgeCountAfterRemoval(edgeCounts, newSecond, oldFirst, oldSecond) > 0) {
            return false;
        }

        decrementEdgeCount(edgeCounts, oldFirst);
        decrementEdgeCount(edgeCounts, oldSecond);
        incrementEdgeCount(edgeCounts, newFirst);
        incrementEdgeCount(edgeCounts, newSecond);
        dsts[first] = d;
        dsts[second] = b;
        return true;
    }

    static void shuffle(int[] values, Random random) {
        for (int i = values.length - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int value = values[i];
            values[i] = values[j];
            values[j] = value;
        }
    }

    private Map<Long, Integer> buildEdgeCounts() {
        Map<Long, Integer> edgeCounts = new HashMap<>();
        for (int u = 0; u < n; u++) {
            for (int i = outPtr[u]; i < outPtr[u + 1]; i++) {
                incrementEdgeCount(edgeCounts, edgeKey(u, outIdx[i]));
            }
        }
        return edgeCounts;
    }

    static Map<Long, Integer> buildEdgeCounts(int[] srcs, int[] dsts) {
        Map<Long, Integer> edgeCounts = new HashMap<>();
        for (int i = 0; i < srcs.length; i++) {
            incrementEdgeCount(edgeCounts, edgeKey(srcs[i], dsts[i]));
        }
        return edgeCounts;
    }

    static int countReciprocalArcs(Map<Long, Integer> edgeCounts) {
        int reciprocalArcs = 0;
        for (Map.Entry<Long, Integer> entry : edgeCounts.entrySet()) {
            int u = sourceFromKey(entry.getKey());
            int v = destinationFromKey(entry.getKey());
            if (u < v) {
                reciprocalArcs += 2 * Math.min(
                        entry.getValue(), edgeCounts.getOrDefault(edgeKey(v, u), 0));
            }
        }
        return reciprocalArcs;
    }

    static int reciprocalContribution(Map<Long, Integer> edgeCounts,
            long first, long second, long third, long fourth) {
        int contribution = reciprocalContribution(edgeCounts, first);
        if (second != first) {
            contribution += reciprocalContribution(edgeCounts, second);
        }
        if (third != first && third != second) {
            contribution += reciprocalContribution(edgeCounts, third);
        }
        if (fourth != first && fourth != second && fourth != third) {
            contribution += reciprocalContribution(edgeCounts, fourth);
        }
        return contribution;
    }

    static int reciprocalContribution(Map<Long, Integer> edgeCounts, long dyad) {
        int u = sourceFromKey(dyad);
        int v = destinationFromKey(dyad);
        if (u == v) {
            return 0;
        }
        return 2 * Math.min(
                edgeCounts.getOrDefault(edgeKey(u, v), 0),
                edgeCounts.getOrDefault(edgeKey(v, u), 0));
    }

    static int projectedReciprocalContribution(
            Map<Long, Integer> edgeCounts,
            long oldFirst, long oldSecond, long newFirst, long newSecond,
            long firstDyad, long secondDyad, long thirdDyad, long fourthDyad) {
        int contribution = projectedReciprocalContribution(
                edgeCounts, oldFirst, oldSecond, newFirst, newSecond, firstDyad);
        if (secondDyad != firstDyad) {
            contribution += projectedReciprocalContribution(
                    edgeCounts, oldFirst, oldSecond, newFirst, newSecond, secondDyad);
        }
        if (thirdDyad != firstDyad && thirdDyad != secondDyad) {
            contribution += projectedReciprocalContribution(
                    edgeCounts, oldFirst, oldSecond, newFirst, newSecond, thirdDyad);
        }
        if (fourthDyad != firstDyad && fourthDyad != secondDyad && fourthDyad != thirdDyad) {
            contribution += projectedReciprocalContribution(
                    edgeCounts, oldFirst, oldSecond, newFirst, newSecond, fourthDyad);
        }
        return contribution;
    }

    static int projectedReciprocalContribution(
            Map<Long, Integer> edgeCounts,
            long oldFirst, long oldSecond, long newFirst, long newSecond, long dyad) {
        int u = sourceFromKey(dyad);
        int v = destinationFromKey(dyad);
        if (u == v) {
            return 0;
        }
        return 2 * Math.min(
                projectedEdgeCount(edgeCounts, edgeKey(u, v), oldFirst, oldSecond, newFirst, newSecond),
                projectedEdgeCount(edgeCounts, edgeKey(v, u), oldFirst, oldSecond, newFirst, newSecond));
    }

    static int edgeCountAfterRemoval(
            Map<Long, Integer> edgeCounts, long key, long oldFirst, long oldSecond) {
        int count = edgeCounts.getOrDefault(key, 0);
        if (key == oldFirst) {
            count--;
        }
        if (key == oldSecond) {
            count--;
        }
        return count;
    }

    static int projectedEdgeCount(
            Map<Long, Integer> edgeCounts, long key,
            long oldFirst, long oldSecond, long newFirst, long newSecond) {
        int count = edgeCountAfterRemoval(edgeCounts, key, oldFirst, oldSecond);
        if (key == newFirst) {
            count++;
        }
        if (key == newSecond) {
            count++;
        }
        return count;
    }

    static void incrementEdgeCount(Map<Long, Integer> edgeCounts, long key) {
        edgeCounts.merge(key, 1, Integer::sum);
    }

    static void decrementEdgeCount(Map<Long, Integer> edgeCounts, long key) {
        int count = edgeCounts.getOrDefault(key, 0);
        if (count <= 0) {
            throw new IllegalStateException("Cannot remove an edge that does not exist");
        }
        if (count == 1) {
            edgeCounts.remove(key);
        } else {
            edgeCounts.put(key, count - 1);
        }
    }

    static long edgeKey(int source, int destination) {
        return ((long) source << 32) | (destination & 0xffffffffL);
    }

    static long dyadKey(int first, int second) {
        return edgeKey(Math.min(first, second), Math.max(first, second));
    }

    static int sourceFromKey(long key) {
        return (int) (key >>> 32);
    }

    static int destinationFromKey(long key) {
        return (int) key;
    }

}
