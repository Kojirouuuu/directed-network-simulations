package sim.network;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * 有向グラフをCSR（Compressed Sparse Row）形式で表現するクラス。
 * この実装では、有向辺(u, v)と(v, u)が同時に存在し、かつ無向辺(u, v)が存在する場合を含みます。
 */
public final class DirectedGraph {
    /** ランダム再接続の前に頂点間で再配置する次数列。 */
    public enum DegreeSide {
        IN,
        OUT
    }

    public final String name; // グラフの名前（RR, ER, BA, CM, ...）
    public final int n; // 頂点数
    public final int m; // 辺数
    private final int[] outPtr; // 長さ n+1 の配列（out-neighbors の開始位置）
    private final int[] outIdx; // 長さ m の配列（out-neighbors の頂点ID）
    private final boolean[] outIsUndirected; // 長さ m の配列（true なら無向由来）
    private final int[] inPtr; // 長さ n+1 の配列（in-neighbors の開始位置）
    private final int[] inIdx; // 長さ m の配列（in-neighbors の頂点ID）
    private final boolean[] inIsUndirected; // 長さ m の配列（true なら無向由来）
    private volatile long cachedFeedForwardLoopCount = -1L;

    private DirectedGraph(String name, int n, int m, int[] outPtr, int[] outIdx, boolean[] outIsUndirected,
            int[] inPtr, int[] inIdx, boolean[] inIsUndirected) {
        this.name = name;
        this.n = n;
        this.m = m;
        this.outPtr = outPtr;
        this.outIdx = outIdx;
        this.outIsUndirected = outIsUndirected;
        this.inPtr = inPtr;
        this.inIdx = inIdx;
        this.inIsUndirected = inIsUndirected;

        validate();
    }

    /**
     * インデックス i の out-neighbor を取得する。
     *
     * @param i インデックス（0 以上 m 未満）
     * @return 隣接頂点のID
     */
    public int getOutNeighbor(int i) {
        if (i < 0 || i >= m) {
            throw new IndexOutOfBoundsException("Index: " + i + ", Size: " + m);
        }
        return outIdx[i];
    }

    /**
     * インデックス i の in-neighbor を取得する。
     *
     * @param i インデックス（0 以上 m 未満）
     * @return 隣接頂点のID
     */
    public int getInNeighbor(int i) {
        if (i < 0 || i >= m) {
            throw new IndexOutOfBoundsException("Index: " + i + ", Size: " + m);
        }
        return inIdx[i];
    }

    /**
     * インデックス i の out-neighbor が無向由来かどうかを判定する。
     *
     * @param i インデックス（0 以上 m 未満）
     * @return 無向由来の場合 true
     */
    public boolean isOutUndirected(int i) {
        if (i < 0 || i >= m) {
            throw new IndexOutOfBoundsException("Index: " + i + ", Size: " + m);
        }
        return outIsUndirected[i];
    }

    /**
     * インデックス i の in-neighbor が無向由来かどうかを判定する。
     *
     * @param i インデックス（0 以上 m 未満）
     * @return 無向由来の場合 true
     */
    public boolean isInUndirected(int i) {
        if (i < 0 || i >= m) {
            throw new IndexOutOfBoundsException("Index: " + i + ", Size: " + m);
        }
        return inIsUndirected[i];
    }

    /**
     * 頂点 u の out-neighbors の範囲を取得する。
     * outIdx の [start, end) が u の out-neighbors を表す。
     *
     * @param u 頂点ID（0 以上 n 未満）
     * @return 範囲 [start, end)
     */
    public IntRange outNeighborRange(int u) {
        if (u < 0 || u >= n) {
            throw new IndexOutOfBoundsException("Vertex: " + u + ", Range: [0, " + n + ")");
        }
        return new IntRange(outPtr[u], outPtr[u + 1]);
    }

    /**
     * 頂点 u の in-neighbors の範囲を取得する。
     * inIdx の [start, end) が u の in-neighbors を表す。
     *
     * @param u 頂点ID（0 以上 n 未満）
     * @return 範囲 [start, end)
     */
    public IntRange inNeighborRange(int u) {
        if (u < 0 || u >= n) {
            throw new IndexOutOfBoundsException("Vertex: " + u + ", Range: [0, " + n + ")");
        }
        return new IntRange(inPtr[u], inPtr[u + 1]);
    }

    /**
     * 有向辺を全て逆向きにした新しいグラフを返す。
     * 元の辺 (u, v) は (v, u) になる。無向辺はそのまま（両方向のまま）扱われる。
     * 頂点数・辺数は変わらない。
     *
     * @return 有向辺を反転した新しい DirectedGraph（名前は name + "_reversed"）
     */
    public DirectedGraph reverseDirectedEdges() {
        return reverseDirectedEdges(name + "_reversed");
    }

    /**
     * 有向辺を全て逆向きにした新しいグラフを返す（名前指定版）。
     * 元の辺 (u, v) は (v, u) になる。無向辺はそのまま扱われる。
     *
     * @param newName 返すグラフの名前（null の場合は name + "_reversed"）
     * @return 有向辺を反転した新しい DirectedGraph
     */
    public DirectedGraph reverseDirectedEdges(String newName) {
        String resultName = newName != null ? newName : name + "_reversed";
        return new DirectedGraph(
                resultName,
                n,
                m,
                inPtr,
                inIdx,
                inIsUndirected,
                outPtr,
                outIdx,
                outIsUndirected);
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
        long cached = cachedFeedForwardLoopCount;
        if (cached >= 0L) {
            return cached;
        }

        FflRewiringState state = new FflRewiringState(this);
        long count = state.countAllFeedForwardLoops();
        cachedFeedForwardLoopCount = count;
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

        FflRewiringState state = new FflRewiringState(this);
        long initialNffl = cachedFeedForwardLoopCount;
        if (initialNffl < 0L) {
            initialNffl = state.countAllFeedForwardLoops();
            cachedFeedForwardLoopCount = initialNffl;
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

    private static int tryDegreePreservingSwap(int[] srcs, int[] dsts, Map<Long, Integer> edgeCounts,
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

    private static boolean tryRandomDegreePreservingSwap(
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

    private static void shuffle(int[] values, Random random) {
        for (int i = values.length - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int value = values[i];
            values[i] = values[j];
            values[j] = value;
        }
    }

    /** Mutable edge-index view used while searching for FFL-increasing swaps. */
    private static final class FflRewiringState {
        private final int vertexCount;
        private final int edgeCount;
        private final int[] outEdgePtr;
        private final int[] inEdgePtr;
        private final int[] sources;
        private final int[] destinations;
        private final int[] inEdgeIndices;
        private final int[] inPositionByEdge;
        private final boolean[] representativeEdge;
        private final Map<Long, Integer> edgeCounts;
        private final long[] cumulativeWedgeWeights;
        private final long totalWedgeWeight;
        private long feedForwardLoopCount;

        private FflRewiringState(DirectedGraph graph) {
            vertexCount = graph.n;
            edgeCount = graph.m;
            outEdgePtr = graph.outPtr;
            inEdgePtr = graph.inPtr;
            sources = new int[edgeCount];
            destinations = Arrays.copyOf(graph.outIdx, edgeCount);
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

        private long countAllFeedForwardLoops() {
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

        private boolean tryIncreasingSwap(Random random) {
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

        private int sampleMiddleVertex(Random random) {
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

        private void swapDestinations(int firstEdge, int secondEdge) {
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

        private void removeEdge(int edge) {
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

        private void addEdge(int edge) {
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

        private int findParallelEdge(int source, int destination, int excludedEdge) {
            for (int edge = outEdgePtr[source]; edge < outEdgePtr[source + 1]; edge++) {
                if (edge != excludedEdge && destinations[edge] == destination) {
                    return edge;
                }
            }
            return -1;
        }

        private long countFeedForwardLoopsContaining(int source, int destination) {
            if (source == destination) {
                return 0L;
            }
            long count = countCommonOutgoing(source, destination);
            count = Math.addExact(count, countCommonIncoming(source, destination));
            return Math.addExact(count, countOutgoingIncomingIntersection(source, destination));
        }

        private long countCommonOutgoing(int first, int second) {
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

        private long countCommonIncoming(int first, int second) {
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

        private long countOutgoingIncomingIntersection(int source, int destination) {
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

        private boolean edgeExists(int source, int destination) {
            return edgeCounts.getOrDefault(edgeKey(source, destination), 0) > 0;
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

    private static Map<Long, Integer> buildEdgeCounts(int[] srcs, int[] dsts) {
        Map<Long, Integer> edgeCounts = new HashMap<>();
        for (int i = 0; i < srcs.length; i++) {
            incrementEdgeCount(edgeCounts, edgeKey(srcs[i], dsts[i]));
        }
        return edgeCounts;
    }

    private static int countReciprocalArcs(Map<Long, Integer> edgeCounts) {
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

    private static int reciprocalContribution(Map<Long, Integer> edgeCounts,
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

    private static int reciprocalContribution(Map<Long, Integer> edgeCounts, long dyad) {
        int u = sourceFromKey(dyad);
        int v = destinationFromKey(dyad);
        if (u == v) {
            return 0;
        }
        return 2 * Math.min(
                edgeCounts.getOrDefault(edgeKey(u, v), 0),
                edgeCounts.getOrDefault(edgeKey(v, u), 0));
    }

    private static int projectedReciprocalContribution(
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

    private static int projectedReciprocalContribution(
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

    private static int edgeCountAfterRemoval(
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

    private static int projectedEdgeCount(
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

    private static void incrementEdgeCount(Map<Long, Integer> edgeCounts, long key) {
        edgeCounts.merge(key, 1, Integer::sum);
    }

    private static void decrementEdgeCount(Map<Long, Integer> edgeCounts, long key) {
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

    private static long edgeKey(int source, int destination) {
        return ((long) source << 32) | (destination & 0xffffffffL);
    }

    private static long dyadKey(int first, int second) {
        return edgeKey(Math.min(first, second), Math.max(first, second));
    }

    private static int sourceFromKey(long key) {
        return (int) (key >>> 32);
    }

    private static int destinationFromKey(long key) {
        return (int) key;
    }

    /**
     * 有向辺リストから DirectedGraph を構築する。
     *
     * @param name グラフの名前
     * @param n 頂点数
     * @param srcs 始点の配列
     * @param dsts 終点の配列
     * @param isUndirected 各辺が無向かどうかの配列
     * @return 構築された DirectedGraph
     */
    public static DirectedGraph fromEdgeListWithUndirectedFlag(String name, int n, int[] srcs, int[] dsts,
            boolean[] isUndirected) {
        if (srcs.length != dsts.length) {
            throw new IllegalArgumentException("Source and destination arrays must have the same length");
        }
        if (isUndirected == null || isUndirected.length != srcs.length) {
            throw new IllegalArgumentException("isUndirected array length must match source array length");
        }

        final int mBase = srcs.length;

        // 展開後の有向辺数を計算（無向辺は逆向きが追加される）
        int undCount = 0;
        for (boolean b : isUndirected) {
            if (b) {
                undCount++;
            }
        }
        final int mExp = mBase + undCount;

        int[] outDeg = new int[n];
        int[] inDeg = new int[n];

        for (int i = 0; i < mBase; i++) {
            int u = srcs[i];
            int v = dsts[i];
            if (u < 0 || u >= n || v < 0 || v >= n) {
                throw new IllegalArgumentException("Invalid edge: (" + u + ", " + v + ")");
            }
            outDeg[u]++;
            inDeg[v]++;

            if (isUndirected[i]) {
                outDeg[v]++;
                inDeg[u]++;
            }
        }

        int[] outPtr = new int[n + 1];
        int[] outPos = new int[n];
        int[] inPtr = new int[n + 1];
        int[] inPos = new int[n];

        for (int u = 0; u < n; u++) {
            outPtr[u + 1] = outPtr[u] + outDeg[u];
            outPos[u] = outPtr[u];
            inPtr[u + 1] = inPtr[u] + inDeg[u];
            inPos[u] = inPtr[u];
        }

        int[] outIdx = new int[mExp];
        int[] inIdx = new int[mExp];
        boolean[] outIsUndirected = new boolean[mExp];
        boolean[] inIsUndirected = new boolean[mExp];

        for (int i = 0; i < mBase; i++) {
            int u = srcs[i];
            int v = dsts[i];
            outIdx[outPos[u]] = v;
            inIdx[inPos[v]] = u;

            if (isUndirected[i]) {
                outIdx[outPos[v]] = u;
                inIdx[inPos[u]] = v;

                outIsUndirected[outPos[u]] = true;
                inIsUndirected[inPos[v]] = true;
                outIsUndirected[outPos[v]] = true;
                inIsUndirected[inPos[u]] = true;

                outPos[v]++;
                inPos[u]++;
            }

            outPos[u]++;
            inPos[v]++;
        }

        return new DirectedGraph(name, n, mExp, outPtr, outIdx, outIsUndirected, inPtr, inIdx, inIsUndirected);
    }

    /**
     * エッジリストファイルを読み、有向グラフを構築する。
     * 各行を "a b" と解釈して有向辺 a→b を追加する。空行および '#' で始まる行は無視する。
     * 頂点IDは出現順に 0..n-1 にマッピングする。
     *
     * @param name グラフ名（null の場合はファイル名を使用）
     * @param path エッジリストファイルのパス
     * @return 構築された DirectedGraph
     * @throws IOException ファイル読み込みエラー
     */
    public static DirectedGraph loadFromEdgeList(String name, Path path) throws IOException {
        if (path == null) {
            throw new IllegalArgumentException("path must be non-null");
        }
        List<int[]> rawEdges = new ArrayList<>();
        Map<String, Integer> idToIndex = new LinkedHashMap<>();

        try (BufferedReader reader = Files.newBufferedReader(path)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                String[] parts = line.split("\\s+");
                if (parts.length < 2) {
                    continue;
                }
                String srcId = parts[0];
                String dstId = parts[1];
                int u = idToIndex.computeIfAbsent(srcId, k -> idToIndex.size());
                int v = idToIndex.computeIfAbsent(dstId, k -> idToIndex.size());
                rawEdges.add(new int[] { u, v });
            }
        }

        int n = idToIndex.size();
        int m = rawEdges.size();
        int[] srcs = new int[m];
        int[] dsts = new int[m];
        boolean[] isUndirected = new boolean[m];

        for (int i = 0; i < m; i++) {
            int[] e = rawEdges.get(i);
            srcs[i] = e[0];
            dsts[i] = e[1];
            isUndirected[i] = false;
        }

        String graphName = name != null ? name : path.getFileName().toString();
        return fromEdgeListWithUndirectedFlag(graphName, n, srcs, dsts, isUndirected);
    }

    /**
     * CSR 形式の整合性をチェックする。計算量は O(n + m)。
     */
    public void validate() {
        requireArrayLen(outPtr, n + 1, "outPtr");
        requireArrayLen(outIdx, m, "outIdx");
        requireArrayLen(inPtr, n + 1, "inPtr");
        requireArrayLen(inIdx, m, "inIdx");

        validatePtr(outPtr, "outPtr");
        validatePtr(inPtr, "inPtr");

        if (outPtr[n] != m) {
            throw new IllegalStateException("outPtr[" + n + "] must be equal to m, but was " + outPtr[n]);
        }
        if (inPtr[n] != m) {
            throw new IllegalStateException("inPtr[" + n + "] must be equal to m, but was " + inPtr[n]);
        }

        validateUndirectedEdges();
    }

    private static void requireArrayLen(int[] arr, int len, String name) {
        if (arr == null) {
            throw new IllegalStateException(name + " must be non-null");
        }
        if (arr.length != len) {
            throw new IllegalStateException(name + " must have length " + len + ", but was " + arr.length);
        }
    }

    private void validatePtr(int[] ptr, String name) {
        if (ptr[0] != 0) {
            throw new IllegalStateException(name + "[0] must be 0, but was " + ptr[0]);
        }
        for (int i = 0; i < n; i++) {
            if (ptr[i + 1] < ptr[i]) {
                throw new IllegalStateException(name + " must be non-decreasing at index " + i);
            }
        }
    }

    private void validateUndirectedEdges() {
        int numberOfUndirectedOut = 0;
        int numberOfUndirectedIn = 0;
        for (int i = 0; i < m; i++) {
            if (outIsUndirected[i]) {
                numberOfUndirectedOut++;
            }
            if (inIsUndirected[i]) {
                numberOfUndirectedIn++;
            }
        }
        if (numberOfUndirectedOut != numberOfUndirectedIn) {
            throw new IllegalStateException(
                    "Number of undirected out-edges (" + numberOfUndirectedOut
                            + ") must equal number of undirected in-edges (" + numberOfUndirectedIn + ")");
        }
    }

    /**
     * 次数の合計が正しいかチェックする。
     */
    public void checkDegreeSums() {
        long outSum = 0;
        long inSum = 0;
        for (int u = 0; u < n; u++) {
            outSum += (outPtr[u + 1] - outPtr[u]);
            inSum += (inPtr[u + 1] - inPtr[u]);
        }
        if (outSum != m) {
            throw new IllegalStateException("Sum of out-degrees must be equal to m, but was " + outSum);
        }
        if (inSum != m) {
            throw new IllegalStateException("Sum of in-degrees must be equal to m, but was " + inSum);
        }
    }

    /**
     * グラフの最大連結成分のサイズを返す。
     *
     * @param includeUndirectedEdges 無向辺を含むかどうか
     * @return 最大連結成分のサイズ
     */
    public int checkConnected(boolean includeUndirectedEdges) {
        if (n == 0) {
            return 0;
        }

        boolean[] visited = new boolean[n];
        ArrayDeque<Integer> q = new ArrayDeque<>();

        int maxSize = 0;
        for (int start = 0; start < n; start++) {
            if (visited[start]) {
                continue;
            }
            visited[start] = true;
            boolean[] curVisited = new boolean[n];
            curVisited[start] = true;

            q.add(start);
            int seen = 1;

            while (!q.isEmpty()) {
                int u = q.poll();

                // out-neighbors を辿る
                IntRange rOut = outNeighborRange(u);
                for (int i = rOut.start; i < rOut.end; i++) {
                    if (!includeUndirectedEdges && isOutUndirected(i)) {
                        continue;
                    }
                    int v = getOutNeighbor(i);
                    if (!curVisited[v]) {
                        visited[v] = true;
                        curVisited[v] = true;
                        q.add(v);
                        seen++;
                    }
                }

                // in-neighbors を辿る
                IntRange rIn = inNeighborRange(u);
                for (int i = rIn.start; i < rIn.end; i++) {
                    if (!includeUndirectedEdges && isInUndirected(i)) {
                        continue;
                    }
                    int v = getInNeighbor(i);
                    if (!curVisited[v]) {
                        visited[v] = true;
                        curVisited[v] = true;
                        q.add(v);
                        seen++;
                    }
                }
            }

            maxSize = Math.max(maxSize, seen);
        }

        return maxSize;
    }

    /**
     * グラフの情報を標準出力に表示する。
     */
    public void printInfo() {
        System.out.println("=== Graph Info ===");
        System.out.println("Name: " + name);
        System.out.println("Vertices (n): " + n);

        int numberOfDirectedArcs = 0;
        int numberOfUndirectedArcs = 0;
        for (int i = 0; i < m; i++) {
            if (outIsUndirected[i]) {
                numberOfUndirectedArcs++;
            } else {
                numberOfDirectedArcs++;
            }
        }
        System.out.println("Directed Arcs (m): " + numberOfDirectedArcs);
        System.out.println("Undirected Arcs (m): " + numberOfUndirectedArcs);
        System.out.println(
                "Note: Even if directed edges (u, v) and (v, u) both exist, this does not mean there exists an undirected edge (u, v).");

        double inAverageDegreeExcludingUndirected = 0;
        double inAverageDegreeIncludingUndirected = 0;
        double outAverageDegreeExcludingUndirected = 0;
        double outAverageDegreeIncludingUndirected = 0;

        int maxInDegree = 0;
        int maxOutDegree = 0;
        int minInDegree = Integer.MAX_VALUE;
        int minOutDegree = Integer.MAX_VALUE;

        for (int u = 0; u < n; u++) {
            int inDirected = 0;
            int inRange = inPtr[u + 1] - inPtr[u];
            for (int i = 0; i < inRange; i++) {
                if (!isInUndirected(inPtr[u] + i)) {
                    inDirected++;
                    inAverageDegreeExcludingUndirected++;
                }
                inAverageDegreeIncludingUndirected++;
            }
            maxInDegree = Math.max(maxInDegree, inDirected);
            minInDegree = Math.min(minInDegree, inDirected);

            int outDirected = 0;
            int outRange = outPtr[u + 1] - outPtr[u];
            for (int i = 0; i < outRange; i++) {
                if (!isOutUndirected(outPtr[u] + i)) {
                    outDirected++;
                    outAverageDegreeExcludingUndirected++;
                }
                outAverageDegreeIncludingUndirected++;
            }
            maxOutDegree = Math.max(maxOutDegree, outDirected);
            minOutDegree = Math.min(minOutDegree, outDirected);
        }
        if (n == 0) {
            minInDegree = 0;
            minOutDegree = 0;
        }
        inAverageDegreeExcludingUndirected /= n;
        inAverageDegreeIncludingUndirected /= n;
        outAverageDegreeExcludingUndirected /= n;
        outAverageDegreeIncludingUndirected /= n;

        System.out.println("");
        System.out.println("Average in-degree (excluding undirected edges): " + inAverageDegreeExcludingUndirected);
        System.out.println("Average in-degree (including undirected edges): " + inAverageDegreeIncludingUndirected);
        System.out.println("Average out-degree (excluding undirected edges): " + outAverageDegreeExcludingUndirected);
        System.out.println("Average out-degree (including undirected edges): " + outAverageDegreeIncludingUndirected);

        System.out.println("");
        System.out.println("Max in-degree: " + maxInDegree);
        System.out.println("Min in-degree: " + minInDegree);
        System.out.println("Max out-degree: " + maxOutDegree);
        System.out.println("Min out-degree: " + minOutDegree);

        // 最大弱連結成分（WCC）
        int maxWccDirectedOnly = checkConnected(false);
        int maxWccIncludingUndirected = checkConnected(true);

        System.out.println("");
        System.out.println("Largest WCC size (directed-only edges): " + maxWccDirectedOnly);
        System.out.println("Largest WCC size (including undirected edges): " + maxWccIncludingUndirected);
    }

    /**
     * Pythonのnetworkxで読み込める形式で辺リストをファイルに書き出す。
     * 出力形式は "u v"（スペース区切り）で、networkx.read_edgelist()で読み込める。
     *
     * @param path 出力先のファイルパス
     * @throws IOException ファイル書き込みエラー
     */
    public void writeEdgeList(Path path) throws IOException {
        Files.createDirectories(path.getParent());
        try (BufferedWriter bw = Files.newBufferedWriter(path, StandardOpenOption.CREATE,
                StandardOpenOption.TRUNCATE_EXISTING);
                PrintWriter out = new PrintWriter(bw)) {
            for (int u = 0; u < n; u++) {
                IntRange rOut = outNeighborRange(u);
                for (int i = rOut.start; i < rOut.end; i++) {
                    int v = getOutNeighbor(i);
                    out.println(u + " " + v);
                }
            }
        }
    }

    /**
     * 整数の範囲 [start, end) を表すクラス。
     */
    public static final class IntRange {
        /** 開始位置（含む） */
        public final int start;

        /** 終了位置（含まない） */
        public final int end;

        public IntRange(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }
}
