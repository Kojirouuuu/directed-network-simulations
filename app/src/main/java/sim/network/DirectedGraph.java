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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

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
    volatile long cachedFeedForwardLoopCount = -1L;

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
        return new GraphOperations(this).countReciprocalArcs();
    }

    /**
     * 全アークに占める相互アークの割合を返す。
     *
     * @return 相互アーク数 / m。空グラフでは 0.0
     */
    public double reciprocity() {
        return new GraphOperations(this).reciprocity();
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
        return new GraphOperations(this).countFeedForwardLoops();
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
        return new GraphOperations(this).increaseFeedForwardLoops(targetNffl, maxIncreaseAttempts, seed);
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
        return new GraphOperations(this).increaseReciprocity(targetReciprocity, maxIncreaseAttempts, neutralSwapAttempts, seed);
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
        return new GraphOperations(this).randomizeByShuffledDegreeSequence(side, seed);
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
        return new GraphOperations(this).expandByRepeatedJointDegreeSequence(multiplier, seed);
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
        return new GraphOperations(this).randomizeByEdgeSwaps(seed);
    }

    record CsrView(int[] outPtr, int[] outIdx, int[] inPtr, boolean[] outIsUndirected) {}

    CsrView csrView() {
        return new CsrView(outPtr, outIdx, inPtr, outIsUndirected);
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
