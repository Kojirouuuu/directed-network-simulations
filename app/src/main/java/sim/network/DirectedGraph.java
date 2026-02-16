package sim.network;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;

/**
 * 有向グラフをCSR（Compressed Sparse Row）形式で表現するクラス。
 * この実装では、有向辺(u, v)と(v, u)が同時に存在し、かつ無向辺(u, v)が存在する場合を含みます。
 */
public final class DirectedGraph {
    public final String name; // グラフの名前（RR, ER, BA, CM, ...）
    public final int n; // 頂点数
    public final int m; // 辺数
    private final int[] outPtr; // 長さ n+1 の配列（out-neighbors の開始位置）
    private final int[] outIdx; // 長さ m の配列（out-neighbors の頂点ID）
    private final boolean[] outIsUndirected; // 長さ m の配列（true なら無向由来）
    private final int[] inPtr; // 長さ n+1 の配列（in-neighbors の開始位置）
    private final int[] inIdx; // 長さ m の配列（in-neighbors の頂点ID）
    private final boolean[] inIsUndirected; // 長さ m の配列（true なら無向由来）
    
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
