package sim.network.topology;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import sim.network.DirectedGraph;
import sim.network.DirectedGraph.IntRange;

/**
 * {@link SchwartzDirectedSF} の不変条件と失敗パスを検証するテスト。
 *
 * <p>
 * 注意: 辺数 Σ K[i] が {@code Integer.MAX_VALUE} を超えるケースは、
 * 数十億ノード規模の {@code int[]} を必要とするため JUnit では再現不能。
 * オーバーフロー検査自体は実装にあるので、コードレビューで担保する。
 */
class SchwartzDirectedSFTest {

    private static long sumInDeg(DirectedGraph g) {
        long s = 0;
        for (int u = 0; u < g.n; u++) {
            IntRange r = g.inNeighborRange(u);
            s += (r.end - r.start);
        }
        return s;
    }

    private static long sumOutDeg(DirectedGraph g) {
        long s = 0;
        for (int u = 0; u < g.n; u++) {
            IntRange r = g.outNeighborRange(u);
            s += (r.end - r.start);
        }
        return s;
    }

    @Test
    void normalCaseSucceeds() {
        int n = 10_000;
        int mIn = 1, kInMax = 100, mOut = 1, kOutMax = 100;
        double gammaIn = 2.5, gammaOut = 3.0;

        for (double corrA : new double[] { 0.0, 0.5, 1.0 }) {
            DirectedGraph g = SchwartzDirectedSF.generate(
                    "test", n, mIn, kInMax, mOut, kOutMax,
                    gammaIn, gammaOut, corrA, 42L);
            assertEquals(n, g.n, "graph vertex count mismatch (corrA=" + corrA + ")");
            assertEquals(sumInDeg(g), sumOutDeg(g),
                    "sum(inDeg) must equal sum(outDeg) (corrA=" + corrA + ")");
        }
    }

    @Test
    void outDegreeWithinBounds() {
        int n = 10_000;
        int mIn = 1, kInMax = 100, mOut = 1, kOutMax = 100;
        DirectedGraph g = SchwartzDirectedSF.generate(
                "test", n, mIn, kInMax, mOut, kOutMax,
                2.5, 3.0, 0.5, 7L);
        for (int u = 0; u < g.n; u++) {
            IntRange r = g.outNeighborRange(u);
            int outDeg = r.end - r.start;
            assertTrue(outDeg >= mOut && outDeg <= kOutMax,
                    "outDeg[" + u + "]=" + outDeg + " outside [" + mOut + ", " + kOutMax + "]");
        }
    }

    @Test
    void inDegreeWithinBounds() {
        int n = 10_000;
        int mIn = 1, kInMax = 100, mOut = 1, kOutMax = 100;
        DirectedGraph g = SchwartzDirectedSF.generate(
                "test", n, mIn, kInMax, mOut, kOutMax,
                2.5, 3.0, 0.5, 9L);
        for (int u = 0; u < g.n; u++) {
            IntRange r = g.inNeighborRange(u);
            int inDeg = r.end - r.start;
            assertTrue(inDeg >= mIn && inDeg <= kInMax,
                    "inDeg[" + u + "]=" + inDeg + " outside [" + mIn + ", " + kInMax + "]");
        }
    }

    @Test
    void infeasibleSumThrowsWhenJiTooHigh() {
        // mIn=10, kInMax=10 → sum(ji) = n*10 = 100 (固定)
        // mOut=1, kOutMax=5 → sum(ki) は最大 n*5 = 50 までしか取れない
        // よって sum(ji) > n*kOutMax で例外
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 10, 10, 10, 1, 5, 2.5, 3.0, 0.0, 1L));
    }

    @Test
    void infeasibleSumThrowsWhenJiTooLow() {
        // mIn=1, kInMax=1 → sum(ji) = n*1 = 10 (固定)
        // mOut=5, kOutMax=10 → sum(ki) は最小 n*5 = 50 から
        // よって sum(ji) < n*mOut で例外
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 10, 1, 1, 5, 10, 2.5, 3.0, 0.0, 1L));
    }

    @Test
    void invalidArgumentsThrow() {
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                null, 100, 1, 10, 1, 10, 2.5, 3.0, 0.5, 1L), "null name");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 0, 1, 10, 1, 10, 2.5, 3.0, 0.5, 1L), "n < 1");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 0, 10, 1, 10, 2.5, 3.0, 0.5, 1L), "mIn < 1");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 0, 10, 2.5, 3.0, 0.5, 1L), "mOut < 1");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 5, 4, 1, 10, 2.5, 3.0, 0.5, 1L), "kInMax < mIn");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 5, 4, 2.5, 3.0, 0.5, 1L), "kOutMax < mOut");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 1, 10, 1.0, 3.0, 0.5, 1L), "gammaIn <= 1");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 1, 10, 2.5, 1.0, 0.5, 1L), "gammaOut <= 1");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 1, 10, 2.5, 3.0, -0.1, 1L), "corrA < 0");
        assertThrows(IllegalArgumentException.class, () -> SchwartzDirectedSF.generate(
                "test", 100, 1, 10, 1, 10, 2.5, 3.0, 1.1, 1L), "corrA > 1");
    }
}
