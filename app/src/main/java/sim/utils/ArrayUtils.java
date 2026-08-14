package sim.utils;

import java.util.Random;

/**
 * 配列操作のユーティリティクラス。
 */
public final class ArrayUtils {
    private ArrayUtils() {
    }

    /**
     * Fisher-Yates シャッフルアルゴリズムで配列をランダムに並び替える。
     *
     * @param arr  シャッフルする配列
     * @param seed 乱数シード
     */
    public static void shuffle(int[] arr, long seed) {
        Random random = new Random(seed);
        for (int i = arr.length - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }

    /**
     * 指定された範囲とステップで整数配列を生成する。
     *
     * @param start 開始値（含む）
     * @param end   終了値（含む）
     * @param step  ステップ
     * @return 生成された配列
     */
    public static int[] arange(int start, int end, int step) {
        int length = (int)((end - start) / step) + 1;
        int[] array = new int[length];
        for (int i = 0; i < length; i++) {
            array[i] = start + i * step;
        }
        return array;
    }
    /**
     * 指定された範囲とステップで実数配列を生成する。
     *
     * @param start 開始値（含む）
     * @param end   終了値（含む）
     * @param step  ステップ
     * @return 生成された配列
     */
    public static double[] arange(double start, double end, double step) {
        int length = (int)((end - start) / step) + 1;
        double[] array = new double[length];
        for (int i = 0; i < length; i++) {
            array[i] = start + i * step;
        }
        return array;
    }

    /**
     * 指定された正の開始値と終了値の間を、等比間隔で分割した実数配列を生成する。
     * 開始値と終了値は配列に含まれる。
     *
     * @param start 開始値（正の有限値）
     * @param end   終了値（正の有限値）
     * @param count 要素数（1以上）
     * @return 生成された配列
     * @throws IllegalArgumentException 開始値または終了値が正の有限値でない場合、
     *                                  または要素数が1未満の場合
     */
    public static double[] geomspace(double start, double end, int count) {
        if (!Double.isFinite(start) || start <= 0.0) {
            throw new IllegalArgumentException("start must be a positive finite value");
        }
        if (!Double.isFinite(end) || end <= 0.0) {
            throw new IllegalArgumentException("end must be a positive finite value");
        }
        if (count < 1) {
            throw new IllegalArgumentException("count must be at least 1");
        }
        if (count == 1) {
            return new double[] { start };
        }

        double[] array = new double[count];
        double logStart = Math.log(start);
        double logStep = (Math.log(end) - logStart) / (count - 1);
        for (int i = 0; i < count; i++) {
            array[i] = Math.exp(logStart + i * logStep);
        }
        // 浮動小数点の丸めにかかわらず端点を厳密に保つ。
        array[0] = start;
        array[count - 1] = end;
        return array;
    }
}
