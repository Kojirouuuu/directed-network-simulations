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
}

