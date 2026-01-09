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
     * @param arr    シャッフルする配列
     * @param random 乱数ジェネレータ
     */
    public static void shuffle(int[] arr, Random random) {
        for (int i = arr.length - 1; i > 0; i--) {
            int j = random.nextInt(i + 1);
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
}

