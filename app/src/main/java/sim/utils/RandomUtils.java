package sim.utils;

import java.util.Random;

/**
 * 乱数生成のユーティリティクラス。
 */
public final class RandomUtils {
    private RandomUtils() {
    }

    /**
     * ポアソン分布に従う乱数を生成する。
     *
     * @param lambda ポアソン分布のパラメータ（平均）
     * @param random 乱数ジェネレータ
     * @return ポアソン分布に従う乱数
     */
    public static int getPoisson(double lambda, Random random) {
        double L = Math.exp(-lambda);
        double p = 1.0;
        int k = 0;

        do {
            k++;
            p *= random.nextDouble();
        } while (p > L);

        return k - 1;
    }
}
