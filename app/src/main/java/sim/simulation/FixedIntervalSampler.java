package sim.simulation;

import java.util.function.DoubleConsumer;

/** 固定時刻格子でシミュレーション状態を記録するための内部ヘルパー。 */
final class FixedIntervalSampler {

    private static final FixedIntervalSampler DISABLED = new FixedIntervalSampler();

    private final boolean enabled;
    private final double dt;
    private final long lastSampleIndex;
    private long nextSampleIndex;

    private FixedIntervalSampler() {
        this.enabled = false;
        this.dt = Double.NaN;
        this.lastSampleIndex = 0L;
        this.nextSampleIndex = 1L;
    }

    private FixedIntervalSampler(double dt, long lastSampleIndex) {
        this.enabled = true;
        this.dt = dt;
        this.lastSampleIndex = lastSampleIndex;
        this.nextSampleIndex = 1L;
    }

    static FixedIntervalSampler disabled() {
        return DISABLED;
    }

    static FixedIntervalSampler create(double dt, double tMax) {
        if (!Double.isFinite(dt) || dt <= 0.0) {
            throw new IllegalArgumentException("dt must be finite and positive");
        }
        if (!Double.isFinite(tMax) || tMax <= 0.0) {
            throw new IllegalArgumentException("tMax must be finite and positive for fixed-interval recording");
        }

        double quotient = tMax / dt;
        double adjustedQuotient = quotient + 8.0 * Math.ulp(quotient);
        double lastIndexAsDouble = Math.floor(adjustedQuotient);
        if (!Double.isFinite(lastIndexAsDouble) || lastIndexAsDouble > Integer.MAX_VALUE - 1L) {
            throw new IllegalArgumentException("dt produces too many time-series rows");
        }
        return new FixedIntervalSampler(dt, (long) lastIndexAsDouble);
    }

    boolean isEnabled() {
        return enabled;
    }

    /** eventTime より前の未記録格子点を、イベント処理前の状態で記録する。 */
    void recordBefore(double eventTime, DoubleConsumer recorder) {
        if (!enabled) {
            return;
        }
        while (hasNext() && sampleTime() < eventTime) {
            recordNext(recorder);
        }
    }

    /** eventTime 以下の未記録格子点を、イベント処理後の状態で記録する。 */
    void recordThrough(double eventTime, DoubleConsumer recorder) {
        if (!enabled) {
            return;
        }
        while (hasNext() && sampleTime() <= eventTime) {
            recordNext(recorder);
        }
    }

    /** シミュレーション終了後、最後の格子点まで最終状態を記録する。 */
    void recordRemaining(DoubleConsumer recorder) {
        if (!enabled) {
            return;
        }
        while (hasNext()) {
            recordNext(recorder);
        }
    }

    private boolean hasNext() {
        return nextSampleIndex <= lastSampleIndex;
    }

    private double sampleTime() {
        return nextSampleIndex * dt;
    }

    private void recordNext(DoubleConsumer recorder) {
        recorder.accept(sampleTime());
        nextSampleIndex++;
    }
}
