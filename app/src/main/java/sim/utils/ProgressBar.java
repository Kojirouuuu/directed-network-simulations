package sim.utils;

/** Console-only formatting; callers retain ownership of progress and thread lifecycle. */
public final class ProgressBar {
    private ProgressBar() {}

    public static void render(long done, long total) {
        int filled = total == 0 ? 100 : (int) Math.min(100, done * 100 / total);
        int percent = total == 0 ? 100 : (int) Math.min(100, done * 100 / total);
        String bar = "#".repeat(filled) + "-".repeat(100 - filled);
        System.out.print("\033[2K\rProgress [%s] %3d%% (%d/%d)".formatted(bar, percent, done, total));
        System.out.flush();
    }
}
