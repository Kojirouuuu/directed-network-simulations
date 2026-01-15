package sim.utils;

import java.nio.file.Files;
import java.nio.file.Path;

/**
 * パス操作のユーティリティクラス。
 */
public final class PathsEx {
    private PathsEx() {}

    /**
     * 指定されたパスが存在する場合、存在しない新しいパスを返す。
     * 新しいパスは "name (1).ext", "name (2).ext", ... のようにインデックス付きのサフィックスを持つ。
     * パスが存在しない場合は、そのまま返す。
     *
     * @param path チェックするパス
     * @return 存在しないパス（元のパスが存在しない場合はそのまま、存在する場合はインデックス付き）
     */
    public static Path resolveIndexed(Path path) {
        if (path == null) {
            return null;
        }
        if (!Files.exists(path)) {
            return path;
        }

        Path dir = path.getParent();
        String filename = path.getFileName().toString();

        int dot = filename.lastIndexOf('.');
        String base;
        String ext;
        // 位置0のドットは拡張子なしの隠しファイルを意味する
        if (dot > 0) {
            base = filename.substring(0, dot);
            ext = filename.substring(dot); // ドットを含む
        } else {
            base = filename;
            ext = "";
        }

        int i = 1;
        while (true) {
            String candidate = String.format("%s (%d)%s", base, i, ext);
            Path p = dir == null ? Path.of(candidate) : dir.resolve(candidate);
            if (!Files.exists(p)) {
                return p;
            }
            i++;
        }
    }
}
