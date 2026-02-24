import numpy as np
from numpy.typing import NDArray

def sample_scale_free_degrees(
    n_samples: int,
    k_min: int,
    k_max: int,
    gamma: float,
    seed: int | None = None,
) -> tuple[float, NDArray[np.int32]]:
    """
    離散パワーロー P(k) ∝ k^{-gamma}, k ∈ [k_min, k_max] から n_samples 個サンプリングし、
    - 各次数 k の出現回数 counts（長さ: k_max-k_min+1）
    - サンプル平均次数 mean_degree
    - サンプルした次数列 degrees（長さ: n_samples）
    を返す。

    Returns
    -------
    counts : (k_max-k_min+1,) int64
        counts[i] は k = (k_min + i) の出現回数
    mean_degree : float
        degrees の平均
    degrees : (n_samples,) int32
        サンプルした次数列
    """
    if n_samples <= 0:
        raise ValueError("n_samples must be positive.")
    if not (1 <= k_min <= k_max):
        raise ValueError("Require 1 <= k_min <= k_max.")
    if gamma <= 0:
        raise ValueError("gamma must be > 0.")

    rng = np.random.default_rng(seed)

    # --- サポート (k_min..k_max) と未正規化重み w(k)=k^{-gamma} ---
    k_support = np.arange(k_min, k_max + 1, dtype=np.float64)
    unnormalized = np.power(k_support, -gamma)

    # 数値安定性：総和が0にならないはずだが念のため
    weight_sum = float(unnormalized.sum())
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        raise ValueError("Invalid weights; check k_min/k_max/gamma.")

    # --- 離散分布のCDF（累積分布関数） ---
    cdf = np.cumsum(unnormalized) / weight_sum
    cdf[-1] = 1.0  # 浮動小数点誤差対策（searchsortedのはみ出し防止）

    # --- 逆変換サンプリング ---
    u = rng.random(n_samples)
    indices = np.searchsorted(cdf, u, side="left")  # 0..len(k_support)-1
    degrees = (k_min + indices).astype(np.int32)

    # --- サンプル統計 ---
    mean_degree = float(degrees.mean())
    # k_min..k_max に揃えたcounts（後処理が楽）
    counts = np.bincount(degrees - k_min, minlength=(k_max - k_min + 1))

    return mean_degree, counts


def theoretical_scale_free_distribution(
    k_min: int,
    k_max: int,
    gamma: float,
) -> tuple[float, NDArray[np.int32]]:
    """
    離散パワーロー P(k) ∝ k^{-gamma}, k ∈ [k_min, k_max] の
    - 正規化確率分布 P(k)
    - 理論平均次数 E[K]
    - 対応するkサポート
    を返す。

    Returns
    -------
    pk : (k_max-k_min+1,) float64
        pk[i] = P(K = k_min + i)
    mean_k : float
        E[K] = Σ_k k P(k)
    k_values : (k_max-k_min+1,) int32
        k_min..k_max
    """
    if not (1 <= k_min <= k_max):
        raise ValueError("Require 1 <= k_min <= k_max.")
    if gamma <= 0:
        raise ValueError("gamma must be > 0.")

    k_values = np.arange(k_min, k_max + 1, dtype=np.float64)

    # 未正規化 w(k)=k^{-gamma}
    w = np.power(k_values, -gamma)
    Z = float(w.sum())
    if not np.isfinite(Z) or Z <= 0:
        raise ValueError("Invalid normalization constant; check parameters.")

    pk = w / Z

    # E[K] = Σ k * P(k) = (Σ k^{1-gamma}) / (Σ k^{-gamma})
    mean_k = float((k_values * pk).sum())

    return mean_k, pk.astype(np.float64)
