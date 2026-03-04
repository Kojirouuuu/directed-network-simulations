"""
有向ネットワーク上の EBCM（directed-infty.c と同等の計算）。

入次数 k_i のみの次数分布 P_i(k_i) を仮定。
- 分布タイプ: "Poi" (Poisson), "Pow" (Power law)
- Binom[k, m] = C(k, m) のテーブル
- EBCMConfig / DegreeDist の構築
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal, Optional

import numpy as np


# -----------------------------------------------------------------------------
# Binom テーブル
# -----------------------------------------------------------------------------

def build_binom(kmax: int) -> np.ndarray:
    """
    二項係数 C(k, m) のテーブルを構築する。
    形状: (kmax+1, kmax+1)。Binom[k, m] = C(k, m), ただし m > k の要素は未使用(0)。

    C の directed-infty.c の build_binom と同一。
    """
    n = kmax + 1
    Binom = np.zeros((n, n), dtype=float)
    for k in range(n):
        Binom[k, 0] = 1.0
        for m in range(1, n):
            Binom[k, m] = Binom[k, m - 1] * (k - m + 1) / m
    return Binom


# -----------------------------------------------------------------------------
# 次数分布 PMF
# -----------------------------------------------------------------------------

def build_poisson_pmf(avg: float, kmax: int) -> np.ndarray:
    """
    Poisson(avg) の PMF P[k] for k=0..kmax。正規化済み。
    """
    P = np.zeros(kmax + 1, dtype=float)
    P[0] = np.exp(-avg)
    for k in range(1, kmax + 1):
        P[k] = P[k - 1] * (avg / k)
    P /= P.sum()
    return P


def build_powerlaw_pmf(gamma: float, kmin: int, kmax: int) -> np.ndarray:
    """
    Power law P(k) ∝ k^{-gamma} for kmin <= k <= kmax。正規化済み。
    k < kmin の要素は 0。
    """
    P = np.zeros(kmax + 1, dtype=float)
    k_vals = np.arange(kmin, kmax + 1, dtype=float)
    Z = np.sum(k_vals ** (-gamma))
    P[kmin : kmax + 1] = k_vals ** (-gamma) / Z
    return P


def build_pmf(
    avg: float,
    gamma: float,
    kmin: int,
    kmax: int,
    dist_type: Literal["Poi", "Pow"],
) -> np.ndarray:
    """
    分布タイプに応じた PMF P[k], k=0..kmax を返す。
    """
    if dist_type == "Poi":
        return build_poisson_pmf(avg, kmax)
    if dist_type == "Pow":
        if kmin < 0:
            raise ValueError("Power law では kmin >= 0 が必要です。")
        return build_powerlaw_pmf(gamma, kmin, kmax)
    raise ValueError(f"未対応の分布タイプ: {dist_type}")


def get_effective_degree_max(
    dist_type: Literal["Poi", "Pow"],
    mean: float,
    kmax: int,
) -> int:
    """
    分布タイプに応じた実効的な最大次数（ループ範囲に使う）。
    Poi: mean + 3*sqrt(mean), Pow: kmax。
    """
    if dist_type == "Poi":
        return int(mean + 3.0 * np.sqrt(mean))
    return kmax


def calculate_mean_powerlaw(kmin: int, kmax: int, gamma: float) -> float:
    """Power law の平均次数 <k> = Σ k * k^{-gamma} / Z。"""
    k_vals = np.arange(kmin, kmax + 1, dtype=float)
    w = k_vals ** (-gamma)
    Z = w.sum()
    return float(np.sum(k_vals * w) / Z)


def create_degree_dist(
    mean: float,
    kmin: int,
    kmax: int,
    gamma: float,
    dist_type: Literal["Poi", "Pow"],
) -> np.ndarray:
    """
    次数分布 P[k], k=0..kmax を返す。長さは kmax+1。
    """
    return build_pmf(mean, gamma, kmin, kmax, dist_type)


# -----------------------------------------------------------------------------
# 設定と次数分布オブジェクト
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class DegreeConfig:
    """各次数分布（ここでは入次数のみ）の設定。"""
    mean: float
    min: int
    max: int
    gamma: float
    type: Literal["Poi", "Pow"]  # "Poi" or "Pow"


@dataclass(frozen=True)
class EBCMConfig:
    """EBCM 全体の次数分布設定（有向のみ: ki のみ）。"""
    ki: DegreeConfig


@dataclass
class DegreeDist:
    """
    構築済みの入次数分布。
    Pi[k] = P(k_i = k), k = 0..ki_max（配列長は ki_max+1）。
    """
    ki_min: int
    ki_max: int
    Pi: np.ndarray
    mean_ki: float
    type_i: Literal["Poi", "Pow"]
    gamma_i: float

    def __post_init__(self) -> None:
        if self.Pi.shape[0] != self.ki_max + 1:
            raise ValueError(
                f"Pi の長さは ki_max+1={self.ki_max + 1} である必要があります（実際: {self.Pi.shape[0]}）。"
            )


def build_degree_dist(cfg: EBCMConfig) -> DegreeDist:
    """
    EBCMConfig から DegreeDist を構築する。
    C の build_all_degree_dist と同等。
    """
    ki = cfg.ki
    type_i = ki.type
    gamma_i = ki.gamma

    ki_min = ki.min if type_i == "Pow" else 0
    ki_max = get_effective_degree_max(type_i, ki.mean, ki.max)

    if type_i == "Pow":
        mean_ki = calculate_mean_powerlaw(ki_min, ki_max, gamma_i)
    else:
        mean_ki = ki.mean

    Pi = create_degree_dist(mean_ki, ki_min, ki_max, gamma_i, type_i)

    return DegreeDist(
        ki_min=ki_min,
        ki_max=ki_max,
        Pi=Pi,
        mean_ki=mean_ki,
        type_i=type_i,
        gamma_i=gamma_i,
    )


def build_ebcm_config(
    dist_type: Literal["Poi", "Pow"],
    mean: float = 12.5,
    kmin: int = 0,
    kmax: int = 1000,
    gamma: float = 2.5,
) -> EBCMConfig:
    """
    よく使う EBCMConfig を簡単に作るヘルパー。
    """
    return EBCMConfig(
        ki=DegreeConfig(mean=mean, min=kmin, max=kmax, gamma=gamma, type=dist_type)
    )


# -----------------------------------------------------------------------------
# ダイナミクスパラメータ
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class DynamicsConfig:
    """ダイナミクスパラメータ（lambda_d, mu, rho0, T）。"""
    T: int
    rho0: float
    lambda_d: float
    mu: float


# -----------------------------------------------------------------------------
# Theta_d, xi_S, Phi の計算（C の Theta_d / xiS_directed / rhs_Phi と同等）
# -----------------------------------------------------------------------------

def theta_d(
    Binom: np.ndarray,
    ki: int,
    kmax: int,
    T: int,
    theta_d_val: float,
) -> float:
    """
    入次数 ki のノードが、閾値 T 未満のメッセージしか受け取らない確率。
    Θ_d(ki, θ_d) = Σ_{m=0}^{min(ki, T-1)} C(ki, m) θ_d^{ki-m} (1-θ_d)^m。
    theta_d=1 のとき 1、T=1 のとき θ_d^ki。
    """
    if theta_d_val == 1.0:
        return 1.0
    if T == 1:
        return float(np.power(theta_d_val, ki))
    if ki == 0:
        return 1.0
    md_max = min(ki, T - 1)
    q = 1.0 - theta_d_val
    s = 0.0
    for md in range(md_max + 1):
        a = ki - md
        s += Binom[ki, md] * (theta_d_val ** a) * (q ** md)
    return float(s)


def theta_d_prime(
    Binom: np.ndarray,
    ki: int,
    kmax: int,
    T: int,
    theta_d_val: float,
) -> float:
    """d/dθ Θ_d(ki, θ)。"""
    if ki == 0:
        return 0.0
    if T == 1:
        return float(ki * np.power(theta_d_val, ki - 1))
    md_max = min(ki, T - 1)
    q = 1.0 - theta_d_val
    s = 0.0
    for md in range(md_max + 1):
        a = ki - md
        b = md
        term = 0.0
        if a >= 1:
            term += a * (theta_d_val ** (a - 1)) * (q ** b)
        if b >= 1:
            term -= b * (theta_d_val ** a) * (q ** (b - 1))
        s += Binom[ki, md] * term
    return float(s)


def theta_d_prime_prime(
    Binom: np.ndarray,
    ki: int,
    kmax: int,
    T: int,
    theta_d_val: float,
) -> float:
    """d²/dθ² Θ_d(ki, θ)。T=1 のときは未使用（呼ばない想定）。"""
    if T == 1:
        raise ValueError("Theta_d_prime_prime は T=1 では使用しません。")
    if ki == 0:
        return 0.0
    md_max = min(ki, T - 1)
    q = 1.0 - theta_d_val
    s = 0.0
    for md in range(md_max + 1):
        a = ki - md
        b = md
        term = 0.0
        if a >= 2:
            term += a * (a - 1) * (theta_d_val ** (a - 2)) * (q ** b)
        if a >= 1 and b >= 1:
            term -= 2.0 * a * b * (theta_d_val ** (a - 1)) * (q ** (b - 1))
        if b >= 2:
            term += b * (b - 1) * (theta_d_val ** a) * (q ** (b - 2))
        s += Binom[ki, md] * term
    return float(s)


def xi_S(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """
    感受性ノードに占める「θ_d で重みづけた」割合（有向のみ）。
    ξ_S = (1 - ρ0) Σ_k P(k) Θ_d(k, θ_d)。S(t) と一致。
    """
    s = 0.0
    for k in range(D.ki_min, D.ki_max + 1):
        pk = D.Pi[k]
        if pk == 0.0:
            continue
        s += pk * theta_d(Binom, k, D.ki_max, dynamics.T, theta_d_val)
    return (1.0 - dynamics.rho0) * float(s)


def xi_S_prime(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """d/dθ ξ_S。"""
    s = 0.0
    for k in range(D.ki_min, D.ki_max + 1):
        pk = D.Pi[k]
        if pk == 0.0:
            continue
        s += pk * theta_d_prime(Binom, k, D.ki_max, dynamics.T, theta_d_val)
    return (1.0 - dynamics.rho0) * float(s)


def xi_S_prime_prime(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """d²/dθ² ξ_S。"""
    if dynamics.T == 1:
        raise ValueError("xi_S_prime_prime は T=1 では使用しません。")
    s = 0.0
    for k in range(D.ki_min, D.ki_max + 1):
        pk = D.Pi[k]
        if pk == 0.0:
            continue
        s += pk * theta_d_prime_prime(Binom, k, D.ki_max, dynamics.T, theta_d_val)
    return (1.0 - dynamics.rho0) * float(s)


def compute_S(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """感受性ノードの割合 S(t) = ξ_S（有向のみでは同一）。"""
    return xi_S(D, dynamics, Binom, theta_d_val)


def phi(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """
    閾値ちょうどで「アクティブ」になる項（R の式で用いる Φ）。
    Φ = Σ_k P(k) C(k, T-1) θ_d^{k-(T-1)} (1-θ_d)^{T-1}。
    theta_d=1 かつ T=1 のとき 1、theta_d=1 かつ T>1 のとき 0。
    """
    T = dynamics.T
    if theta_d_val == 1.0:
        return 1.0 if T == 1 else 0.0
    s = 0.0
    for k in range(D.ki_min, D.ki_max + 1):
        pk = D.Pi[k]
        if pk == 0.0:
            continue
        if T == 1:
            s += pk * (theta_d_val ** k)
            continue
        if k == 0:
            continue
        q = 1.0 - theta_d_val
        term = (
            Binom[k, T - 1]
            * (theta_d_val ** (k - (T - 1)))
            * (q ** (T - 1))
        )
        s += pk * term
    return float(s)


def g_d(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """dθ_d/dt = g_d = -λ_d (θ_d - ξ_S) + μ(1 - θ_d)。"""
    xi = xi_S(D, dynamics, Binom, theta_d_val)
    return (
        -dynamics.lambda_d * (theta_d_val - xi)
        + dynamics.mu * (1.0 - theta_d_val)
    )


def g_d_prime(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """d/dθ g_d。根探索・安定性解析用。"""
    xi_p = xi_S_prime(D, dynamics, Binom, theta_d_val)
    return -dynamics.lambda_d * (1.0 - xi_p) - dynamics.mu


def g_d_prime_prime(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    theta_d_val: float,
) -> float:
    """d²/dθ² g_d。"""
    xi_pp = xi_S_prime_prime(D, dynamics, Binom, theta_d_val)
    return dynamics.lambda_d * xi_pp


# -----------------------------------------------------------------------------
# ルンゲクッタ法による S, A, R の時間発展（ebcm.py の lambda_u=0 相当）
# -----------------------------------------------------------------------------

def rhs_directed(
    t: float,
    y: np.ndarray,
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
) -> np.ndarray:
    """
    y = [θ_d, R] に対する ODE 右辺。
    dθ_d/dt = g_d,  dR/dt = μ A,  A = 1 - S - R。
    S は θ_d から compute_S で求める（陽に dS/dt は積分しない）。
    """
    theta_d_val = float(y[0])
    R = float(y[1])
    dtheta_d = g_d(D, dynamics, Binom, theta_d_val)
    S = compute_S(D, dynamics, Binom, theta_d_val)
    A = max(0.0, 1.0 - S - R)
    dR = dynamics.mu * A
    return np.array([dtheta_d, dR], dtype=float)


def rk4_step_vec(
    f: Callable[[float, np.ndarray], np.ndarray],
    t: float,
    y: np.ndarray,
    dt: float,
) -> np.ndarray:
    """1 ステップ RK4。f(t, y) は dy/dt を返す。"""
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_directed_ebcm(
    D: DegreeDist,
    dynamics: DynamicsConfig,
    Binom: np.ndarray,
    t_max: float,
    dt: float,
    y0: Optional[np.ndarray] = None,
) -> dict[str, np.ndarray]:
    """
    ルンゲクッタ 4 次法で θ_d, R を t=0..t_max まで積分し、
    S, A, R の時系列を返す。ebcm.py の simulate_ebcm と同様（有向のみで y=[θ_d, R]）。

    返り値:
      t, theta_d, S, A, R（各 1 次元配列、長さ n_steps）
    """
    n_steps = int(np.ceil(t_max / dt)) + 1
    t = np.linspace(0.0, t_max, n_steps)

    if y0 is None:
        # θ_d(0)=1, R(0)=0（本文の設定と一致）
        y = np.array([1.0, 0.0], dtype=float)
    else:
        y = np.array(y0, dtype=float).copy()

    theta_d_series = np.empty(n_steps, dtype=float)
    S_series = np.empty(n_steps, dtype=float)
    A_series = np.empty(n_steps, dtype=float)
    R_series = np.empty(n_steps, dtype=float)

    def rhs(t_: float, y_: np.ndarray) -> np.ndarray:
        return rhs_directed(t_, y_, D, dynamics, Binom)

    for i in range(n_steps):
        theta_d_val = float(y[0])
        R = float(y[1])
        S = compute_S(D, dynamics, Binom, theta_d_val)
        A = 1.0 - S - R

        theta_d_series[i] = theta_d_val
        S_series[i] = S
        A_series[i] = A
        R_series[i] = R

        if i < n_steps - 1:
            y = rk4_step_vec(rhs, t[i], y, dt)
            y[0] = float(np.clip(y[0], 0.0, 1.0))
            y[1] = float(np.clip(y[1], 0.0, 1.0))

    return {
        "t": t,
        "theta_d": theta_d_series,
        "S": S_series,
        "A": A_series,
        "R": R_series,
    }
