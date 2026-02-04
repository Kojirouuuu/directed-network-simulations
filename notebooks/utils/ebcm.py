"""
半有向（有向＋無向）ネットワーク上の SAR の EBCM。
非冗長メモリと閾値採用を仮定。

モデル（閉じた ODE）:
  dθ_d/dt = -λ_d [θ_d - ξ_S^(d)(θ_d,θ_u)] + μ(1-θ_d)
  dθ_u/dt = -λ_u [θ_u - ξ_S^(u)(θ_d,θ_u)] + μ(1-θ_u)

ただし
  ξ_S^(d) = (1-ρ0) Σ_{k_i,k_o,k_u} (k_o P / <k_o>) Θ( k_i, k_u ; θ_d, θ_u )
  ξ_S^(u) = (1-ρ0) Σ_{k_i,k_o,k_u} (k_u P / <k_u>) Θ( k_i, k_u-1 ; θ_d, θ_u )

Θ は有向入辺 (k_i) と無向辺 (k_u) 上の二項受信で
総受信メッセージ m_d+m_u < T となる確率。

本コードの対応:
  - 因数分解された次数分布 P(k_i,k_o,k_u)=P_i P_o P_u（デフォルト・本文と一致）
  - （任意）疎な同時 PMF を (k_i,k_o,k_u,p) のリストで指定可能

R(t) の時間発展:
  dR/dt = μ A(t),  A(t)=1-S(t)-R(t)
dS/dt は陽には計算しない。
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Tuple, Union, Optional, Sequence

import numpy as np
from scipy.stats import poisson
from scipy.special import gammaln


# -----------------------------
# 次数分布（1次元）
# -----------------------------

class DistributionType(Enum):
    poisson = "poisson"
    power_law = "power_law"
    constant = "constant"


@dataclass(frozen=True)
class DegreeDistribution:
    distribution_type: DistributionType
    N: int
    distribution: np.ndarray          # k=0..N-1 に対する P(k)
    valid_deg_range: Tuple[int, int]  # P(k)>0 となる (k_min, k_max)

    @property
    def k_values(self) -> np.ndarray:
        return np.arange(0, self.N)

    @property
    def pk(self) -> np.ndarray:
        return self.distribution

    @property
    def k_min(self) -> int:
        return self.valid_deg_range[0]

    @property
    def k_max(self) -> int:
        return self.valid_deg_range[1]

    def mean_k(self) -> float:
        k = self.k_values
        return float(np.sum(k * self.pk))


def degree_distribution(
    distribution_type: Union[DistributionType, str],
    N: int,
    k_ave: int,
    gamma: float = 3.0,
    k_min: int = 1,
) -> DegreeDistribution:
    """
    正規化された離散次数分布 P(k), k=0..N-1 を構築する。
    Poisson: Poisson(k_ave)、3*k_ave で打ち切り。
    Power Law: P(k) ∝ k^{-gamma}
    Constant: δ(k - k_ave)
    """
    if isinstance(distribution_type, str):
        distribution_type = DistributionType(distribution_type)

    k = np.arange(0, N)

    if distribution_type == DistributionType.poisson:
        pk = poisson.pmf(k, k_ave)
        pk[k > 3 * k_ave] = 0.0

    elif distribution_type == DistributionType.power_law:
        pk = np.zeros_like(k, dtype=float)
        ks = np.arange(k_min, N)
        pk[ks] = ks ** (-gamma)

    elif distribution_type == DistributionType.constant:
        pk = np.zeros_like(k, dtype=float)
        if not (0 <= k_ave < N):
            raise ValueError("RR では k_ave は [0, N-1] の範囲である必要があります。")
        pk[k_ave] = 1.0

    else:
        raise ValueError(f"未知の次数分布タイプ: {distribution_type}")

    s = pk.sum()
    if s <= 0:
        raise ValueError("次数分布の和が 0 です。")
    pk = pk / s

    nonzero = np.where(pk > 0)[0]
    k_min_eff = int(nonzero.min()) if nonzero.size else 0
    k_max_eff = int(nonzero.max()) if nonzero.size else 0

    return DegreeDistribution(
        distribution_type=distribution_type,
        N=N,
        distribution=pk,
        valid_deg_range=(k_min_eff, k_max_eff),
    )


# -----------------------------
# 半有向次数モデル
# -----------------------------

@dataclass(frozen=True)
class SemiDirectedDegreeFactorized:
    """
    因数分解された次数分布:
      P(k_i,k_o,k_u) = P_i(k_i) P_o(k_o) P_u(k_u)

    本文の仮定と一致。
    """
    pk_in: DegreeDistribution
    pk_out: DegreeDistribution
    pk_und: DegreeDistribution

    def mean_ki(self) -> float:
        return self.pk_in.mean_k()

    def mean_ko(self) -> float:
        return self.pk_out.mean_k()

    def mean_ku(self) -> float:
        return self.pk_und.mean_k()


@dataclass(frozen=True)
class SparseJointDegreePMF:
    """
    任意: (k_i,k_o,k_u,p) で与える疎な同時 PMF。
    p の和は 1 であること。
    """
    ki: np.ndarray
    ko: np.ndarray
    ku: np.ndarray
    p: np.ndarray

    def __post_init__(self) -> None:
        if not (self.ki.shape == self.ko.shape == self.ku.shape == self.p.shape):
            raise ValueError("ki, ko, ku, p は同じ shape である必要があります。")
        s = float(np.sum(self.p))
        if not np.isfinite(s) or abs(s - 1.0) > 1e-8:
            raise ValueError(f"同時 PMF の和は 1 である必要があります（得られた値: {s}）。")

    def mean_ko(self) -> float:
        return float(np.sum(self.ko * self.p))

    def mean_ku(self) -> float:
        return float(np.sum(self.ku * self.p))


@dataclass(frozen=True)
class EBCMSemiDirectedConfig:
    """
    半有向 SAR EBCM のパラメータ。
    """
    rho0: float
    mu: float
    T: int
    lambda_d: float
    lambda_u: float

    # 以下のいずれか一方を指定すること:
    deg_factorized: Optional[SemiDirectedDegreeFactorized] = None
    deg_joint_sparse: Optional[SparseJointDegreePMF] = None

    def __post_init__(self) -> None:
        if (self.deg_factorized is None) == (self.deg_joint_sparse is None):
            raise ValueError("deg_factorized と deg_joint_sparse のどちらか一方だけを指定してください。")
        if not (0.0 <= self.rho0 <= 1.0):
            raise ValueError("rho0 は [0,1] の範囲である必要があります。")
        if self.mu < 0:
            raise ValueError("mu は >= 0 である必要があります。")
        if self.T <= 0:
            raise ValueError("T は >= 1 である必要があります。")
        if self.lambda_d < 0 or self.lambda_u < 0:
            raise ValueError("lambda_d と lambda_u は >= 0 である必要があります。")


# -----------------------------
# 二項分布用ヘルパー（m_max = T-1 で打ち切り）
# -----------------------------

def _binom_pmf_trunc_table(n: np.ndarray, theta: float, m_max: int) -> np.ndarray:
    """
    m=0..m_max に対する PMF テーブルを返す:
      P(M=m) = C(n,m) (1-theta)^m theta^(n-m)
    M ~ Bin(n, p=1-theta)。

    形状: (len(n), m_max+1)
    """
    n = np.asarray(n, dtype=int)
    n = np.maximum(n, 0)
    m = np.arange(m_max + 1, dtype=int)[None, :]      # (1, M)
    nn = n[:, None]                                   # (N, 1)

    out = np.zeros((n.size, m_max + 1), dtype=float)
    valid = m <= nn
    if not np.any(valid):
        return out

    # 数値安定のため theta=0 または 1 の端は明示的に扱う。
    if theta <= 0.0:
        # theta=0 => p=1 なので M=n に確定。
        for i, ni in enumerate(n):
            if 0 <= ni <= m_max:
                out[i, ni] = 1.0
        return out

    if theta >= 1.0:
        # theta=1 => p=0 なので M=0 に確定。
        out[:, 0] = 1.0
        return out

    log_theta = np.log(theta)
    log_1m = np.log1p(-theta)

    # log C(n,m) + m log(1-theta) + (n-m) log(theta)
    # 安定化のため gammaln を使用
    logC = gammaln(nn + 1.0) - gammaln(m + 1.0) - gammaln(nn - m + 1.0)
    logpmf = logC + m * log_1m + (nn - m) * log_theta

    # 無効な要素は -inf にし exp=0 とする
    logpmf = np.where(valid, logpmf, -np.inf)
    out = np.exp(logpmf)
    return out


def _binom_cdf_trunc_table(n: np.ndarray, theta: float, m_max: int) -> np.ndarray:
    """
    r=0..m_max に対する CDF テーブル:
      CDF[r] = P(M <= r)
    打ち切り PMF を使用。
    形状: (len(n), m_max+1)
    """
    pmf = _binom_pmf_trunc_table(n, theta, m_max)
    return np.cumsum(pmf, axis=1)


# -----------------------------
# EBCM 本体の計算
# -----------------------------

@dataclass
class _PrecomputedWeights:
    """
    (k_i, k_u) の台上の行列:
      p_ik_u:  sum_{k_o} P(k_i,k_o,k_u)                      (S(t) 用)
      w_d:     sum_{k_o} [k_o P(k_i,k_o,k_u)] / <k_o>        (ξ_S^(d) 用)
      w_u:     sum_{k_o} [k_u P(k_i,k_o,k_u)] / <k_u>        (ξ_S^(u) 用)
    """
    ki_vals: np.ndarray
    ku_vals: np.ndarray
    p_ik_u: np.ndarray
    w_d: np.ndarray
    w_u: np.ndarray


def _build_weights_factorized(deg: SemiDirectedDegreeFactorized) -> _PrecomputedWeights:
    ki_min, ki_max = deg.pk_in.valid_deg_range
    ku_min, ku_max = deg.pk_und.valid_deg_range

    ki_vals = np.arange(ki_min, ki_max + 1, dtype=int)
    ku_vals = np.arange(ku_min, ku_max + 1, dtype=int)

    Pi = deg.pk_in.pk[ki_min:ki_max + 1]
    Pu = deg.pk_und.pk[ku_min:ku_max + 1]

    p_ik_u = Pi[:, None] * Pu[None, :]

    # 因数分解 P_i P_o P_u の場合:
    # ξ_S^(d) の重みは Σ k_o P_o / <k_o> = 1 より p_ik_u に簡約される。
    w_d = p_ik_u.copy()

    ku_mean = deg.mean_ku()
    if ku_mean > 0:
        w_u = Pi[:, None] * ((ku_vals * Pu) / ku_mean)[None, :]
    else:
        # 実質的に無向辺なし
        w_u = np.zeros_like(p_ik_u)

    # 整合性: (ki,ku) 上で p_ik_u の和が 1 になるように
    p_sum = float(np.sum(p_ik_u))
    if abs(p_sum - 1.0) > 1e-8:
        p_ik_u = p_ik_u / p_sum
        w_d = w_d / float(np.sum(w_d))
        if ku_mean > 0:
            w_u = w_u / float(np.sum(w_u))

    return _PrecomputedWeights(ki_vals=ki_vals, ku_vals=ku_vals, p_ik_u=p_ik_u, w_d=w_d, w_u=w_u)


def _build_weights_sparse_joint(joint: SparseJointDegreePMF) -> _PrecomputedWeights:
    # k_i と k_u の台を構築
    ki_vals = np.unique(joint.ki.astype(int))
    ku_vals = np.unique(joint.ku.astype(int))
    ki_vals.sort()
    ku_vals.sort()

    ki_index = {k: i for i, k in enumerate(ki_vals)}
    ku_index = {k: i for i, k in enumerate(ku_vals)}

    p_ik_u = np.zeros((ki_vals.size, ku_vals.size), dtype=float)
    w_d = np.zeros_like(p_ik_u)
    w_u = np.zeros_like(p_ik_u)

    ko_mean = joint.mean_ko()
    ku_mean = joint.mean_ku()
    if ko_mean <= 0:
        raise ValueError("有向辺の重み付けには <k_o> が正である必要があります。")
    if ku_mean < 0:
        raise ValueError("<k_u> は非負である必要があります。")

    for ki, ko, ku, p in zip(joint.ki, joint.ko, joint.ku, joint.p):
        i = ki_index[int(ki)]
        j = ku_index[int(ku)]
        p_ik_u[i, j] += float(p)
        w_d[i, j] += float(ko) * float(p) / ko_mean
        if ku_mean > 0:
            w_u[i, j] += float(ku) * float(p) / ku_mean

    # 数値的な正規化（わずかなずれを吸収）
    p_ik_u /= float(np.sum(p_ik_u))
    w_d /= float(np.sum(w_d))
    if ku_mean > 0 and float(np.sum(w_u)) > 0:
        w_u /= float(np.sum(w_u))

    return _PrecomputedWeights(ki_vals=ki_vals, ku_vals=ku_vals, p_ik_u=p_ik_u, w_d=w_d, w_u=w_u)


@dataclass
class SemiDirectedSAREBCM:
    config: EBCMSemiDirectedConfig

    def __post_init__(self) -> None:
        if self.config.deg_factorized is not None:
            self._w = _build_weights_factorized(self.config.deg_factorized)
        else:
            assert self.config.deg_joint_sparse is not None
            self._w = _build_weights_sparse_joint(self.config.deg_joint_sparse)

        self._m_max = self.config.T - 1  # T-1 で打ち切り

    def _theta_matrices(self, theta_d: float, theta_u: float) -> tuple[np.ndarray, np.ndarray]:
        """
        (k_i,k_u) グリッド上で Θ 行列を計算:
          Theta_full[i,j]   = Θ(k_i, k_u;    theta_d, theta_u)
          Theta_cav_u[i,j]  = Θ(k_i, k_u-1;  theta_d, theta_u)   （無向キャビティ用）
        """
        ki = self._w.ki_vals
        ku = self._w.ku_vals

        pmf_d = _binom_pmf_trunc_table(ki, theta_d, self._m_max)          # (nki, T)
        cdf_u = _binom_cdf_trunc_table(ku, theta_u, self._m_max)          # (nku, T)
        cdf_u_rev = cdf_u[:, ::-1]                                       # 列は r=T-1..0
        Theta_full = pmf_d @ cdf_u_rev.T                                 # (nki, nku)

        ku_cav = np.maximum(ku - 1, 0)
        cdf_u_cav = _binom_cdf_trunc_table(ku_cav, theta_u, self._m_max)
        Theta_cav_u = pmf_d @ (cdf_u_cav[:, ::-1].T)

        return Theta_full, Theta_cav_u

    def S_xi(self, theta_d: float, theta_u: float) -> tuple[float, float, float]:
        """
        次を計算:
          S(t)          = (1-ρ0) Σ p_ik_u * Θ_full
          ξ_S^(d)(t)    = (1-ρ0) Σ w_d    * Θ_full
          ξ_S^(u)(t)    = (1-ρ0) Σ w_u    * Θ_cav_u
        """
        Theta_full, Theta_cav_u = self._theta_matrices(theta_d, theta_u)

        one_minus_rho0 = 1.0 - self.config.rho0
        S = one_minus_rho0 * float(np.sum(self._w.p_ik_u * Theta_full))
        xi_d = one_minus_rho0 * float(np.sum(self._w.w_d * Theta_full))
        xi_u = one_minus_rho0 * float(np.sum(self._w.w_u * Theta_cav_u))

        return S, xi_d, xi_u

    def rhs_theta(self, theta_d: float, theta_u: float) -> tuple[float, float, float]:
        """
        (θ_d, θ_u) の右辺。便宜上 S(t) も返す。
        """
        S, xi_d, xi_u = self.S_xi(theta_d, theta_u)

        dtheta_d = -self.config.lambda_d * (theta_d - xi_d) + self.config.mu * (1.0 - theta_d)
        dtheta_u = -self.config.lambda_u * (theta_u - xi_u) + self.config.mu * (1.0 - theta_u)

        return dtheta_d, dtheta_u, S

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        y = [θ_d, θ_u, R] に対する ODE 系。
        R は dR/dt = μ A（A = 1 - S - R）で時間発展。
        dS/dt を陽に計算しない。
        """
        theta_d, theta_u, R = float(y[0]), float(y[1]), float(y[2])
        dtheta_d, dtheta_u, S = self.rhs_theta(theta_d, theta_u)

        A = max(0.0, 1.0 - S - R)
        dR = self.config.mu * A

        return np.array([dtheta_d, dtheta_u, dR], dtype=float)

    def F_map(self, theta_d: float, theta_u: float) -> tuple[float, float]:
        """
        定常自己整合写像 F = (Fd, Fu) を返す。
          Fd = (λd ξ_S^(d) + μ) / (λd + μ)
          Fu = (λu ξ_S^(u) + μ) / (λu + μ)
        """
        _, xi_d, xi_u = self.S_xi(theta_d, theta_u)

        denom_d = self.config.lambda_d + self.config.mu
        denom_u = self.config.lambda_u + self.config.mu

        # λ+μ=0 の退化ケースは「方程式が拘束しない」ので NaN を返して明示
        Fd = ((self.config.lambda_d * xi_d + self.config.mu) / denom_d) if denom_d > 0 else np.nan
        Fu = ((self.config.lambda_u * xi_u + self.config.mu) / denom_u) if denom_u > 0 else np.nan

        return float(Fd), float(Fu)

    def fixed_point_residual(self, theta_d: float, theta_u: float) -> tuple[float, float]:
        """
        残差 g = (Fd-θd, Fu-θu) を返す。
        """
        Fd, Fu = self.F_map(theta_d, theta_u)
        return float(Fd - theta_d), float(Fu - theta_u)

    def fixed_point_residual_vec(self, x: np.ndarray) -> np.ndarray:
        """
        root/solve 用: x=[θd,θu] -> [Fd-θd, Fu-θu]
        """
        td, tu = float(x[0]), float(x[1])
        gd, gu = self.fixed_point_residual(td, tu)
        return np.array([gd, gu], dtype=float)



# -----------------------------
# 簡易積分器（RK4）
# -----------------------------

def rk4_step_vec(f, t: float, y: np.ndarray, dt: float) -> np.ndarray:
    k1 = f(t, y)
    k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def simulate_ebcm(
    model: SemiDirectedSAREBCM,
    t_max: float,
    dt: float,
    y0: Optional[np.ndarray] = None,
) -> dict[str, np.ndarray]:
    """
    y=[θ_d,θ_u,R] を t=0..t_max まで積分（RK4）。
    θ_d, θ_u, S, A, R の時系列を返す。
    """
    n_steps = int(np.ceil(t_max / dt)) + 1
    t = np.linspace(0.0, t_max, n_steps)

    if y0 is None:
        # θ_d(0)=θ_u(0)=1, R(0)=0（本文の設定と一致）
        y = np.array([1.0, 1.0, 0.0], dtype=float)
    else:
        y = np.array(y0, dtype=float).copy()

    theta_d_series = np.empty(n_steps, dtype=float)
    theta_u_series = np.empty(n_steps, dtype=float)
    S_series = np.empty(n_steps, dtype=float)
    A_series = np.empty(n_steps, dtype=float)
    R_series = np.empty(n_steps, dtype=float)

    for i in range(n_steps):
        theta_d, theta_u, R = float(y[0]), float(y[1]), float(y[2])
        S, _, _ = model.S_xi(theta_d, theta_u)
        A = 1.0 - S - R

        theta_d_series[i] = theta_d
        theta_u_series[i] = theta_u
        S_series[i] = S
        A_series[i] = A
        R_series[i] = R

        if i < n_steps - 1:
            y = rk4_step_vec(model.rhs, t[i], y, dt)
            # 有効範囲にクランプ（数値安定のため）
            y[0] = float(np.clip(y[0], 0.0, 1.0))
            y[1] = float(np.clip(y[1], 0.0, 1.0))
            y[2] = float(np.clip(y[2], 0.0, 1.0))

    return {
        "t": t,
        "theta_d": theta_d_series,
        "theta_u": theta_u_series,
        "S": S_series,
        "A": A_series,
        "R": R_series,
    }

def newton_solve_theta(
    model,
    x0,
    tol=1e-10,
    max_iter=50,
    eps=1e-6,
    damping=True,
    clamp01=True,
    verbose=True,
):
    """
    Newton 法で g(θ_d, θ_u)=0 を解く。
    g は model.fixed_point_residual_vec(x)、x=[θ_d, θ_u]。

    初期値 x0、終了条件 tol（||g||_2）、最大反復 max_iter。
    ヤコビアンは中心差分（eps）。damping=True でバックトラック直線探索、
    clamp01=True で θ_d, θ_u を [0,1] にクランプ。verbose で反復ログ出力。

    戻り値: 解 x (shape (2,)) と info 辞書（converged, iter, residual_norm 等）。
    """
    x = np.array(x0, dtype=float).reshape(2,)
    if clamp01:
        x = np.clip(x, 0.0, 1.0)

    def g(x):
        return model.fixed_point_residual_vec(x)

    def jacobian_central(x, eps):
        # J[:,j] = (g(x+eps e_j) - g(x-eps e_j)) / (2 eps)
        J = np.zeros((2, 2), dtype=float)
        for j in range(2):
            dx = np.zeros(2, dtype=float)
            dx[j] = eps
            xp = x + dx
            xm = x - dx
            if clamp01:
                xp = np.clip(xp, 0.0, 1.0)
                xm = np.clip(xm, 0.0, 1.0)
            gp = g(xp)
            gm = g(xm)
            J[:, j] = (gp - gm) / (2.0 * eps)
        return J

    gval = g(x)
    ng = float(np.linalg.norm(gval, ord=2))

    # if verbose:
    #     print(f"[Newton] iter=0 x={x} ||g||={ng:.3e} g={gval}")

    for it in range(1, max_iter + 1):
        if ng < tol:
            return x, {
                "converged": True,
                "iter": it - 1,
                "residual_norm": ng,
                "residual": gval,
            }

        J = jacobian_central(x, eps)

        # Solve J * step = -g
        try:
            step = np.linalg.solve(J, -gval)
        except np.linalg.LinAlgError:
            # fallback: least squares if singular/ill-conditioned
            step, *_ = np.linalg.lstsq(J, -gval, rcond=None)

        # Optional damping (backtracking line search)
        alpha = 1.0
        x_new = x + alpha * step
        if clamp01:
            x_new = np.clip(x_new, 0.0, 1.0)

        if damping:
            # accept if reduces ||g||; else shrink alpha
            for _ in range(20):
                g_new = g(x_new)
                ng_new = float(np.linalg.norm(g_new, ord=2))
                if np.isfinite(ng_new) and ng_new < ng:
                    break
                alpha *= 0.5
                x_new = x + alpha * step
                if clamp01:
                    x_new = np.clip(x_new, 0.0, 1.0)
            else:
                # damping failed to improve
                g_new = g(x_new)
                ng_new = float(np.linalg.norm(g_new, ord=2))
        else:
            g_new = g(x_new)
            ng_new = float(np.linalg.norm(g_new, ord=2))

        x, gval, ng = x_new, g_new, ng_new

        if verbose:
            condJ = np.linalg.cond(J) if np.all(np.isfinite(J)) else np.inf
            # print(
            #     f"[Newton] iter={it} alpha={alpha:.3e} x={x} "
            #     f"||g||={ng:.3e} cond(J)={condJ:.3e} g={gval}"
            # )

    return x, {
        "converged": False,
        "iter": max_iter,
        "residual_norm": ng,
        "residual": gval,
    }