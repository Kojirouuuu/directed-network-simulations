"""
ネットワークタイプに応じたパス構成などの switch ロジックを一元化するユーティリティ。
Java の SwitchUtils.buildNetworkPath と同等のパス形式を提供します。
"""

from pathlib import Path
from typing import Optional


def _require(value: Optional[object], network_type: str, param: str) -> object:
    if value is None:
        raise ValueError(f"{network_type} requires {param}")
    return value


def build_network_path(
    network_type: str,
    N: int,
    *,
    kd_ave: Optional[int] = None,
    ku_ave: Optional[float] = None,
    ku_min: Optional[int] = None,
    ku_max: Optional[int] = None,
    k_in_min: Optional[int] = None,
    k_in_max: Optional[int] = None,
    k_out_min: Optional[int] = None,
    k_out_max: Optional[int] = None,
    kd_min: Optional[int] = None,
    kd_max: Optional[int] = None,
    m0: Optional[int] = None,
    m: Optional[int] = None,
    gamma: Optional[float] = None,
    swap_num: Optional[int] = None,
    gamma_in: Optional[float] = None,
    gamma_out: Optional[float] = None,
    corr_a: Optional[float] = None,
) -> Path:
    """
    ネットワークパス（シミュレーションディレクトリ以下の相対パス）を構築する。
    使用されないパラメータは None でよい。ネットワークタイプに応じて必要なものだけ参照される。

    Args:
        network_type: ネットワークタイプ（DirectedCM, DirectedCMInPow, DirectedCMOutPow,
            PowPow, SameInOut, SchwartzDirectedSF, CM, ER, BA, DirectedBA, ego-Twitter, rev-ego-Twitter, gplus, rev-gplus, higgs-social, rev-higgs-social）
        N: 頂点数
        kd_ave: DirectedCM 用（None 可）
        ku_ave: ER 用（None 可）
        ku_min, ku_max: CM 用（None 可）
        k_in_min, k_in_max: DirectedCMInPow, SchwartzDirectedSF 用（None 可）
        k_out_min, k_out_max: DirectedCMOutPow, SchwartzDirectedSF 用（None 可）
        kd_min, kd_max: PowPow, SameInOut 用（None 可）。CM では ku_min, ku_max として使用
        m0, m: DirectedBA, BA 用（None 可）
        gamma: DirectedCMInPow, DirectedCMOutPow, PowPow, SameInOut 用（None 可）
        swap_num: PowPow 用（None 可）
        gamma_in, gamma_out, corr_a: SchwartzDirectedSF 用（None 可）

    Returns:
        ネットワークパス（例: DirectedCMInPow/N=500000/gamma=2.50/kInMin=5/kInMax=707）

    Raises:
        ValueError: 不明なネットワークタイプまたは必要なパラメータが None の場合
    """
    if network_type == "DirectedCM":
        _require(kd_ave, "DirectedCM", "kd_ave")
        network_specific = f"kdAve={kd_ave}"
    elif network_type == "DirectedCMInPow":
        _require(k_in_min, "DirectedCMInPow", "k_in_min")
        _require(k_in_max, "DirectedCMInPow", "k_in_max")
        _require(gamma, "DirectedCMInPow", "gamma")
        network_specific = f"gamma={gamma:.2f}/kInMin={k_in_min}/kInMax={k_in_max}"
    elif network_type == "DirectedCMOutPow":
        _require(k_out_min, "DirectedCMOutPow", "k_out_min")
        _require(k_out_max, "DirectedCMOutPow", "k_out_max")
        _require(gamma, "DirectedCMOutPow", "gamma")
        network_specific = f"gamma={gamma:.2f}/kOutMin={k_out_min}/kOutMax={k_out_max}"
    elif network_type == "PowPow":
        _require(kd_min, "PowPow", "kd_min")
        _require(kd_max, "PowPow", "kd_max")
        _require(gamma, "PowPow", "gamma")
        _require(swap_num, "PowPow", "swap_num")
        network_specific = f"gamma={gamma:.2f}/kdMin={kd_min}/kdMax={kd_max}/swapNum={swap_num}"
    elif network_type == "SameInOut":
        _require(kd_min, "SameInOut", "kd_min")
        _require(kd_max, "SameInOut", "kd_max")
        _require(gamma, "SameInOut", "gamma")
        network_specific = f"gamma={gamma:.2f}/kdMin={kd_min}/kdMax={kd_max}"
    elif network_type == "SchwartzDirectedSF":
        _require(gamma_in, "SchwartzDirectedSF", "gamma_in")
        _require(gamma_out, "SchwartzDirectedSF", "gamma_out")
        _require(k_in_min, "SchwartzDirectedSF", "k_in_min")
        _require(k_in_max, "SchwartzDirectedSF", "k_in_max")
        _require(k_out_min, "SchwartzDirectedSF", "k_out_min")
        _require(k_out_max, "SchwartzDirectedSF", "k_out_max")
        _require(corr_a, "SchwartzDirectedSF", "corr_a")
        network_specific = (
            f"gammaIn={gamma_in:.2f}/gammaOut={gamma_out:.2f}"
            f"/kInMin={k_in_min}/kInMax={k_in_max}"
            f"/kOutMin={k_out_min}/kOutMax={k_out_max}"
            f"/corrA={corr_a:.2f}"
        )
    elif network_type == "CM":
        _require(ku_min, "CM", "ku_min")
        _require(ku_max, "CM", "ku_max")
        _require(gamma, "CM", "gamma")
        network_specific = f"gamma={gamma:.2f}/kuMin={ku_min}/kuMax={ku_max}"
    elif network_type == "ER":
        _require(ku_ave, "ER", "ku_ave")
        network_specific = f"kuAve={ku_ave:.2f}"
    elif network_type in ("DirectedBA", "BA"):
        _require(m0, network_type, "m0")
        _require(m, network_type, "m")
        network_specific = f"m0={m0}/m={m}"
    elif network_type == "ego-Twitter":
        network_specific = ""
    elif network_type == "rev-ego-Twitter":
        network_specific = ""
    elif network_type == "gplus":
        network_specific = ""
    elif network_type == "rev-gplus":
        network_specific = ""
    elif network_type == "higgs-social":
        network_specific = ""
    elif network_type == "rev-higgs-social":
        network_specific = ""
    else:
        raise ValueError(f"Unknown network type: {network_type}")

    return Path(f"{network_type}/N={N}/{network_specific}")
