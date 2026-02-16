import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd
from tqdm import tqdm

def has_any_file_recursive(dir_path: str) -> bool:
    return any(f.is_file() for f in Path(dir_path).rglob("*"))

def configure_output_path(
    N: int,
    abst: str,
    threshold: int,
    network_type: str,
    k: float,
) -> str:
    """出力パスを構築する。
    
    Args:
        N: ネットワークサイズ
        threshold: 閾値
        network_type: ネットワークタイプ
        
    Returns:
        出力ディレクトリの絶対パス
    """

    if network_type == "DirectedCMOutPow":
        k = int(k)
        k_path = f"kOutMin={k}"
    elif network_type == "DirectedCMInPow":
        k = int(k)
        k_path = f"kInMin={k}"
    elif network_type == "ER":
        k_path = f"z={k:.2f}"
    else:
        raise ValueError(f"Unknown network type: {network_type}")
    
    output_path = os.path.abspath(
        os.path.join(
            '..',
            f'app/out/fastsar/{abst}/{network_type}/threshold={threshold}/N={N}/{k_path}'
        )
    )
    if not has_any_file_recursive(output_path):
        output_path = output_path.replace("app/out", "out")

    return output_path

