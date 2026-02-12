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
    threshold: int,
    network_type: str,
    k: int,
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
        k_path = f"kOutMin={k}"
    elif network_type == "DirectedCMInPow":
        k_path = f"kInMin={k}"
    else:
        raise ValueError(f"Unknown network type: {network_type}")
    
    output_path = os.path.abspath(
        os.path.join(
            '..',
            f'app/out/fastsar/{network_type}/threshold={threshold}/N={N}/{k_path}'
        )
    )
    if not has_any_file_recursive(output_path):
        output_path = output_path.replace("app/out", "out")

    return output_path


def _load_and_clean_dataframe(file_path: str, batch_index: int, is_final: bool = False) -> Optional[pd.DataFrame]:
    """結果ファイルを読み込み、不完全な反復を削除する。
    
    Args:
        file_path: CSVファイルのパス
        batch_index: バッチインデックス（ログ用）
        is_final: 最終状態データかどうか
        
    Returns:
        クリーンアップされたDataFrame、ファイルが存在しない場合はNone
    """
    if not os.path.exists(file_path):
        return None
    
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"ファイル読み込みエラー (batch={batch_index}): {e}")
        return None
    
    if df.empty:
        return None
    
    # 不完全な反復を検出して削除（時系列データの場合のみ）
    if is_final and 'itr' in df.columns:
        max_itr = int(df['itr'].max())
        for itr in range(max_itr):
            count_itr = (df['itr'] == itr).sum()
            count_next = (df['itr'] == itr + 1).sum()
            
            if count_itr != count_next and count_next > 0:
                # 不完全な反復を削除
                df = df[df['itr'] != itr + 1]
                print(f"削除しました: batch={batch_index}, itr={itr}")
                break
    
    return df


def _load_all_dataframes(output_path: str, batch_size: int, is_final: bool = False) -> Tuple[pd.DataFrame, List[str]]:
    """すべての結果ファイルを読み込み、結合する。
    
    Args:
        output_path: 出力ディレクトリのパス
        batch_size: バッチサイズ
        is_final: 最終状態データかどうか
        
    Returns:
        結合されたDataFrameと読み込まれたファイル名のリスト
    """
    
    df_all = pd.DataFrame()
    existing_files: List[str] = []
    
    for i in range(batch_size):
        index = f"{i:02d}"
        file_name = f'results_{index}.csv'
        file_path = os.path.join(output_path, file_name)
        
        df = _load_and_clean_dataframe(file_path, i, is_final)
        if df is None:
            continue
        
        # バッチインデックスを調整（最初のバッチ以外、時系列データの場合のみ）
        if not df_all.empty and 'itr' in df.columns:
            max_itr = int(df_all['itr'].max())
            df['itr'] = df['itr'] + max_itr + 1
        
        df_all = pd.concat([df_all, df], ignore_index=True)
        existing_files.append(file_name)
    
    print(f"読み込んだファイル: {existing_files}")
    return df_all, existing_files


def _compute_result_final_state(
    df_all: pd.DataFrame,
    dfs_by_itr: Dict[int, pd.DataFrame],
) -> Dict[str, Any]:
    """最終状態データから高次元配列を構築する。
    
    Args:
        df_all: すべてのデータを含むDataFrame
        dfs_by_itr: 反復ごとにグループ化されたDataFrameの辞書
        
    Returns:
        高次元配列とパラメータ値の辞書
    """
    # パラメータ列を自動検出
    param_columns = ['itr', 'rho0', 'lambdaDirected', 'lambdaNondirected', 'mu']
    value_columns = ['A', 'R', 'initialAdoptedTime', 'finalAdoptedTime']
    
    # 存在するパラメータ列のみを使用
    available_param_columns = [col for col in param_columns if col in df_all.columns]
    available_value_columns = [col for col in value_columns if col in df_all.columns]
    
    # 各パラメータのユニークな値を取得
    param_values = {}
    param_to_idx = {}
    for col in available_param_columns:
        if col == 'itr':
            continue  # itrは別途処理
        unique_vals = sorted(df_all[col].unique())
        param_values[col] = np.array(unique_vals)
        param_to_idx[col] = {val: idx for idx, val in enumerate(unique_vals)}
    
    # 配列の次元を決定: (itr, rho0, lambdaDirected, lambdaNondirected, mu, ...)
    dims = [len(dfs_by_itr)]
    for col in available_param_columns:
        if col != 'itr':
            dims.append(len(param_values[col]))
    
    # 高次元配列を初期化
    arrays = {}
    for value_col in available_value_columns:
        arrays[value_col] = np.full(tuple(dims), np.nan)
    
    # データを配列に格納
    for itr_idx, (itr, df) in enumerate(tqdm(dfs_by_itr.items(), desc="配列構築")):
        for _, row in df.iterrows():
            indices = [itr_idx]
            for col in available_param_columns:
                if col != 'itr':
                    val = row[col]
                    indices.append(param_to_idx[col][val])
            
            for value_col in available_value_columns:
                if value_col in row:
                    arrays[value_col][tuple(indices)] = row[value_col]
    
    # delta_timeを計算（存在する場合）
    if 'initialAdoptedTime' in arrays and 'finalAdoptedTime' in arrays:
        arrays['delta_time'] = arrays['finalAdoptedTime'] - arrays['initialAdoptedTime']
    
    result = {
        'arrays': arrays,
        'param_values': param_values,
        'param_columns': available_param_columns,
        'value_columns': available_value_columns,
    }
    
    return result


def _compute_result_time_series(
    df_all: pd.DataFrame,
    dfs_by_itr: Dict[int, pd.DataFrame],
) -> Dict[str, Any]:
    """時系列データから高次元配列を構築する。
    
    Args:
        df_all: すべてのデータを含むDataFrame
        dfs_by_itr: 反復ごとにグループ化されたDataFrameの辞書
        
    Returns:
        高次元配列とパラメータ値の辞書
    """
    # パラメータ列を自動検出
    param_columns = ['itr', 'rho0', 'lambdaDirected', 'lambdaNondirected', 'mu', 'time']
    value_columns = ['A', 'R', 'S']
    
    # 存在するパラメータ列のみを使用
    available_param_columns = [col for col in param_columns if col in df_all.columns]
    available_value_columns = [col for col in value_columns if col in df_all.columns]
    
    # 各パラメータのユニークな値を取得（time以外）
    param_values = {}
    param_to_idx = {}
    for col in available_param_columns:
        if col == 'itr' or col == 'time':
            continue
        unique_vals = sorted(df_all[col].unique())
        param_values[col] = np.array(unique_vals)
        param_to_idx[col] = {val: idx for idx, val in enumerate(unique_vals)}
    
    # timeのユニークな値を取得（全データから）
    if 'time' in available_param_columns:
        time_values = sorted(df_all['time'].unique())
        param_values['time'] = np.array(time_values)
        # timeは浮動小数点なので、許容誤差を考慮してマッピング
        time_to_idx = {}
        for idx, t in enumerate(time_values):
            time_to_idx[t] = idx
        param_to_idx['time'] = time_to_idx
    
    # 配列の次元を決定: (itr, rho0, lambdaDirected, lambdaNondirected, mu, time, ...)
    dims = [len(dfs_by_itr)]
    for col in available_param_columns:
        if col != 'itr':
            dims.append(len(param_values[col]))
    
    # 高次元配列を初期化
    arrays = {}
    for value_col in available_value_columns:
        arrays[value_col] = np.full(tuple(dims), np.nan)
    
    # データを配列に格納
    for itr_idx, (itr, df) in enumerate(tqdm(dfs_by_itr.items(), desc="配列構築")):
        for _, row in df.iterrows():
            indices = [itr_idx]
            for col in available_param_columns:
                if col != 'itr':
                    val = row[col]
                    if col == 'time':
                        # timeの場合は最も近い値を探す（浮動小数点誤差を考慮）
                        if val in param_to_idx[col]:
                            indices.append(param_to_idx[col][val])
                        else:
                            # 最も近い時間インデックスを見つける
                            time_idx = np.argmin(np.abs(param_values[col] - val))
                            indices.append(time_idx)
                    else:
                        indices.append(param_to_idx[col][val])
            
            for value_col in available_value_columns:
                if value_col in row:
                    arrays[value_col][tuple(indices)] = row[value_col]
    
    result = {
        'arrays': arrays,
        'param_values': param_values,
        'param_columns': available_param_columns,
        'value_columns': available_value_columns,
    }
    
    return result


def configure_result(output_path: str, batch_size: int, is_final: bool = False) -> Dict[str, Any]:
    """結果ファイルを読み込み、高次元配列として管理する。
    
    Args:
        output_path: 出力ディレクトリのパス
        batch_size: バッチサイズ
        is_final: 最終状態データかどうか（Falseの場合は時系列データ）
        
    Returns:
        高次元配列とパラメータ値を含む辞書
        - arrays: 各値（A, R, S等）の高次元配列
        - param_values: 各パラメータのユニークな値の配列
        - param_columns: パラメータ列名のリスト
        - value_columns: 値列名のリスト
    """
    df_all, existing_files = _load_all_dataframes(output_path, batch_size, is_final)
    
    if df_all.empty:
        raise ValueError("読み込めるデータファイルがありませんでした")
    
    if 'itr' not in df_all.columns:
        raise ValueError("CSVファイルに'itr'列が存在しません")
    
    dfs_by_itr = {itr: sub_df for itr, sub_df in df_all.groupby('itr')}
    
    if is_final:
        result = _compute_result_final_state(df_all, dfs_by_itr)
    else:
        result = _compute_result_time_series(df_all, dfs_by_itr)
    
    # パラメータ値の情報を出力
    print(f"読み込んだ反復数: {len(dfs_by_itr)}")
    for col, vals in result['param_values'].items():
        print(f"{col}: {len(vals)}個のユニークな値")
    print(f"値列: {result['value_columns']}")
    
    return result

