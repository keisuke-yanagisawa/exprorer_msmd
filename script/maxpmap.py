#!/usr/bin/python3

from typing import List

import numpy as np
from gridData import Grid

from script.utilities.logger import logger

VERSION = "1.0.0"


def _check(gs: List[Grid]) -> None:
    """グリッドのサイズとポジションの互換性を確認する内部関数
    
    Raises:
        ValueError: グリッドのサイズやポジションが異なる場合
    """
    if len(gs) == 1:
        return
        
    reference = gs[0]
    ref_shape = reference.grid.shape
    ref_origin = reference.origin
    ref_delta = reference.delta
    
    for g in gs[1:]:
        if g.grid.shape != ref_shape:
            raise ValueError("Grids have different sizes")
        if not np.allclose(g.origin, ref_origin):
            raise ValueError("Grids have different origins")
        if not np.allclose(g.delta, ref_delta):
            raise ValueError("Grids have different deltas")


def grid_max(gs: List[Grid]) -> Grid:
    """複数のグリッドデータから各点の最大値を計算する
    
    Args:
        gs: Gridオブジェクトのリスト
        
    Returns:
        Grid: 最大値を持つ新しいGridオブジェクト
        
    Raises:
        ValueError: グリッドリストが空の場合、またはグリッドのサイズやポジションが異なる場合
    """
    if not gs:
        raise ValueError("Empty grid list")
        
    _check(gs)
    ret = gs[0]
    ret.grid = np.max([g.grid for g in gs], axis=0)
    return ret


def gen_max_pmap(inpaths: List[str], outpath: str) -> str:
    """複数のpmapファイルから最大値のpmapファイルを生成する
    
    Args:
        inpaths: 入力pmapファイル（dx形式）のパスのリスト
        outpath: 出力pmapファイル（dx形式）のパス
        
    Returns:
        str: 出力ファイルのパス
        
    Raises:
        ValueError: 入力ファイルリストが空の場合、またはグリッドのサイズやポジションが異なる場合
    """
    if not inpaths:
        raise ValueError("No input files provided")
        
    gs = [Grid(n) for n in inpaths]
    max_pmap = grid_max(gs)
    max_pmap.export(outpath, type="double")
    return outpath
