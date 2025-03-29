#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
プロファイル関連のクラスと関数を提供するモジュール
"""

import glob
import os
from dataclasses import dataclass
from typing import Final, List, Optional

import numpy as np
import numpy.typing as npt
from gridData import Grid
from tqdm import tqdm


@dataclass
class PreProfile:
    """
    各残基に対するbulk_probasumとprobasum_gridを持つ

    bulk_probasumはその残基が均一に存在する場合の単位体積あたりの存在確率を、フレーム全体で合計したもの
    probasum_gridはMSMDシミュレーションで得られた、その残基が存在する確率のグリッドデータ
    """
    ala_bulk_probasum: float
    ala_probasum_grid: Grid
    cys_bulk_probasum: float
    cys_probasum_grid: Grid
    asp_bulk_probasum: float
    asp_probasum_grid: Grid
    glu_bulk_probasum: float
    glu_probasum_grid: Grid
    phe_bulk_probasum: float
    phe_probasum_grid: Grid
    gly_bulk_probasum: float
    gly_probasum_grid: Grid
    his_bulk_probasum: float
    his_probasum_grid: Grid
    ile_bulk_probasum: float
    ile_probasum_grid: Grid
    lys_bulk_probasum: float
    lys_probasum_grid: Grid
    leu_bulk_probasum: float
    leu_probasum_grid: Grid
    met_bulk_probasum: float
    met_probasum_grid: Grid
    asn_bulk_probasum: float
    asn_probasum_grid: Grid
    pro_bulk_probasum: float
    pro_probasum_grid: Grid
    gln_bulk_probasum: float
    gln_probasum_grid: Grid
    arg_bulk_probasum: float
    arg_probasum_grid: Grid
    ser_bulk_probasum: float
    ser_probasum_grid: Grid
    thr_bulk_probasum: float
    thr_probasum_grid: Grid
    val_bulk_probasum: float
    val_probasum_grid: Grid
    trp_bulk_probasum: float
    trp_probasum_grid: Grid
    tyr_bulk_probasum: float
    tyr_probasum_grid: Grid

    def __add__(self, other: "PreProfile") -> "PreProfile":
        """
        2つのPreProfileを加算する
        """
        return PreProfile(
            ala_bulk_probasum=self.ala_bulk_probasum + other.ala_bulk_probasum,
            ala_probasum_grid=self.ala_probasum_grid + other.ala_probasum_grid,
            cys_bulk_probasum=self.cys_bulk_probasum + other.cys_bulk_probasum,
            cys_probasum_grid=self.cys_probasum_grid + other.cys_probasum_grid,
            asp_bulk_probasum=self.asp_bulk_probasum + other.asp_bulk_probasum,
            asp_probasum_grid=self.asp_probasum_grid + other.asp_probasum_grid,
            glu_bulk_probasum=self.glu_bulk_probasum + other.glu_bulk_probasum,
            glu_probasum_grid=self.glu_probasum_grid + other.glu_probasum_grid,
            phe_bulk_probasum=self.phe_bulk_probasum + other.phe_bulk_probasum,
            phe_probasum_grid=self.phe_probasum_grid + other.phe_probasum_grid,
            gly_bulk_probasum=self.gly_bulk_probasum + other.gly_bulk_probasum,
            gly_probasum_grid=self.gly_probasum_grid + other.gly_probasum_grid,
            his_bulk_probasum=self.his_bulk_probasum + other.his_bulk_probasum,
            his_probasum_grid=self.his_probasum_grid + other.his_probasum_grid,
            ile_bulk_probasum=self.ile_bulk_probasum + other.ile_bulk_probasum,
            ile_probasum_grid=self.ile_probasum_grid + other.ile_probasum_grid,
            lys_bulk_probasum=self.lys_bulk_probasum + other.lys_bulk_probasum,
            lys_probasum_grid=self.lys_probasum_grid + other.lys_probasum_grid,
            leu_bulk_probasum=self.leu_bulk_probasum + other.leu_bulk_probasum,
            leu_probasum_grid=self.leu_probasum_grid + other.leu_probasum_grid,
            met_bulk_probasum=self.met_bulk_probasum + other.met_bulk_probasum,
            met_probasum_grid=self.met_probasum_grid + other.met_probasum_grid,
            asn_bulk_probasum=self.asn_bulk_probasum + other.asn_bulk_probasum,
            asn_probasum_grid=self.asn_probasum_grid + other.asn_probasum_grid,
            pro_bulk_probasum=self.pro_bulk_probasum + other.pro_bulk_probasum,
            pro_probasum_grid=self.pro_probasum_grid + other.pro_probasum_grid,
            gln_bulk_probasum=self.gln_bulk_probasum + other.gln_bulk_probasum,
            gln_probasum_grid=self.gln_probasum_grid + other.gln_probasum_grid,
            arg_bulk_probasum=self.arg_bulk_probasum + other.arg_bulk_probasum,
            arg_probasum_grid=self.arg_probasum_grid + other.arg_probasum_grid,
            ser_bulk_probasum=self.ser_bulk_probasum + other.ser_bulk_probasum,
            ser_probasum_grid=self.ser_probasum_grid + other.ser_probasum_grid,
            thr_bulk_probasum=self.thr_bulk_probasum + other.thr_bulk_probasum,
            thr_probasum_grid=self.thr_probasum_grid + other.thr_probasum_grid,
            val_bulk_probasum=self.val_bulk_probasum + other.val_bulk_probasum,
            val_probasum_grid=self.val_probasum_grid + other.val_probasum_grid,
            trp_bulk_probasum=self.trp_bulk_probasum + other.trp_bulk_probasum,
            trp_probasum_grid=self.trp_probasum_grid + other.trp_probasum_grid,
            tyr_bulk_probasum=self.tyr_bulk_probasum + other.tyr_bulk_probasum,
            tyr_probasum_grid=self.tyr_probasum_grid + other.tyr_probasum_grid,
        )


@dataclass
class Profile:
    """
    各残基に対する、bulkとの存在確率比を保持する
    """
    A: Grid
    C: Grid
    D: Grid
    E: Grid
    F: Grid
    G: Grid
    H: Grid
    I: Grid
    K: Grid
    L: Grid
    M: Grid
    N: Grid
    P: Grid
    Q: Grid
    R: Grid
    S: Grid
    T: Grid
    V: Grid
    W: Grid
    Y: Grid

    def __init__(self, preprofile: PreProfile):
        """
        PreProfileからProfileを作成する
        """
        self.A = preprofile.ala_probasum_grid / preprofile.ala_bulk_probasum
        self.C = preprofile.cys_probasum_grid / preprofile.cys_bulk_probasum
        self.D = preprofile.asp_probasum_grid / preprofile.asp_bulk_probasum
        self.E = preprofile.glu_probasum_grid / preprofile.glu_bulk_probasum
        self.F = preprofile.phe_probasum_grid / preprofile.phe_bulk_probasum
        self.G = preprofile.gly_probasum_grid / preprofile.gly_bulk_probasum
        self.H = preprofile.his_probasum_grid / preprofile.his_bulk_probasum
        self.I = preprofile.ile_probasum_grid / preprofile.ile_bulk_probasum
        self.K = preprofile.lys_probasum_grid / preprofile.lys_bulk_probasum
        self.L = preprofile.leu_probasum_grid / preprofile.leu_bulk_probasum
        self.M = preprofile.met_probasum_grid / preprofile.met_bulk_probasum
        self.N = preprofile.asn_probasum_grid / preprofile.asn_bulk_probasum
        self.P = preprofile.pro_probasum_grid / preprofile.pro_bulk_probasum
        self.Q = preprofile.gln_probasum_grid / preprofile.gln_bulk_probasum
        self.R = preprofile.arg_probasum_grid / preprofile.arg_bulk_probasum
        self.S = preprofile.ser_probasum_grid / preprofile.ser_bulk_probasum
        self.T = preprofile.thr_probasum_grid / preprofile.thr_bulk_probasum
        self.V = preprofile.val_probasum_grid / preprofile.val_bulk_probasum
        self.W = preprofile.trp_probasum_grid / preprofile.trp_bulk_probasum
        self.Y = preprofile.tyr_probasum_grid / preprofile.tyr_bulk_probasum


def generate_grid(grid, PITCH: float, origin: npt.ArrayLike) -> Grid:
    """
    与えられたグリッド情報を基に、gridData.Gridを作成する

    Parameters
    ----------
    grid : numpy.ndarray
        3次元のnumpy.array
    PITCH : float
        gridの幅。x, y, zが等しいことを仮定している
    origin : numpy.ndarray
        gridの開始点。中心点ではなく、グリッドにおける(0, 0, 0)、すなわち一番左上の座標を入力する。

    Returns
    -------
    Grid
        gridData.Grid
    """
    return Grid(
        grid=grid,
        origin=origin,
        delta=[PITCH, PITCH, PITCH],
    )


def profile_cnt_sum(base_dir: str, probe: str = "A00", aa: str = "GLY", EPS: float = 0.1) -> Grid:
    """
    多数のsingle count profileの総和を取ったグリッドオブジェクトを作成する。

    ここで作成されたprofileは、後でlogを取ったりするので、値が完全に0になるのは望ましくない。
    そこで、全ての座標について、EPSだけ値を加算しておく。

    Parameters
    ----------
    base_dir : str
        ベースディレクトリ
    probe : str, optional
        プローブID, by default "A00"
    aa : str, optional
        アミノ酸名, by default "GLY"
    EPS : float, optional
        0を避けるための微小値, by default 0.1

    Returns
    -------
    Grid
        グリッドオブジェクト
    """
    # single profileのパスを全て取得
    grid_paths = [path for path in glob.glob(f"{base_dir}/*/{probe}/system*/simulation/*_{probe}_{aa}_environment.dx")]
    
    if len(grid_paths) == 0:
        raise ValueError(f"there is no grids: {probe}, {aa}")

    # 全てのグリッドの値を足し合わせる
    grid = Grid(grid_paths[0])
    
    for path in tqdm(grid_paths[1:]):
        try:
            tmp = Grid(path)
        except AttributeError as e:
            print("failed to open grid file", path)
            raise e

        grid.grid += tmp.grid
    
    grid.grid += EPS  # log(0) にならないように微小の値を足す
    
    return grid


def bulk_sum(base_dir: str, probe: str = "A00", aa: str = "GLY") -> float:
    """
    多数のbulk countの総和を取る。

    全てのbulk countはtxtファイルに保管されているので、そのファイルを全て開いて、合計値を出す。

    Parameters
    ----------
    base_dir : str
        ベースディレクトリ
    probe : str, optional
        プローブID, by default "A00"
    aa : str, optional
        アミノ酸名, by default "GLY"

    Returns
    -------
    float
        バルクカウントの総和
    """
    # single profileのパスを全て取得
    bulk_paths = [path for path in glob.glob(f"{base_dir}/*/{probe}/system*/simulation/*_{probe}_{aa}_bulk_cnt.txt")]
    
    if len(bulk_paths) == 0:
        raise ValueError(f"there is no bulk cnt txt file: {probe}, {aa}")

    # ファイルの中身をfloatとして読み込んで、総和を取る
    bulk_val = 0
    
    for path in bulk_paths:
        val = float(open(path).read().strip())
        bulk_val += val

    return bulk_val


def is_same_file_cnt(base_dir: str, probe: str = "A00", aa: str = "GLY") -> bool:
    """
    profileファイルとbulk countファイルの数が一致しているかどうかを確認する。

    Parameters
    ----------
    base_dir : str
        ベースディレクトリ
    probe : str, optional
        プローブID, by default "A00"
    aa : str, optional
        アミノ酸名, by default "GLY"

    Returns
    -------
    bool
        一致していればTrue、そうでなければFalse
    """
    grid_paths = [path for path in glob.glob(f"{base_dir}/*/{probe}/system*/simulation/*_{probe}_{aa}_environment.dx")]
    bulk_paths = [path for path in glob.glob(f"{base_dir}/*/{probe}/system*/simulation/*_{probe}_{aa}_bulk_cnt.txt")]

    return len(grid_paths) == len(bulk_paths)


def create_single_profile(traj_path: str, top_path: str, ref_probe_path: str, probe: str, aa: str, output_dir: str = "."):
    """
    単一のプロファイルを作成する

    Parameters
    ----------
    traj_path : str
        トラジェクトリファイルのパス
    top_path : str
        トポロジーファイルのパス
    ref_probe_path : str
        参照プローブファイルのパス
    probe : str
        プローブID
    aa : str
        アミノ酸名
    output_dir : str, optional
        出力ディレクトリ, by default "."
    """
    import pytraj as pt

    # トラジェクトリの読み込み
    traj = pt.iterload(
        filename=traj_path,
        top=top_path,
    )
    
    # 最後の1000フレームだけ使用
    traj = traj[-1000:]
    
    # 参照プローブの読み込み
    ref_probe = pt.load(ref_probe_path)
    
    # プローブの残基IDリストを取得
    probe_resis = set()
    for residue in traj[f":{probe}"].topology.residues:
        probe_resis.add(residue.index)
    
    # ボックスの体積を計算
    box_volume = np.prod(traj.unitcells[0][:3])
    
    # グリッドデータの初期化
    grid_data = None
    bulk_cnt = 0  # bulk空間に存在する残基の数
    
    # 各プローブ残基について処理
    from script.utilities.analysis.matching_score import preprocessing
    
    for probe_resi in tqdm(probe_resis):
        # トラジェクトリの前処理
        ret = preprocessing(traj, ref_probe, probe_resi)
        
        # グリッドの計算
        data = pt.grid(
            ret,
            command=f"80 10 80 10 80 10 ((@/C&@CB)&(!:{probe})&(:{aa}))",
            dtype="dataset",
        ).to_ndarray()
        
        if grid_data is None:
            grid_data = data
        else:
            grid_data += data / box_volume
        
        bulk_cnt += data.sum()
    
    # グリッドの作成と保存
    SIZE, PITCH = 80, 1
    grid = generate_grid(grid_data, PITCH, [-SIZE * PITCH / 2, -SIZE * PITCH / 2, -SIZE * PITCH / 2])
    
    os.makedirs(output_dir, exist_ok=True)
    grid.export(f"{output_dir}/{probe}_{aa}_environment.dx")
    
    # バルクカウントの保存
    with open(f"{output_dir}/{probe}_{aa}_bulk_cnt.txt", "w") as fout:
        fout.write(f"{bulk_cnt}")