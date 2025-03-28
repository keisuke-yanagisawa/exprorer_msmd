#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DXファイルの内容を2次元スライスとして可視化するスクリプト
"""

import sys

import gridData
import matplotlib.pyplot as plt
import numpy as np


def visualize_dx_slices(dx_file, output_prefix):
    """
    DXファイルの内容を3つの主要平面（xy, xz, yz）でスライスして可視化する
    
    Parameters:
    -----------
    dx_file : str
        DXファイルのパス
    output_prefix : str
        出力画像ファイルの接頭辞
    """
    # DXファイルを読み込む
    grid = gridData.Grid(dx_file)
    
    # グリッドの形状とサイズを取得
    grid_shape = grid.grid.shape
    
    # 中心のインデックスを計算
    center_x = grid_shape[0] // 2
    center_y = grid_shape[1] // 2
    center_z = grid_shape[2] // 2
    
    # 3つの主要平面でスライス
    xy_slice = grid.grid[:, :, center_z]
    xz_slice = grid.grid[:, center_y, :]
    yz_slice = grid.grid[center_x, :, :]
    
    # 可視化
    plt.figure(figsize=(18, 6))
    
    # xy平面
    plt.subplot(1, 3, 1)
    plt.imshow(xy_slice.T, origin='lower', cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title(f'XY Plane (Z={center_z})')
    plt.xlabel('X')
    plt.ylabel('Y')
    
    # xz平面
    plt.subplot(1, 3, 2)
    plt.imshow(xz_slice.T, origin='lower', cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title(f'XZ Plane (Y={center_y})')
    plt.xlabel('X')
    plt.ylabel('Z')
    
    # yz平面
    plt.subplot(1, 3, 3)
    plt.imshow(yz_slice.T, origin='lower', cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Value')
    plt.title(f'YZ Plane (X={center_x})')
    plt.xlabel('Y')
    plt.ylabel('Z')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_slices.png')
    plt.close()
    
    print(f"スライス画像を保存しました: {output_prefix}_slices.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("使用方法: python visualize_dx.py <dx_file> [output_prefix]")
        sys.exit(1)
    
    dx_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else dx_file.split('.')[0]
    
    visualize_dx_slices(dx_file, output_prefix)