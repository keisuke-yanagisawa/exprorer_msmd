#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
画像ファイルの内容を確認するスクリプト
"""

import sys

import numpy as np
from PIL import Image


def check_image(image_path):
    """
    画像ファイルの内容を確認する
    
    Parameters
    ----------
    image_path : str
        画像ファイルのパス
    """
    # 画像を読み込む
    img = Image.open(image_path)
    
    # 画像の情報を表示
    print(f"画像ファイル: {image_path}")
    print(f"サイズ: {img.size[0]} x {img.size[1]} ピクセル")
    print(f"フォーマット: {img.format}")
    print(f"モード: {img.mode}")
    
    # 画像の一部を表示（ピクセル値）
    img_array = np.array(img)
    print("\n画像の一部（左上の10x10ピクセル）:")
    print(img_array[0:10, 0:10])
    
    # 画像の色分布を表示
    if img.mode == "RGB":
        r, g, b = img.split()
        print("\n色分布:")
        print(f"赤: 平均={np.mean(r):.2f}, 最小={np.min(r)}, 最大={np.max(r)}")
        print(f"緑: 平均={np.mean(g):.2f}, 最小={np.min(g)}, 最大={np.max(g)}")
        print(f"青: 平均={np.mean(b):.2f}, 最小={np.min(b)}, 最大={np.max(b)}")

def main():
    """
    メイン関数
    """
    if len(sys.argv) < 2:
        print("使用方法: python check_image.py <画像ファイルパス>")
        return
    
    image_path = sys.argv[1]
    check_image(image_path)

if __name__ == "__main__":
    main()