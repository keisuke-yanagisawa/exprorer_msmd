#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
画像ファイルを表示し、分子の向きを確認するスクリプト
"""

import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont


def analyze_image(image_path):
    """
    画像ファイルを解析し、分子の向きを確認する

    Parameters
    ----------
    image_path : str
        画像ファイルのパス

    Returns
    -------
    dict
        画像の解析結果
    """
    # 画像を読み込む
    img = Image.open(image_path)
    
    # 画像の情報を取得
    width, height = img.size
    format = img.format
    mode = img.mode
    
    # 画像をNumPy配列に変換
    img_array = np.array(img)
    
    # 赤色と青色のピクセルを検出（分子のハイライト部分）
    red_pixels = np.sum((img_array[:, :, 0] > 200) & (img_array[:, :, 1] < 100) & (img_array[:, :, 2] < 100))
    blue_pixels = np.sum((img_array[:, :, 0] < 100) & (img_array[:, :, 1] < 100) & (img_array[:, :, 2] > 200))
    
    # 結果を返す
    return {
        "path": image_path,
        "size": (width, height),
        "format": format,
        "mode": mode,
        "red_pixels": red_pixels,
        "blue_pixels": blue_pixels,
        "total_pixels": width * height
    }

def create_comparison_image(image_paths, output_path):
    """
    複数の画像を並べて比較画像を作成する

    Parameters
    ----------
    image_paths : list of str
        画像ファイルのパスのリスト
    output_path : str
        出力画像ファイルのパス
    """
    # 画像を読み込む
    images = [Image.open(path) for path in image_paths]
    
    # 画像のサイズを取得
    widths, heights = zip(*(img.size for img in images))
    
    # 比較画像のサイズを計算
    total_width = sum(widths)
    max_height = max(heights)
    
    # 新しい画像を作成（白色は(255, 255, 255)）
    comparison_img = Image.new('RGB', (total_width, max_height + 30), (255, 255, 255))
    
    # 画像を並べて配置
    x_offset = 0
    for i, img in enumerate(images):
        comparison_img.paste(img, (x_offset, 0))
        
        # 画像の下にファイル名を表示
        draw = ImageDraw.Draw(comparison_img)
        filename = os.path.basename(image_paths[i])
        draw.text((x_offset + 10, max_height + 10), filename, fill=(0, 0, 0))
        
        x_offset += img.width
    
    # 比較画像を保存
    comparison_img.save(output_path)
    
    print(f"比較画像を {output_path} に保存しました")

def main():
    """
    メイン関数
    """
    # 画像ファイルのパス
    image_paths = [
        "./output/superimposed_trimethylbenzene.png",
        "./output/superimposed_dimethylbenzene.png",
        "./output/superimposed_methylbenzene.png"
    ]
    
    # 各画像を解析
    for path in image_paths:
        if os.path.exists(path):
            result = analyze_image(path)
            print(f"\n画像ファイル: {result['path']}")
            print(f"  サイズ: {result['size'][0]} x {result['size'][1]} ピクセル")
            print(f"  フォーマット: {result['format']}")
            print(f"  モード: {result['mode']}")
            print(f"  赤色ピクセル数: {result['red_pixels']} ({result['red_pixels'] / result['total_pixels'] * 100:.2f}%)")
            print(f"  青色ピクセル数: {result['blue_pixels']} ({result['blue_pixels'] / result['total_pixels'] * 100:.2f}%)")
        else:
            print(f"\n画像ファイル {path} が見つかりません")
    
    # 比較画像を作成
    output_path = "./output/comparison.png"
    create_comparison_image(image_paths, output_path)

if __name__ == "__main__":
    main()