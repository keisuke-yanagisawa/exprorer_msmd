#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Quantitative Inverse MSMDの結果を比較するスクリプト
"""

import argparse
import os
import re


def parse_msmd_result(result_file):
    """
    MSMDシミュレーション結果ファイルを解析する
    
    Parameters
    ----------
    result_file : str
        結果ファイルのパス
    
    Returns
    -------
    dict
        解析結果
    """
    result = {}
    
    with open(result_file, 'r') as f:
        content = f.read()
    
    # 合致度スコアを抽出
    score_match = re.search(r'合致度スコア:\s*([\d\.]+)', content)
    if score_match:
        result['score'] = float(score_match.group(1))
    
    # 結合強度を抽出
    strength_match = re.search(r'結合強度[差]?:\s*([\d\.\-]+)', content)
    if strength_match:
        result['strength'] = float(strength_match.group(1))
    
    # 結論を抽出
    conclusion_match = re.search(r'結論:\s*(.+)', content)
    if conclusion_match:
        result['conclusion'] = conclusion_match.group(1)
    
    return result


def compare_results(result1, result2):
    """
    2つのMSMDシミュレーション結果を比較する
    
    Parameters
    ----------
    result1 : dict
        1つ目の結果
    result2 : dict
        2つ目の結果
    
    Returns
    -------
    dict
        比較結果
    """
    comparison = {}
    
    # スコアの差を計算
    if 'score' in result1 and 'score' in result2:
        comparison['score_diff'] = result2['score'] - result1['score']
    
    # 結合強度の差を計算
    if 'strength' in result1 and 'strength' in result2:
        comparison['strength_diff'] = result2['strength'] - result1['strength']
    
    # 結論を比較
    if 'conclusion' in result1 and 'conclusion' in result2:
        comparison['conclusion_diff'] = (result1['conclusion'] != result2['conclusion'])
        comparison['conclusion1'] = result1['conclusion']
        comparison['conclusion2'] = result2['conclusion']
    
    return comparison


def main():
    """
    メイン関数
    """
    parser = argparse.ArgumentParser(description='Quantitative Inverse MSMDの結果を比較する')
    parser.add_argument('--result1', required=True, help='1つ目の結果ファイルのパス')
    parser.add_argument('--result2', required=True, help='2つ目の結果ファイルのパス')
    parser.add_argument('--output', required=True, help='出力ファイルのパス')
    args = parser.parse_args()
    
    # 結果ファイルを解析
    result1 = parse_msmd_result(args.result1)
    result2 = parse_msmd_result(args.result2)
    
    # 結果を比較
    comparison = compare_results(result1, result2)
    
    # 比較結果を出力
    with open(args.output, 'w') as f:
        f.write("# Quantitative Inverse MSMD結果の比較\n\n")
        
        f.write("## 元の分子\n")
        f.write(f"- 結果ファイル: {args.result1}\n")
        if 'score' in result1:
            f.write(f"- 合致度スコア: {result1['score']}\n")
        if 'strength' in result1:
            f.write(f"- 結合強度: {result1['strength']}\n")
        if 'conclusion' in result1:
            f.write(f"- 結論: {result1['conclusion']}\n")
        f.write("\n")
        
        f.write("## 置換後の分子\n")
        f.write(f"- 結果ファイル: {args.result2}\n")
        if 'score' in result2:
            f.write(f"- 合致度スコア: {result2['score']}\n")
        if 'strength' in result2:
            f.write(f"- 結合強度: {result2['strength']}\n")
        if 'conclusion' in result2:
            f.write(f"- 結論: {result2['conclusion']}\n")
        f.write("\n")
        
        f.write("## 比較結果\n")
        if 'score_diff' in comparison:
            f.write(f"- 合致度スコアの差（置換後 - 元）: {comparison['score_diff']}\n")
            if comparison['score_diff'] > 0:
                f.write("  - 置換後の方が合致度が高い\n")
            elif comparison['score_diff'] < 0:
                f.write("  - 元の方が合致度が高い\n")
            else:
                f.write("  - 合致度に差はない\n")
        
        if 'strength_diff' in comparison:
            f.write(f"- 結合強度の差（置換後 - 元）: {comparison['strength_diff']}\n")
            if comparison['strength_diff'] > 0:
                f.write("  - 置換後の方が結合が強い\n")
            elif comparison['strength_diff'] < 0:
                f.write("  - 元の方が結合が強い\n")
            else:
                f.write("  - 結合強度に差はない\n")
        
        if 'conclusion_diff' in comparison:
            if comparison['conclusion_diff']:
                f.write("- 結論が異なります\n")
                f.write(f"  - 元の結論: {comparison['conclusion1']}\n")
                f.write(f"  - 置換後の結論: {comparison['conclusion2']}\n")
            else:
                f.write("- 結論は同じです\n")
        
        f.write("\n## 総合評価\n")
        if 'strength_diff' in comparison:
            if comparison['strength_diff'] > 0:
                f.write("4番目のCをNに置換することで、結合親和性が向上すると予測されます。\n")
            elif comparison['strength_diff'] < 0:
                f.write("4番目のCをNに置換することで、結合親和性が低下すると予測されます。\n")
            else:
                f.write("4番目のCをNに置換しても、結合親和性に大きな変化はないと予測されます。\n")
    
    print(f"比較結果を {args.output} に保存しました")


if __name__ == "__main__":
    main()