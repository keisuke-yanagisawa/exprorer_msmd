# coding: utf-8

"""
回転を処理するscipy.spatial.transform.Rotation モジュールに関係する関数群。

version: 1.1.1
last update: 17 Feb, 2021
Authors: Keisuke Yanagisawa
"""

import warnings

import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R


def __get_a_nonzero_vector(coords, eps=1e-5, priority=None):
    """
    return an coord index which is not at the origin.
    if all coordinates are at the origin, return None.
    """
    if priority is None:
        priority = list(range(coords.shape[0]))

    for idx in priority:
        coord = coords[idx]
        if norm(coord) > eps:
            return idx
    return None


def __get_a_li_vector(coords, ref, eps=1e-5, priority=None):
    """
    return an coord index which is linearly independent with provided coord (ref).
    if there is no linearly independent coord, it returns None.
    """
    if priority is None:
        priority = list(range(coords.shape[0]))

    for idx in priority:
        coord = coords[idx]
        if np.abs(np.abs(coord.dot(ref)) - norm(coord) * norm(ref)) > eps:
            return idx
    return None


def __standardize_direction(coords, eps=1e-5, priority=None):
    """
    return coordinates which are rotated to specific direction.
    """
    if priority is None:
        priority = list(range(coords.shape[0]))

    src = np.zeros((2, 3))
    dst = np.zeros((2, 3))

    # determine a coord which will be put on x-axis
    first = __get_a_nonzero_vector(coords, eps, priority=priority)
    if first is None:
        return coords  # do nothing
    dst[0] = [norm(coords[first]), 0, 0]
    src[0] = coords[first]

    # determine a coord which will be put on x-y plane
    second = __get_a_li_vector(coords, src[0], eps, priority=priority)
    if second is not None:
        cos_theta = coords[first].dot(coords[second]) / norm(coords[first]) / norm(coords[second])
        sin_theta = np.sqrt(1 - cos_theta**2)
        dst[1] = [norm(coords[second]) * cos_theta, norm(coords[second]) * sin_theta, 0]
        src[1] = coords[second]

    # calculate rotation matrix
    rot, rmsd = R.align_vectors(dst, src)
    if rmsd > eps:
        warnings.warn("R.align_vectors(dst, src) outputs non-zero rmsd")
        warnings.warn(f"rmsd is: {rmsd}")

    # apply rotation matrix & return rotated coords
    return rot.apply(coords)


# output: transformed coordinates, displacement, rotation matrix
def _standardize_coords(coords, priority=None, origin=None):
    """
    複数の3次元座標（点群）を入力として受け取り、
    その点群の座標に応じて向きを1つに定める。
    点群の順番が異なっている場合は異なる結果となるため注意が必要。
    """
    if priority is None:
        priority = list(range(len(coords)))

    tmp = coords - coords[priority[0]]  # put `priority[0]` to origin point temporaly
    tmp = __standardize_direction(tmp, priority=priority)

    if origin is None:
        origin_coord = tmp.mean(axis=0)
    else:
        origin_coord = tmp[origin]
    tmp = tmp - origin_coord  # re-put `origin_coord` to origin point

    rot, rmsd = R.align_vectors(tmp + origin_coord, coords - coords[priority[0]])
    if rmsd > 1e-2:
        warnings.warn(f"rmsd = {rmsd:.5f}")

    return tmp, -coords[priority[0]], rot, -origin_coord


class CoordinateStandardizer:
    """
    sklearn風のAPIを備えた点群座標変換クラス。
    fitで点群変換を決定し、transformで点群を別の座標軸に置き換える。
    """

    def __init__(self):
        self.is_fitted = False
        self.d1 = None
        self.rot = None
        self.d2 = None

    def fit(self, coords, priority=None, origin=None):
        self.is_fitted = True
        _, self.d1, self.rot, self.d2 = _standardize_coords(coords, priority, origin)

    def transform(self, coords):
        if not self.is_fitted:
            raise RuntimeError(
                """This CoordinateStandardizer instance is not fitted yet.
            Call 'fit' with appropriate arguments before using this standardizer."""
            )

        tmp = coords + self.d1
        tmp = self.rot.apply(tmp)
        tmp += self.d2
        return tmp

    def fit_transform(self, coords, priority=None, origin=None):
        self.fit(coords, priority, origin)
        return self.transform(coords)
