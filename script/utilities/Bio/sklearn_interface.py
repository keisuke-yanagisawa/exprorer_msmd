import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from sklearn.utils.validation import check_is_fitted
from sklearn.base import TransformerMixin, BaseEstimator
import numpy.typing as npt

"""
BioPythonの関数をsklearnのモデルのように利用する関数/クラス群。

last update: 21 Jun, 2021
Authors: Keisuke Yanagisawa
"""

__all__ = [
    "SuperImposer"
]


class SuperImposer(TransformerMixin, BaseEstimator):
    """
    構造重ね合わせを行うBioPythonのクラスを
    scikit-learnのインターフェースでwrapしたクラス。
    """

    def __init__(self):
        pass

    def _reset(self):
        if hasattr(self, "rot_"):
            del self.rot_
            del self.tran_

    def _superimpose(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> None:
        sup = SVDSuperimposer()
        sup.set(reference_coords, coords)
        sup.run()
        self.rot_, self.tran_ = sup.get_rotran()

    def fit(self, coords: npt.ArrayLike, reference_coords: npt.ArrayLike) -> None:
        """
        与えられた2つの点群をなるべく重ねるような並行・回転移動を算出します。

        与えられた2つの点群はそれぞれ対応関係があることを仮定します。
        すなわち、それぞれの0番目の要素同士がなるべく重なるように、
        1番目の要素同士がなるべく重なるように…と重ね合わせを行います。

        Parameters
        ----------
        coords : list
            重ね合わせのために移動させる点群
        reference_coords : list
            重ね合わせ先の点群

        Returns
        -------
        SuperImposer
            fit済みのオブジェクト
        """
        self._reset()
        self._superimpose(coords, reference_coords)
        return self

    def transform(self, coords: npt.ArrayLike) -> npt.ArrayLike:
        """
        fit()で計算された並進・回転に基づいて
        与えられた点群を移動させます。

        Parameters
        ----------
        coords : list
            移動させる点群
        """
        check_is_fitted(self)
        coords = np.array(coords)
        return np.dot(coords, self.rot_) + self.tran_

    def inverse_transform(self, coords: npt.ArrayLike) -> npt.ArrayLike:
        """
        逆方向の移動を行います。

        Parameters
        ----------
        coords : list
            transform()した後の点群

        Returns
        -------
        np.array
            transform()する前の点群座標
        """
        coords = np.array(coords)
        check_is_fitted(self)
        return np.dot(coords - self.tran_, np.linalg.inv(self.rot_))
