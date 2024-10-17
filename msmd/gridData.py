from typing import Optional

from gridData import Grid


def generate_grid(grid, PITCH: int, origin: Optional[list[int]] = None, center: Optional[list[int]] = None) -> Grid:
    """与えられたグリッド情報を基に、gridData.Gridを作成する。

    Parameters
    -----------
    grid: np.ndarray
        3次元のnumpy.array
    PITCH: int
        gridの幅。x, y, zが等しい、立方体構造であることを仮定している
    origin: Optional[list[int]]
        gridの開始点。中心点ではなく、グリッドにおける(0, 0, 0)、すなわち一番左上の座標を入力する。
        originとcenterはいずれか一方のみを指定する。
    center: Optional[list[int]]
        gridの中心点。
        originとcenterはいずれか一方のみを指定する。

    Returns
    --------
    g: gridData.Grid
        生成されたGridオブジェクト

    Raises
    -------
    ValueError
        originとcenterのどちらも指定されていない場合、およびoriginとcenterの両方が指定されていた場合

        入力されたgridが立方体でない場合
    """

    if origin is None and center is None:
        raise ValueError("Either origin or center must be specified.")
    if origin is not None and center is not None:
        raise ValueError("Both origin and center cannot be specified at the same time.")

    if grid.shape[0] != grid.shape[1] or grid.shape[0] != grid.shape[2]:
        raise ValueError("The grid must be a cube.")

    SIZE = grid.shape[0]
    g = Grid()
    g.delta = [PITCH, PITCH, PITCH]
    g.origin = (
        origin
        if origin is not None
        else [center[0] - SIZE * PITCH / 2, center[1] - SIZE * PITCH / 2, center[2] - SIZE * PITCH / 2]
    )
    g.n = [SIZE, SIZE, SIZE]
    g.grid = grid
    return g
