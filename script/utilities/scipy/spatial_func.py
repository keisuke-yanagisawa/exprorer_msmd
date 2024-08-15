import numpy as np
import numpy.typing as npt
from scipy.spatial import KDTree
import warnings


def estimate_volume(points: npt.NDArray[np.float_], radii: npt.NDArray[np.float_], granularity=10):
    """
    Calculate the volume of a set of spheres.
    Estimate by the number of occupied voxels.

    Possible error rate of estimation is ~10% without raising warnings.

    input:
      points: list of points
      radii: radius of each sphere
    output:
      volume: volume of the set of spheres
    """

    # generate a grid
    x_max, y_max, z_max = np.max(points, axis=0)
    x_min, y_min, z_min = np.min(points, axis=0)
    x_pitch = np.min([(x_max - x_min) / granularity, 1])
    y_pitch = np.min([(y_max - y_min) / granularity, 1])
    z_pitch = np.min([(z_max - z_min) / granularity, 1])

    if np.min(radii) < np.max([x_pitch, y_pitch, z_pitch]) * 2:  # 2 is a magic number
        warnings.warn(
            "The granularity is too small compared to the radii of spheres. "
            "The volumes of small spheres are underestimated. "
            "Consider increasing the granularity. ", 
            RuntimeWarning)

    # enumerate possible points to be occupied
    xs = np.arange(x_min - np.max(radii), x_max + np.max(radii) + 1, x_pitch)
    ys = np.arange(y_min - np.max(radii), y_max + np.max(radii) + 1, y_pitch)
    zs = np.arange(z_min - np.max(radii), z_max + np.max(radii) + 1, z_pitch)
    if len(xs) * len(ys) * len(zs) > 1e7:
        warnings.warn(
            "The number of occupied voxels is too large. "
            "It will spend a long time to enumerate. ",
            RuntimeWarning)
    xx, yy, zz = np.meshgrid(xs, ys, zs)
    grid_points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T
    # create a kdtree for fast nearest neighbor search
    tree = KDTree(grid_points)

    # count the number of occupied voxels
    occupied = np.zeros(len(grid_points), dtype=bool)
    for point, radius in zip(points, radii):
        occupied[tree.query_ball_point(point, radius, p=2)] = True
    return occupied.sum() * x_pitch * y_pitch * z_pitch
