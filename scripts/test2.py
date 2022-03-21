"""
int CHermiteBC::locate(double x)
{
    Eigen::VectorXd& grid = space.grid;
    int right = grid.size() - 1;
    int left = 0;
    int mid = left + (right - left) / 2;

    while (left < right){
        if (grid[mid] < x) {
            left = mid;
        } else if (grid[mid] > x) {
            right = mid;
        }
    }

    return right;
}
"""
import numpy as np


def binary_search(x: float, grid):
    right = len(grid) - 1
    left = 0
    mid = left + (right - left) // 2

    i = 0
    while right - left > 1:
        print(f"left: {left}, right: {right}, mid: {mid}")
        if grid[mid] < x:
            left = mid
        elif grid[mid] > x:
            right = mid
        else:
            return mid

        mid = left + (right - left) // 2

    return left


if __name__ == "__main__":
    grid = [0.1, 0.5, 0.65, 10, 15]
    ns = list(range(0, len(grid)))
    dct = {n: g for n, g in zip(ns, grid)}
    print(dct)
    x = 16
    print(x)
    out = binary_search(x, grid)
    print(out)
