import numpy as np

def forward_iteration_matrix_2d_3x3():
    N = 9  # 3x3 grid
    forward_matrix = np.zeros((N, N))
    idx = np.arange(N).reshape((3, 3), order='F')

    for i in range(3):
        for j in range(3):
            current_index = idx[i, j]

            if i + 1 < 3:
                forward_matrix[current_index, idx[i + 1, j]] = 1

            if j + 1 < 3:
                forward_matrix[current_index, idx[i, j + 1]] = 1

    return forward_matrix

# Example usage
forward_matrix_2d_3x3 = forward_iteration_matrix_2d_3x3()

# Print the 3x3 forward iteration matrix
print("Forward Iteration Matrix (2D - 3x3):")
print(forward_matrix_2d_3x3)

# https://github.com/KarelVanHoey/MD_PistonRing