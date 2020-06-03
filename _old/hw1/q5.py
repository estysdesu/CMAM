#%%
# Setup
import numpy as np

A = (0, 0, 0, 1)
B = (2, 0, 0, 1)
C = (1, 2, 0, 1)

T1 = [[2, 0, 0, 0], [0, 3, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
T2 = [[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
T3 = [
    [np.cos(np.deg2rad(90)), np.sin(np.deg2rad(90)), 0, 0],
    [-np.sin(np.deg2rad(90)), np.cos(np.deg2rad(90)), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
]
T4 = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [1, 2, 0, 1]]

np.matrix([A, B, C]) * np.matrix(T1) * np.matrix(T2) * np.matrix(T3) * np.matrix(T4)


# %%
