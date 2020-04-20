#%%
# Setup
import numpy as np

AB = (1, -1, -1)
CD = (3, 1, -3)

#%%
# Check parallel
cross_AB_CD = np.cross(AB, CD)

#%%
# Projection minimum distance (https://math.stackexchange.com/questions/302598/how-to-prove-that-two-lines-in-3d-are-not-parallel-and-do-not-intersect-also-h)
A = (1, 2, 3)
C = (1, 3, 4)
AC = [c - a for a, c in zip(A, C)]
min_dist = np.linalg.norm(np.dot(AC, cross_AB_CD)) / np.linalg.norm(cross_AB_CD)
print(f"Shortest distance between lines AB and CD: {min_dist}")


# %%
