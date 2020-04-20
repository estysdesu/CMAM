#%%
# Setup
import numpy as np

B = (1, 0, -1)
C = (-4, -2, 1)
D = (-2, 5, 3)

#%%
# Calculate area of traingle BCD
BC = [c - b for b, c in zip(B, C)]
BD = [d - b for b, d in zip(B, D)]
A_BCD = 0.5 * np.linalg.norm(np.cross(BC, BD))
print(f"Area of triangle BCD: {A_BCD}")

#%%
# Find plane of triangle BCD
n = np.cross(BC, BD)
right = np.dot(n, B)

# %%
# Find height of pyramid
P = (3, -2, 5)
H = np.abs(n[0] * P[0] + n[1] * P[1] + n[2] * P[2] - right) / np.sqrt(
    np.square(n[0]) + np.square(n[1]) + np.square(n[2])
)

# %%
