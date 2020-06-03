#%%
# Setup
import numpy as np

F = (0, 0, 0)
G = (-2, -1, 1)
H = (1, -3, -2)

#%%
# Find normal vector to find equation of plane
FG = [g - f for f, g in zip(F, G)]
GH = [h - g for g, h in zip(G, H)]
n = np.cross(FG, GH)
print(f"normal vec: {n}")
# u = n / np.linalg.norm(n)
# print(f"unit normal vec: {u}")

# %%
# Determine if point Q is inside triangle FGH
P = (7.0 / 3, 11.0 / 3, 14.0 / 3)
FP = [p - f for f, p in zip(F, P)]
A_FPG = 0.5 * np.linalg.norm(np.cross(FG, FP))

GP = [p - g for g, p in zip(G, P)]
A_GPH = 0.5 * np.linalg.norm(np.cross(GH, GP))

HF = [f - h for h, f in zip(H, F)]
HP = [p - h for h, p in zip(H, P)]
A_HPF = 0.5 * np.linalg.norm(np.cross(HF, HP))

FH = [h - f for f, h in zip(F, H)]
A_FGH = 0.5 * np.linalg.norm(np.cross(FG, FH))

print(
    f"Area of triangle FGH: {A_FGH}\nSum of area of triangles FPG, GPH, HPF: {A_FPG+A_GPH+A_HPF}"
)

# %%
