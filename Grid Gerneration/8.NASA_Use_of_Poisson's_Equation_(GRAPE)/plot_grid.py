#!/usr/bin/env python3
"""
Visualization of Periodic Hill structured grid.
Plots all xi-lines and eta-lines to match Fröhlich et al. (2005) Figure 1.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Load grid
data = np.load("periodic_hill_grid.npz")
xg = data['xg']
yg = data['yg']
Nx = int(data['Nx'])
Ny = int(data['Ny'])

print(f"Grid loaded: Nx={Nx}, Ny={Ny}")
print(f"x range: [{xg.min():.2f}, {xg.max():.2f}]")
print(f"y range: [{yg.min():.2f}, {yg.max():.2f}]")

# ========== Figure 1: Full grid (all xi and eta lines) ==========
fig, ax = plt.subplots(1, 1, figsize=(18, 5))

# Plot xi-lines (constant j, varying i) — these are "horizontal" lines
for j in range(Ny):
    ax.plot(xg[:, j], yg[:, j], 'b-', linewidth=0.3, alpha=0.7)

# Plot eta-lines (constant i, varying j) — these are "vertical" lines
for i in range(Nx):
    ax.plot(xg[i, :], yg[i, :], 'r-', linewidth=0.3, alpha=0.7)

# Fill the hill region
x_bot = xg[:, 0]
y_bot = yg[:, 0]
ax.fill_between(x_bot, 0, y_bot, color='0.85', zorder=2)
ax.plot(x_bot, y_bot, 'k-', linewidth=1.0, zorder=3)
# Top wall
ax.plot([0, 252], [85.008, 85.008], 'k-', linewidth=1.0)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Periodic Hill Grid (TTM Poisson, Nx=%d, Ny=%d)' % (Nx, Ny))
ax.set_aspect('equal')
ax.set_xlim(-5, 257)
ax.set_ylim(-2, 90)

plt.tight_layout()
plt.savefig('periodic_hill_grid.png', dpi=200, bbox_inches='tight')
print("Saved: periodic_hill_grid.png")

# ========== Figure 2: Zoom near first hill ==========
fig2, ax2 = plt.subplots(1, 1, figsize=(12, 6))

for j in range(Ny):
    ax2.plot(xg[:, j], yg[:, j], 'b-', linewidth=0.4, alpha=0.7)
for i in range(Nx):
    ax2.plot(xg[i, :], yg[i, :], 'r-', linewidth=0.4, alpha=0.7)

ax2.fill_between(x_bot, 0, y_bot, color='0.85', zorder=2)
ax2.plot(x_bot, y_bot, 'k-', linewidth=1.2, zorder=3)
ax2.plot([0, 252], [85.008, 85.008], 'k-', linewidth=1.0)

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_title('Periodic Hill Grid - Zoom Near Hill')
ax2.set_aspect('equal')
ax2.set_xlim(-5, 130)
ax2.set_ylim(-2, 50)

plt.tight_layout()
plt.savefig('periodic_hill_grid_zoom.png', dpi=200, bbox_inches='tight')
print("Saved: periodic_hill_grid_zoom.png")

# ========== Figure 3: Zoom into hill crest ==========
fig3, ax3 = plt.subplots(1, 1, figsize=(10, 8))

for j in range(Ny):
    ax3.plot(xg[:, j], yg[:, j], 'b-', linewidth=0.5, alpha=0.8)
for i in range(Nx):
    ax3.plot(xg[i, :], yg[i, :], 'r-', linewidth=0.5, alpha=0.8)

ax3.fill_between(x_bot, 0, y_bot, color='0.85', zorder=2)
ax3.plot(x_bot, y_bot, 'k-', linewidth=1.5, zorder=3)

ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_title('Periodic Hill Grid - Hill Crest Detail')
ax3.set_aspect('equal')
ax3.set_xlim(-5, 70)
ax3.set_ylim(-2, 55)

plt.tight_layout()
plt.savefig('periodic_hill_grid_crest.png', dpi=200, bbox_inches='tight')
print("Saved: periodic_hill_grid_crest.png")

# ========== Verification: Check eta-line straightness ==========
print()
print("=== Eta-line Straightness Verification ===")
# Eta-line at i=0 (x=0): should be straight (vertical)
x_spread_left = xg[0, :].max() - xg[0, :].min()
print(f"  i=0   (x=0):   x-spread = {x_spread_left:.6e}  {'STRAIGHT' if x_spread_left < 1e-10 else 'CURVED'}")

# Eta-line at i=Nx-1 (x=Lx): should be straight (vertical)
x_spread_right = xg[-1, :].max() - xg[-1, :].min()
print(f"  i={Nx-1} (x=Lx): x-spread = {x_spread_right:.6e}  {'STRAIGHT' if x_spread_right < 1e-10 else 'CURVED'}")

# Sample some interior eta-lines — these should be curved
for idx in [10, 40, 80, 120, 150]:
    x_spread = xg[idx, :].max() - xg[idx, :].min()
    mean_x = xg[idx, :].mean()
    print(f"  i={idx:3d} (mean x≈{mean_x:.1f}): x-spread = {x_spread:.4f}  {'STRAIGHT' if x_spread < 0.01 else 'CURVED'}")

print()
print("All plots saved!")
