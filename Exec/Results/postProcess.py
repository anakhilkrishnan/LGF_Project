import yt
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import exp1
import os

# ==========================================
# User Defined Inputs (Matching ParmParse)
# ==========================================
n_cell = 128
dom_lo = -1.0
dom_hi = 1.0
variance = 1.0   # From SourceField.H
plotfile = "plt00000"

# Ensure scipy is installed: pip install scipy
# ==========================================

# 1. Load AMReX Data
print(f"Loading {plotfile}...")
ds = yt.load(plotfile)

# Extract the full numerical grid. 
# ds.domain_dimensions is [128, 128, 1], squeezing it removes the dummy Z-axis
grid = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
phi_num = grid['phi'].squeeze().v  # .v strips the yt units to give a raw numpy array

# 2. Generate Physical Coordinate Grid (Cell-Centered)
# AMReX data is evaluated at the center of each cell.
dx = (dom_hi - dom_lo) / n_cell
x_coords = np.linspace(dom_lo + dx/2.0, dom_hi - dx/2.0, n_cell)
y_coords = np.linspace(dom_lo + dx/2.0, dom_hi - dx/2.0, n_cell)

# indexing='ij' ensures X and Y arrays match yt's [x_index, y_index] memory layout
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')

# 3. Compute Exact Analytical Solution
R2 = X**2 + Y**2

# Apply a tiny offset to avoid log(0) and exp1(0) if a cell center hits exactly (0,0)
R2 = np.where(R2 == 0, 1e-15, R2)

# 2D Exact Solution: phi(r) = (v/4) * [ ln(r^2/v) + E1(r^2/v) ]
# Note: Because the LGF convolution acts on the infinite domain, the arbitrary 
# integration constant C naturally vanishes.
phi_exact = (variance / 4.0) * (np.log(R2 / variance) + exp1(R2 / variance))

# 4. Compute Absolute Error
abs_error = np.abs(phi_num - phi_exact)
max_error = np.max(abs_error)
print(f"Maximum Absolute Error: {max_error:.4e}")

# ==========================================
# 5. Plotting Results
# ==========================================
# We use .T (transpose) on the arrays because yt stores them as [x, y], 
# but plt.imshow expects image data as [y, x].

extent = [dom_lo, dom_hi, dom_lo, dom_hi]
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))

# Plot A: Numerical
im0 = axs[0].imshow(phi_num.T, extent=extent, origin='lower', cmap='viridis')
axs[0].set_title(r"AMReX Numerical $\phi$")
fig.colorbar(im0, ax=axs[0])

# Plot B: Exact Analytical
im1 = axs[1].imshow(phi_exact.T, extent=extent, origin='lower', cmap='viridis')
axs[1].set_title(r"Exact Analytical $\phi$")
fig.colorbar(im1, ax=axs[1])

# Plot C: Error Map
im2 = axs[2].imshow(abs_error.T, extent=extent, origin='lower', cmap='magma')
axs[2].set_title(f"Absolute Error\n(Max: {max_error:.2e})")
fig.colorbar(im2, ax=axs[2])

plt.tight_layout()
plt.savefig("error_maps_2D.png", dpi=300)
print(f"Saved 2D visual error map to: {os.path.abspath('error_maps_2D.png')}")

# ==========================================
# 6. Optional: 1D Cross-Section Plot
# ==========================================
mid_idx = n_cell // 2  # The index corresponding to y = 0
plt.figure(figsize=(8, 5))

# Plotting the horizontal centerline
plt.plot(x_coords, phi_num[:, mid_idx], 'ro', label='AMReX Numerical', markersize=5)
plt.plot(x_coords, phi_exact[:, mid_idx], 'k-', label='Exact Analytical', linewidth=2)

plt.title("1D Cross-Section at $y = 0$")
plt.xlabel("$x$ coordinate")
plt.ylabel(r"Potential ($\phi$)")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig("slice_1D.png", dpi=300)
print(f"Saved 1D slice comparison to: {os.path.abspath('slice_1D.png')}")