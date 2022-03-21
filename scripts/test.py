import matplotlib.pyplot as plt
import build.hermite_splines as splines
import numpy as np


left = -2.0
center = 1.0
right = 2.0

print("initialization...")
grid = np.linspace(left, right, 200)


print("cheking phi function...")
phi_vals = [splines.phi(x, left, center, right) for x in grid]
phi1_vals = [splines.phi1(x, left, center, right) for x in grid]
phi2_vals = [splines.phi2(x, left, center, right) for x in grid]

print("phi is correct!")

print("cheking psi function...")
psi_vals = [splines.psi(x, left, center, right) for x in grid]
psi1_vals = [splines.psi1(x, left, center, right) for x in grid]
psi2_vals = [splines.psi2(x, left, center, right) for x in grid]

print("psi is correct!")


# plt.figure(1)
plt.plot(grid, phi_vals, grid, phi1_vals, grid, phi2_vals)
plt.legend(["phi", "1der", "2der"])
plt.grid(True)
plt.show()

# plt.figure(2)
plt.plot(grid, psi_vals, grid, psi1_vals, grid, psi2_vals)
plt.legend(["psi", "1der", "2der"])
plt.grid(True)
plt.show()
