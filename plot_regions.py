import numpy as np
import matplotlib.pyplot as plt

# 1. Create a shared grid coordinate system
x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
X, Y = np.meshgrid(x, y)

# 2. Define two different 2D datasets
Z1 = np.sin(X**2 + Y**2)
Z2 = np.exp(-((X-1)**2 + (Y-1)**2))

fig, ax = plt.subplots()

# 3. Plot the background layer (fully opaque)
mesh1 = ax.pcolormesh(X, Y, Z1, cmap='Blues', shading='auto')

# 4. Overlay the foreground layer (semi-transparent)
mesh2 = ax.pcolormesh(X, Y, Z2, cmap='Reds', alpha=0.5, shading='auto')

plt.colorbar(mesh1, ax=ax, label='Layer 1')
plt.colorbar(mesh2, ax=ax, label='Layer 2')
plt.show()