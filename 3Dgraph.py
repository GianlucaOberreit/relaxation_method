import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# List of filenames
#n = [123, 100]
#filenames = ['10.npy', '25.npy', '50.npy', '100.npy']
filenames = ['21.npy', '100.npy']
#filenames = [f'test_{i}.npy' for i in n]

# Create a new figure
fig = plt.figure(figsize=(12, 8))

for i, filename in enumerate(filenames):
    # Load data from file
    data = np.load(filename)
    
    # Define the x, y coordinates
    x = np.linspace(0, 1, data.shape[0])
    y = np.linspace(0, 1, data.shape[1])
    X, Y = np.meshgrid(x, y)
    
    # Create a subplot
    ax = fig.add_subplot(2, 2, i + 1, projection='3d')
    
    # Plot the contour
    contour = ax.contour3D(X, Y, data, 500, cmap='plasma')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'L = {filename}')
    
    # Add a color bar
    fig.colorbar(contour, ax=ax, shrink=0.5, aspect=5)

plt.tight_layout()
plt.show()
