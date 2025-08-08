import matplotlib.pyplot as plt
import numpy as np

# Example data
x = np.linspace(-10, 10, 50)  # 50 points from -10 to 10 on the x-axis
y = np.random.normal(0, 5, 50)  # Random y values centered around 0

# Create the plot with a larger size for better visibility
plt.figure(figsize=(10, 6))

# Add background color for y > 0 area
plt.axhspan(0, max(y) + 1, facecolor='red', alpha=0.2)

# Add background color for y < 0 area
plt.axhspan(min(y) - 1, 0, facecolor='blue', alpha=0.2)

# Plot several layers of points with decreasing opacity
sizes = [120, 180, 240]  # Sizes for different layers
alphas = [0.3, 0.15, 0.05]  # Decreasing alpha for gradient effect

# Plot from greater to smaller for layering
for s, alpha in zip(sizes, alphas):
    plt.scatter(x, y, color='none', edgecolors='white', linewidths=20, s=s, alpha=alpha, zorder=3)

# Main scatter plot
plt.scatter(x, y, color='black', edgecolors='white', linewidths=3, s=100, zorder=5)

# Add a label and show the plot
plt.title('Scatter Plot with Gradient-like Transition')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)

plt.show()