import numpy as np
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
import matplotlib
from brokenaxes import brokenaxes

# Example data (replace these with your actual data)
data_age = np.random.randint(20, 70, 100)
data_y_name = np.random.normal(50, 30, 100)
x_fit = np.linspace(20, 70, 100)
y_fit = np.polyval(np.polyfit(data_age, data_y_name, 1), x_fit)
residuals = data_y_name - np.polyval(np.polyfit(data_age, data_y_name, 1), data_age)
r_squared = 0.85  # Example R-squared value

fig = plt.figure(figsize=(8, 6))  # Create a figure
# Create broken axes by specifying the ranges for the y-axis
bax = brokenaxes(ylims=((10, 30), (60, 80)), hspace=0.05)
bax.set_yscale('log')  # 设置为对数坐标轴

# Plot the first dataset on broken axes
bax.scatter(data_age, data_y_name, s=10, label='Data Points 1', color='blue')
bax.plot(x_fit, y_fit, color='red', label=f'R-squared 1: {r_squared:.3f}')
bax.fill_between(x_fit, y_fit - 0.5 * np.std(residuals), y_fit + 0.5 * np.std(residuals), color='gray', alpha=0.2)

# Plot the second dataset on the same broken axes
bax.scatter(data_age, data_y_name, s=10, label='Data Points 2', color='green', marker='x')
bax.plot(x_fit, y_fit, color='orange', label=f'R-squared 2: {r_squared:.3f}', linestyle='--')
bax.fill_between(x_fit, y_fit - 0.5 * np.std(residuals), y_fit + 0.5 * np.std(residuals), color='yellow', alpha=0.2)

# Common labels, legend, and title
bax.set_xlabel('Age (year)')
bax.set_ylabel('y_name')
bax.legend()
# Optionally: tkinter.mainloop() ax.set_title(f'{y_name} vs Age with Polynomial Fit', fontsize=20)

plt.tight_layout()  # Adjust the layout to prevent overlapping labels
plt.show()  # Display the figure