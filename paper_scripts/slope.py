import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# Sample data points
x = np.linspace(0, 1, 20)  # start from 0.1 to avoid taking the third root of zero
y = np.cbrt(x)  # Third root of x

# Fit a spline to the data points. Here, s determines the smoothness, lower values fit the data more tightly
spline = UnivariateSpline(x, y, s=0.02)

# Get the first derivative of the spline
spline_derivative = spline.derivative()

# Find x values where the slope is approximately 1
# Creating a dense x range to evaluate and find when the derivative is closest to 1
x_dense = np.linspace(min(x), max(x), 300)
y_dense = spline(x_dense)
dy_dx = spline_derivative(x_dense)

# Finding the point closest to the desired slope of 1
closest_index = np.argmin(np.abs(dy_dx - 1))
x_at_slope_one = x_dense[closest_index]
y_at_slope_one = y_dense[closest_index]

print(f"The curve has a slope of approximately 1 at x = {x_at_slope_one:.2f}, y = {y_at_slope_one:.2f}")

# Plotting the results
plt.figure(figsize=(8, 6))
plt.plot(x, y, 'o', label='Data points')
plt.plot(x_dense, y_dense, label='Fitted spline')
plt.plot(x_dense, dy_dx, label='Derivative of spline')
plt.axvline(x=x_at_slope_one, color='r', linestyle='--', label=f'Slope = 1 at x = {x_at_slope_one:.2f}')
plt.title('Curve and its Derivative')
plt.xlabel('Parameter Value')
plt.ylabel('F1-Score')
plt.legend()
plt.grid(True)

# Set equal scaling
plt.axis('equal')  # This line ensures that the scale of the x and y axes are the same

plt.show()
