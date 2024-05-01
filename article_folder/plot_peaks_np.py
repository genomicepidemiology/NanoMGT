import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt


def process_data(csv_file, threshold=5):  # threshold in percentage
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])

    for (mrd, batch), group in grouped:
        if str(mrd) == "0.04" and batch == 9:
            np_values = group['np'].values
            f1_scores = group['average_f1'].values

            # Fit a spline to the data points
            spline = UnivariateSpline(np_values, f1_scores, s=None)  # Let the spline choose the best `s`
            derivative = spline.derivative()

            # Generate dense np values for a finer analysis
            np_dense = np.linspace(np_values.min(), np_values.max(), 300)
            f1_dense = spline(np_dense)  # Values of the spline
            derivative_values = derivative(np_dense)  # Values of the derivative

            # Calculate the percentage change in derivative values
            pct_change = np.abs(np.diff(derivative_values) / derivative_values[:-1]) * 100

            # Find the index where the percentage change drops below the threshold
            significant_change = np.where(pct_change < threshold)[0]
            if significant_change.size > 0:
                diminishing_point = significant_change[0] + 1  # +1 because diff reduces the array size by 1
            else:
                diminishing_point = len(np_dense) - 1  # No diminishing point found within threshold, take the last

            # Plotting
            plt.figure(figsize=(12, 6))
            plt.subplot(1, 2, 1)
            plt.plot(np_values, f1_scores, 'bo', label='Original Data')
            plt.plot(np_dense, f1_dense, 'r-', label='Spline Fit')
            plt.axvline(x=np_dense[diminishing_point], color='g', label='Diminishing Returns Start')
            plt.xlabel('np Values')
            plt.ylabel('F1 Scores')
            plt.title('F1 Scores and Spline Fit')
            plt.legend()

            plt.subplot(1, 2, 2)
            plt.plot(np_dense[:-1], pct_change, 'k-', label='Percentage Change')
            plt.axvline(x=np_dense[diminishing_point], color='g', label='Diminishing Returns Start')
            plt.xlabel('np Values')
            plt.ylabel('Percentage Change in Derivative')
            plt.title('Derivative Changes')
            plt.legend()

            plt.tight_layout()
            plt.show()

            return np_dense[diminishing_point]  # Return the np value of diminishing returns start


def main():
    csv_file = 'np_f1_scores.csv'
    threshold = 5  # Percentage change threshold
    diminishing_np_value = process_data(csv_file, threshold)
    print("Diminishing returns start at np value:", diminishing_np_value)


if __name__ == "__main__":
    main()
