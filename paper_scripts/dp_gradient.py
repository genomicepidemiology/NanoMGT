import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])
    results = []

    for (mrd, batch), group in grouped:
        dp_values = group['dp'].values  # Adjusted to use 'dp' values
        f1_scores = group['average_f1'].values

        # Normalize dp_values and f1_scores to range [0, 1]
        dp_values_normalized = (dp_values - np.min(dp_values)) / (np.max(dp_values) - np.min(dp_values))
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / (np.max(f1_scores) - np.min(f1_scores))

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(dp_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense dp values for a finer analysis in the normalized range
        dp_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(dp_dense_normalized)
        derivative_values_normalized = derivative(dp_dense_normalized)

        # Calculate the target derivative for a 30-degree angle
        target_slope = np.tan(np.radians(30))

        # Initialize the last valid index
        last_valid_index = 0

        # Iterate over all points to find derivatives close to the target slope
        for idx in range(len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.01:  # Close to 30-degree slope
                if f1_dense_normalized[idx] > f1_dense_normalized[last_valid_index]:  # Ensure F1-score is increasing
                    last_valid_index = idx

        # Get the dp value at the last valid index found
        target_dp_value = dp_values.min() + dp_dense_normalized[last_valid_index] * (dp_values.max() - dp_values.min())
        results.append((mrd, batch, target_dp_value))

    return results

def main():
    csv_file = 'dp_f1_scores.csv'  # Adjusted to use the correct CSV file for 'dp' values
    results = process_data(csv_file)
    for result in results:
        mrd, batch, peak_value = result
        print(f"The dp value with a 30-degree slope for MRD: {mrd:.2f}, Batch: {batch} is at dp: {peak_value}")

if __name__ == "__main__":
    main()
