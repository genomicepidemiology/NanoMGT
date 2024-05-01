import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])
    results = []

    for (mrd, batch), group in grouped:
        np_values = group['np'].values
        f1_scores = group['average_f1'].values

        # Normalize np_values and f1_scores to range [0, 1]
        np_values_normalized = (np_values - np.min(np_values)) / (np.max(np_values) - np.min(np_values))
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / (np.max(f1_scores) - np.min(f1_scores))

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(np_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense np values for a finer analysis in the normalized range
        np_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(np_dense_normalized)
        derivative_values_normalized = derivative(np_dense_normalized)

        # Calculate the target derivative for a 30-degree angle
        target_slope = np.tan(np.radians(20))

        # Collect all np values where the derivative is close to the 30-degree slope
        valid_np_values = []
        for idx in range(len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.02 and f1_dense_normalized[idx] > f1_dense_normalized[0]:
                valid_np_value = np_values.min() + np_dense_normalized[idx] * (np_values.max() - np_values.min())
                valid_np_values.append(valid_np_value)

        # Calculate the lowest gradient slope angle in degrees
        min_slope_angle = np.degrees(np.arctan(np.min(derivative_values_normalized)))

        # Get the first and last F1 score from the normalized data
        first_f1_score = f1_dense_normalized[0]
        last_f1_score = f1_dense_normalized[-1]

        # Append results
        results.append({
            'mrd': mrd,
            'batch': batch,
            'first_f1_score': first_f1_score,
            'last_f1_score': last_f1_score,
            'lowest_slope_angle': min_slope_angle,
            'valid_np_values': valid_np_values
        })

    return results

def main():
    csv_file = 'np_f1_scores.csv'
    results = process_data(csv_file)
    for result in results:
        print(f"MRD: {result['mrd']:.2f}, Batch: {result['batch']}, First F1 Score: {result['first_f1_score']:.2f}, Last F1 Score: {result['last_f1_score']:.2f}, Lowest Slope Angle: {result['lowest_slope_angle']:.2f} degrees")
        if result['valid_np_values']:
            print(f"   Valid np values at slopes close to 30 degrees: {result['valid_np_values']}")
        else:
            print("   No valid np values found close to 30 degrees.")

if __name__ == "__main__":
    main()
