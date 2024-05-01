import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])
    results = []

    for (mrd, batch), group in grouped:
        pp_values = group['pp'].values  # Changed from cor_values to pp_values
        f1_scores = group['average_f1'].values

        # Normalize pp_values and f1_scores to range [0, 1]
        pp_values_normalized = (pp_values - np.min(pp_values)) / (np.max(pp_values) - np.min(pp_values))
        f1_scores_normalized = (f1_scores - np.min(f1_scores)) / (np.max(f1_scores) - np.min(f1_scores))

        # Fit a spline to the normalized data points
        spline = UnivariateSpline(pp_values_normalized, f1_scores_normalized, s=None)
        derivative = spline.derivative()

        # Generate dense pp values for a finer analysis in the normalized range
        pp_dense_normalized = np.linspace(0, 1, 450)
        f1_dense_normalized = spline(pp_dense_normalized)
        derivative_values_normalized = derivative(pp_dense_normalized)

        # Calculate the target derivative for a 20-degree angle
        target_slope = np.tan(np.radians(20))

        # Collect all pp values where the derivative is close to the 20-degree slope
        valid_pp_values = []
        for idx in range(len(derivative_values_normalized)):
            if abs(derivative_values_normalized[idx] - target_slope) < 0.02 and f1_dense_normalized[idx] > f1_dense_normalized[0]:
                valid_pp_value = pp_values.min() + pp_dense_normalized[idx] * (pp_values.max() - pp_values.min())
                valid_pp_values.append(valid_pp_value)

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
            'valid_pp_values': valid_pp_values
        })

    return results

def main():
    csv_file = 'pp_f1_scores.csv'  # Changed to use the correct CSV file for 'pp' values
    results = process_data(csv_file)
    for result in results:
        print(f"MRD: {result['mrd']:.2f}, Batch: {result['batch']}, First F1 Score: {result['first_f1_score']:.2f}, Last F1 Score: {result['last_f1_score']:.2f}, Lowest Slope Angle: {result['lowest_slope_angle']:.2f} degrees")
        if result['valid_pp_values']:
            print(f"   Valid pp values at slopes close to 20 degrees: {result['valid_pp_values']}")
        else:
            print("   No valid pp values found close to 20 degrees.")

if __name__ == "__main__":
    main()
