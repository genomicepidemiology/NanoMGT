import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline

def moving_average(data, window_size):
    """Calculate the moving average using a simple rolling window approach."""
    weights = np.ones(window_size) / window_size
    return np.convolve(data, weights, mode='valid')

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['maf', 'batch'])

    for (maf, batch), group in grouped:
        pp_values = group['pp'].values
        f1_scores = group['average_f1'].values

        # Fit a spline to the data points
        spline = UnivariateSpline(pp_values, f1_scores, s=None)  # Let the spline choose the best `s`
        derivative = spline.derivative()

        # Generate dense pp values for a finer analysis
        pp_dense = np.linspace(pp_values.min(), pp_values.max(), 300)
        derivative_values = derivative(pp_dense)

        # Calculate moving average of the derivative
        window_size = 5  # Adjust window size as needed
        ma_derivative_values = moving_average(derivative_values, window_size)

        # Adjust pp_dense to match the length of the moving average output
        ma_pp_values = pp_dense[len(pp_dense) - len(ma_derivative_values):]

        # Initialize to track the best f1-score and associated pp value
        max_f1 = f1_scores[0]  # Start with the first f1-score
        best_pp = pp_values[0]  # Start with the first pp value
        significant_improvement_threshold = 0.05  # 5% improvement considered significant

        for pp, f1 in zip(ma_pp_values, np.interp(ma_pp_values, pp_values, f1_scores)):
            if f1 > max_f1 * (1 + significant_improvement_threshold):  # Check if there is a significant improvement in f1-score
                max_f1 = f1  # Update the maximum f1-score
                best_pp = pp  # Update the best pp

        print(f"For MAF: {maf:.2f}, Batch: {batch}, the pp value with the highest significant observed F1-score is {best_pp:.2f}")

def main():
    csv_file = 'pp_f1_scores.csv'
    process_data(csv_file)

if __name__ == "__main__":
    main()
