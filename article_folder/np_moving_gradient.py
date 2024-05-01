import pandas as pd
import numpy
from scipy.interpolate import UnivariateSpline

def moving_average(data, window_size):
    """Calculate the moving average using a simple rolling window approach."""
    weights = numpy.ones(window_size) / window_size
    return numpy.convolve(data, weights, mode='valid')

def process_data(csv_file):
    data = pd.read_csv(csv_file)
    grouped = data.groupby(['mrd', 'batch'])

    for (mrd, batch), group in grouped:
        novelty_penalty_values = group['np'].values  # Renamed variable to avoid confusion
        f1_scores = group['average_f1'].values

        # Fit a spline to the data points
        spline = UnivariateSpline(novelty_penalty_values, f1_scores, s=None)
        derivative = spline.derivative()

        # Generate dense novelty_penalty values for a finer analysis
        novelty_penalty_dense = numpy.linspace(novelty_penalty_values.min(), novelty_penalty_values.max(), 300)
        derivative_values = derivative(novelty_penalty_dense)

        # Calculate moving average of the derivative
        window_size = 5
        ma_derivative_values = moving_average(derivative_values, window_size)

        # Adjust novelty_penalty_dense to match the length of the moving average output
        ma_novelty_penalty_values = novelty_penalty_dense[len(novelty_penalty_dense) - len(ma_derivative_values):]

        # Initialize to track the best f1-score and associated novelty_penalty value
        max_f1 = f1_scores[0]
        best_novelty_penalty = novelty_penalty_values[0]
        significant_improvement_threshold = 0.05

        for novelty_penalty, f1 in zip(ma_novelty_penalty_values, numpy.interp(ma_novelty_penalty_values, novelty_penalty_values, f1_scores)):
            if f1 > max_f1 * (1 + significant_improvement_threshold):
                max_f1 = f1
                best_novelty_penalty = novelty_penalty

        print(f"For MRD: {mrd:.2f}, Batch: {batch}, the novelty_penalty value with the highest significant observed F1-score is {best_novelty_penalty:.2f}")

def main():
    csv_file = 'np_f1_scores.csv'
    process_data(csv_file)

if __name__ == "__main__":
    main()
