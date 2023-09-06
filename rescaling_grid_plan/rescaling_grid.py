import argparse
import csv
import random
import pydicom
import copy

def read_weights_from_csv(csv_file_path):
    weights = []
    with open(csv_file_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            weights.append(float(row[0]))
    return weights

def print_comparison(layer, scale_factor, original_weights, modified_weights, num_values):
    print(f"\nLayer {layer}  Scale Factor {scale_factor}")
    print("Original | Modified")
    print("---------|---------")

    if len(original_weights) >= num_values:
        sample_indices = random.sample(range(len(original_weights)), num_values)
    else:
        sample_indices = range(len(original_weights))

    sampled_original = [original_weights[i] for i in sample_indices]
    sampled_new = [modified_weights[i] for i in sample_indices]

    for original, modified in zip(sampled_original, sampled_new):
        print(f"{original:8.4f} | {modified:8.4f}")

def rescale():
    parser = argparse.ArgumentParser(description='Modify DICOM file weights.')
    parser.add_argument('-i', '--input', required=True, help='Path to input DICOM file')
    parser.add_argument('-w', '--weights', required=True, help='Path to weights CSV file')
    parser.add_argument('-o', '--output', required=True, help='Path to output DICOM file')
    parser.add_argument('-p', '--print', type=int, default=None, help='Number of random values to print for comparison')
    parser.add_argument('-pd', '--plan_dose', type=float, default=None, help='Plan dose')
    parser.add_argument('-rd', '--rescale_dose', type=float, default=None, help='Rescaled dose')
    args = parser.parse_args()

    dicom_data = pydicom.dcmread(args.input)
    new_dicom_data = copy.deepcopy(dicom_data)
    scale_factors = read_weights_from_csv(args.weights)
    ion_control_point_sequence = new_dicom_data.IonBeamSequence[0].IonControlPointSequence

    plan_rescale_ratio = 1.0  # Initialize to 1 so it doesn't affect multiplication if not set
    if args.plan_dose is not None and args.rescale_dose is not None:
        plan_rescale_ratio = args.rescale_dose / args.plan_dose

    for i, scale_factor in zip(range(0, len(ion_control_point_sequence), 2), scale_factors):
        weights = ion_control_point_sequence[i].ScanSpotMetersetWeights
        # Modify the weights with both scale_factor and plan_rescale_ratio
        new_weights = [w * scale_factor * plan_rescale_ratio for w in weights]
        ion_control_point_sequence[i].ScanSpotMetersetWeights = new_weights

        if args.print:
            print_comparison(i // 2 + 1, scale_factor, weights, new_weights, args.print)

    new_dicom_data.save_as(args.output)
    print(f"Weight rescaled paln is saved as {args.output}")

if __name__ == "__main__":
    rescale()
