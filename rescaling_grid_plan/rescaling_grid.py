import argparse
import logging
import csv
import random
import pydicom
import copy

logger = logging.getLogger(__name__)


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
    parser.add_argument('input', help='Path to DICOM file to be modified')
    parser.add_argument('-o', '--output', default="output.dcm", required=False, help='Path to output DICOM file')
    parser.add_argument('-w', '--weights', required=False, help='Optional path CSV file with weight per energy layer',
                        default=None)
    parser.add_argument('-p', '--print', type=int, default=None, help='Number of random values to print for comparison')
    parser.add_argument('-v', '--verbosity', action='count', help="increase verbosity", default=0)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--scale_factor', type=float, default=1.0, help='Scale factor to be applied everywhere')
    group.add_argument('-d', '--dose', type=float, default=None, help='New target dose')

    args = parser.parse_args()

    if args.verbosity == 1:
        logging.basicConfig(level=logging.INFO)
    elif args.verbosity > 1:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig()

    dcm = pydicom.dcmread(args.input)
    dcm_new = copy.deepcopy(dcm)

    weights = None
    if args.weights:
        weights = read_weights_from_csv(args.weights)

    if (args.dose):
        beam_dose = dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamDose
        scale_factor = args.dose / beam_dose
        logger.debug(f"Beam dose: {beam_dose} Gy")
        logger.debug(f"Target dose: {args.dose} Gy")
        logger.debug(f"Scale factor: {scale_factor}")
    else:
        scale_factor = args.scale_factor

    logger.debug(f"Scale Factor: {scale_factor}")
    # final_original_cumulative_weight = dcm.IonBeamSequence[0].FinalCumulativeMetersetWeight
    # number_of_control_points = dcm.IonBeamSequence[0].NumberOfControlPoints
    # number_of_energy_layers = int(number_of_control_points / 2)

    ion_control_point_sequence = dcm_new.IonBeamSequence[0].IonControlPointSequence  # all double energy layers

    original_cumulative_weight = 0.0  # per energy layer
    new_cumulative_weight = 0.0  # per energy layer

    for i, icp in enumerate(ion_control_point_sequence):

        logger.debug(f" --------- Processing energy layer {i}")

        icp.CumulativeMetersetWeight = new_cumulative_weight

        weights = icp.ScanSpotMetersetWeights
        new_weights = [w * scale_factor for w in weights]
        icp.ScanSpotMetersetWeights = new_weights

        original_cumulative_weight += sum(weights)
        new_cumulative_weight += sum(new_weights)  # Calculate the new cumulative weight

        logger.debug(f"Layer {i} Cumulative Weight Before: {original_cumulative_weight}")
        logger.debug(f"Layer {i} Cumulative Weight After: {new_cumulative_weight}")

        if args.print:
            print_comparison(i, scale_factor, weights, new_weights, args.print)

    dcm_new.IonBeamSequence[0].FinalCumulativeMetersetWeight = new_cumulative_weight
    logger.info(f"Final Cumulative Weight before rescaling: {original_cumulative_weight}")
    logger.info(f"Final Cumulative Weight after rescaling: {new_cumulative_weight}")

    dcm_new.save_as(args.output)
    logger.debug(f"Weight rescaled plan is saved as {args.output}")


if __name__ == "__main__":
    rescale()
