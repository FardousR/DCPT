import sys
import argparse
import logging
import copy
import random
import pydicom
from datetime import datetime


logger = logging.getLogger(__name__)

MU_MIN = 1.0  # at least this many MU in a single spot


def read_weights(csv_file_path):
    weights = []
    with open(csv_file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            weights.append(float(line))
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


def main(args=None):
    parser = argparse.ArgumentParser(description='Modify DICOM file weights.')
    parser.add_argument('input', help='Path to DICOM file to be modified')
    parser.add_argument('-o', '--output', default="output.dcm", required=False, help='Path to output DICOM file')
    parser.add_argument('-l', '--label', default=None, required=False, help='Label string for DICOM file')
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

    if args.label:
        dcm_new.RTPlanLabel = args.label
    now = datetime.now()
    dcm_new.RTPlanDate = now.strftime("%Y%m%d")
    dcm_new.RTPlanTime = now.strftime("%H%M%S")

    # not sure if beam dose is always given in dicom?
    original_beam_dose = dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamDose
    if (args.dose):
        scale_factor = args.dose / original_beam_dose
    else:
        scale_factor = args.scale_factor

    new_beam_dose = original_beam_dose * scale_factor

    final_original_cumulative_weight = dcm.IonBeamSequence[0].FinalCumulativeMetersetWeight
    number_of_control_points = dcm.IonBeamSequence[0].NumberOfControlPoints
    number_of_energy_layers = int(number_of_control_points / 2)

    csv_weights = None
    if args.weights:
        # csv_weights = read_weights_from_csv(args.weights)
        csv_weights = read_weights(args.weights)
        csv_weigths_len = len(csv_weights)
        if csv_weigths_len != number_of_energy_layers:
            raise Exception(f"CSV file energy layers {csv_weigths_len} must \
                            match number of energy layers in dicom file {number_of_energy_layers}.")

    original_beam_meterset = dcm.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset
    new_beam_meterset = original_beam_meterset * scale_factor
    new_meterset_per_weight = new_beam_meterset / final_original_cumulative_weight

    ion_control_point_sequence = dcm_new.IonBeamSequence[0].IonControlPointSequence  # all double energy layers

    original_cumulative_weight = 0.0  # per energy layer
    new_cumulative_weight = 0.0  # per energy layer
    points_discarded = 0

    for i, icp in enumerate(ion_control_point_sequence):

        logger.debug(f" --------- Processing energy layer {i}")

        icp.CumulativeMetersetWeight = new_cumulative_weight

        weights = icp.ScanSpotMetersetWeights
        if csv_weights:
            _ie = int(i * 0.5)
            csv_weight = csv_weights[_ie]
            logger.debug(f"CSV weight energy layer {_ie} {csv_weight}")
        else:
            csv_weight = 1.0

        new_weights = [0.0] * len(weights)

        for j, w in enumerate(weights):
            value = w * csv_weight
            if value > 0.0 and value * new_meterset_per_weight < MU_MIN:
                logger.debug(f"Discarding point with weight {value:.2f} and {value*new_meterset_per_weight:.2f} [MU]")
                points_discarded += 1
                value = 0.0
            new_weights[j] = value

        icp.ScanSpotMetersetWeights = new_weights

        original_cumulative_weight += sum(weights)
        new_cumulative_weight += sum(new_weights)  # Calculate the new cumulative weight

        logger.debug(f"Layer {i} Cumulative Weight Before: {original_cumulative_weight}")
        logger.debug(f"Layer {i} Cumulative Weight After: {new_cumulative_weight}")

        if args.print:
            print_comparison(i, scale_factor, weights, new_weights, args.print)

    # set remaining meta data
    dcm_new.IonBeamSequence[0].FinalCumulativeMetersetWeight = new_cumulative_weight
    dcm_new.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset = new_beam_meterset
    dcm_new.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamDose = new_beam_dose

    dcm_new.save_as(args.output)
    hline = 72 * '-'
    logger.info("                                  Original           New   ")
    logger.info(hline)
    logger.info(f"Final Cumulative Weight   : {original_cumulative_weight:14.2f}  {new_cumulative_weight:14.2f}  ")
    logger.info(f"Beam Meterset             : {original_beam_meterset:14.2f}  {new_beam_meterset:14.2f}  [MU] ")
    logger.info(f"Beam Dose                 : {original_beam_dose:14.2f}  {new_beam_dose:14.2f}  [Gy(RBE)]] ")
    logger.info(hline)
    logger.info(f"Scale Factor : {scale_factor:.4f}")
    logger.info(f"Weight rescaled plan is saved as : '{args.output}'")
    if points_discarded > 0:
        logger.warning(f" *** Discarded {points_discarded} spots which were below {MU_MIN:.2f} [MU] ***")


if __name__ == '__main__':

    sys.exit(main(sys.argv[1:]))
