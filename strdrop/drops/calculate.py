import json
import os
import Levenshtein
import numpy
import pandas as pd

from cyvcf2 import VCF, Variant
from pathlib import Path
from typing import List, Tuple

import logging

logger = logging.getLogger(__name__)


def get_allele(variant: Variant, pos: int, step: int) -> str | None:
    """Return allele for variant sample number pos and allele step (0, 1)"""
    if variant.genotypes[pos][step] < 1:
        return variant.REF

    idx = variant.genotypes[pos][step] - 1
    return variant.REF if idx > len(variant.ALT) - 1 else variant.ALT[idx]


def get_variant_edr_sd(variant: Variant, ind_nr: int = 0) -> Tuple[float, float]:
    """Return edit ratio and sequencing depth for variant, individual"""
    a1 = get_allele(variant, ind_nr, 0)
    a2 = get_allele(variant, ind_nr, 1)

    edit_ratio = Levenshtein.ratio(a1, a2)

    ref_sd = 0
    alt_sd = 0
    sd_values = variant.format("SD")[ind_nr]
    if len(sd_values) == 1:
        alt_sd = max(int(sd_values[0]), 0)
    if len(sd_values) == 2:
        alt_sd = max(int(sd_values[1]), 0)
        ref_sd = max(int(sd_values[0]), 0)

    sd = ref_sd + alt_sd

    return (edit_ratio, sd)


def parse_sds_training(
    file: Path, sequencing_depths: dict = None, edit_ratios: dict = None, chrom: dict = None
) -> bool:
    """Parse SDs from VCF. Return None if file was not found. Return count of individuals if found."""
    if not os.path.isfile(file):
        return False

    if chrom is None:
        chrom = {}
    if edit_ratios is None:
        edit_ratios = {}
    if sequencing_depths is None:
        sequencing_depths = {}

    training_vcf = VCF(file)
    for variant in training_vcf:
        trid = variant.INFO.get("TRID")
        if trid not in chrom:
            chrom[trid] = variant.CHROM

        (edit_ratio, sd) = get_variant_edr_sd(variant)

        edit_ratios.setdefault(trid, []).append(edit_ratio)
        sequencing_depths.setdefault(trid, []).append(sd)

    return True


def parse_sds_test(
    file: Path, sequencing_depths: dict = None, edit_ratios: dict = None, chrom: dict = None
) -> int | None:
    """Parse SDs from VCF. Return None if file was not found. Return count of individuals if found."""
    if not os.path.isfile(file):
        return None

    if chrom is None:
        chrom = {}

    training_vcf = VCF(file)
    individuals = training_vcf.samples
    nr_inds = len(individuals)

    for variant in training_vcf:
        trid = variant.INFO.get("TRID")
        if trid not in chrom:
            chrom[trid] = variant.CHROM

        edit_ratio = numpy.zeros(nr_inds)
        sd = numpy.zeros(nr_inds)

        for pos in range(nr_inds):
            (edit_ratio[pos], sd[pos]) = get_variant_edr_sd(variant, pos)

        edit_ratios[trid] = edit_ratio
        sequencing_depths[trid] = sd

    return nr_inds


def write_training_data(training_set: Path, data: List[dict]):
    """Write training set dictionaries to json file"""
    with open(training_set, "w") as f:
        json.dump(data, f)


def read_training_data(training_set: Path):
    """Read training set dictionaries to json file"""
    with open(training_set) as f:
        data = json.load(f)
    return data


def parse_training_data(training_set) -> Tuple[dict, dict]:
    """Open or parse training data files.
    Accepts a dir with data files, utilising only one sample per file. Returns training data dicts, keyed by locus id,
    with a list of values per locus.

    Alternatively, accepts a json file with the same kind of data.
    """
    if not os.path.isfile(training_set):
        training_data = {}
        training_edit_ratio = {}
        n_training_cases = 0
        for training_file in os.listdir(training_set):
            tf = os.path.join(training_set, training_file)
            if parse_sds_training(tf, training_data, training_edit_ratio):
                n_training_cases = n_training_cases + 1
        for trid in training_data.keys():
            training_data[trid] = sorted(training_data[trid])
            training_edit_ratio[trid] = sorted(training_edit_ratio[trid])
    else:
        training_data, training_edit_ratio = read_training_data(training_set)

    return (training_data, training_edit_ratio)


def get_total_set_p_edr_for_case(
    training_data: dict, nr_inds: int, test_data: dict, test_edit_ratio: dict, annotation: dict
) -> numpy.array:
    """Set P and EDR. Returns sample totals for counts, in a (N_samples, 1) numpy array."""

    case_total = numpy.zeros((nr_inds, 1))

    for trid in test_data:
        if trid not in training_data:
            continue

        annotation[trid] = {"p": numpy.zeros(nr_inds), "edit_ratio": numpy.zeros(nr_inds)}
        td = pd.Series(sorted(training_data[trid]))
        for pos in range(nr_inds):
            count_value = (td[td < test_data[trid][pos]]).sum()
            total_value = td.sum()
            case_total[pos] += test_data[trid][pos]
            p = count_value / total_value if total_value > 0 else 0
            annotation[trid]["p"][pos] = p
            annotation[trid]["edit_ratio"][pos] = test_edit_ratio[trid][pos]

    return case_total


def call_test_file(input_file: Path, xy: List, training_data: dict, alpha, edit, fraction) -> dict:
    """Parse test (case of interest) VCF. This is allowed to be a multisample VCF.
    Return annotation dict, containing per locus information. Each per locus value is a numpy array with the
    dimension (number_of_samples, 1).
    """
    annotation = {}
    test_data = {}
    test_edit_ratio = {}
    test_chrom = {}
    nr_inds = parse_sds_test(input_file, test_data, test_edit_ratio, test_chrom)

    p_threshold = alpha / len(test_data.keys())

    case_total = get_total_set_p_edr_for_case(
        training_data, nr_inds, test_data, test_edit_ratio, annotation
    )

    samples = VCF(input_file).samples
    case_total_n_trids = len(test_data.keys())
    case_average_depth = case_total / case_total_n_trids
    logger.info(f"Case average depth {case_average_depth}")

    for trid in test_data:
        if trid not in training_data:
            logger.warning(f"Skipping {trid}: not present in training data")
            continue

        # normalise depth ratio with sample average depth
        annotation[trid]["depth_ratio"] = numpy.zeros(nr_inds)
        annotation[trid]["coverage_warning"] = numpy.zeros(nr_inds, dtype=bool)
        annotation[trid]["coverage_drop"] = numpy.zeros(nr_inds, dtype=bool)

        for pos, sample in enumerate(samples):
            locus_depth = test_data[trid][pos] / case_average_depth[pos]
            annotation[trid]["depth_ratio"][pos] = locus_depth

            if (annotation[trid]["p"][pos] < p_threshold) and (test_edit_ratio[trid][pos] > edit):
                logger.info(
                    f"{trid} locus overall low with {test_data[trid][pos]} (P={annotation[trid]['p'][pos]}) and ratio is less over edit distance cutoff {test_edit_ratio[trid][pos]} (ind {pos})."
                )
                annotation[trid]["coverage_warning"][pos] = True

            fraction_cutoff = fraction
            sample_is_xy = sample in xy if xy is not None else False

            if sample_is_xy and ("X" in test_chrom[trid] or "Y" in test_chrom[trid]):
                fraction_cutoff = fraction - 0.5 if fraction > 0.5 else 0.05

            if (
                locus_depth < fraction_cutoff
                and trid in test_edit_ratio
                and test_edit_ratio[trid][pos] > edit
            ):
                logger.info(
                    f"{trid} locus coverage low with {test_data[trid][pos]}, below {fraction} of case average and edit distance ratio is over cutoff {test_edit_ratio[trid][pos]} (ind {pos})."
                )
                annotation[trid]["coverage_warning"][pos] = True

            if (
                locus_depth < fraction_cutoff
                and (annotation[trid]["p"][pos] < p_threshold)
                and trid in test_edit_ratio
                and test_edit_ratio[trid][pos] > edit
            ):
                logger.warning(f"Calling coverage drop for {trid}, ind {pos}.")
                annotation[trid]["coverage_drop"][pos] = True

    return annotation
