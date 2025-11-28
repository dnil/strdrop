import json
import os
import Levenshtein
import pandas as pd

from cyvcf2 import VCF, Variant
from pathlib import Path
from typing import List, Tuple

import logging
logger = logging.getLogger(__name__)


def get_allele(variant:Variant, pos:int, step:int) -> str | None:
    """Return allele for variant sample number pos and allele step (0, 1)"""
    if variant.genotypes[pos][step] == 0:
        return variant.REF
    
    idx = variant.genotypes[pos][step] - 1
    return variant.ALT[idx]


def get_variant_edr_sd(variant: Variant, ind_nr: int) -> Tuple[float, float]:
    """Return edit ratio and sequencing depth for variant, individual"""
    a1 = get_allele(variant, ind_nr, 0)
    a2 = get_allele(variant, ind_nr, 1)

    edit_ratio = Levenshtein.ratio(a1, a2)

    ref_sd = 0
    alt_sd = 0
    sd_values = variant.format('SD')[0]
    if len(sd_values) == 1:
        alt_sd = max(int(sd_values[0]), 0)
    if len(sd_values) == 2:
        alt_sd = max(int(sd_values[1]), 0)
        ref_sd = max(int(sd_values[0]), 0)

    sd = ref_sd + alt_sd

    return (edit_ratio, sd)

def parse_sds(file: Path, training_data:dict = {}, edit_ratios:dict={}, chrom:dict={}, ind_nr: int = 0) -> bool:
    """Parse SDs from VCF. Return False if file was not found."""
    if not os.path.isfile(file):
        return False
        
    training_vcf = VCF(file)
    for variant in training_vcf:
        (edit_ratio, sd) = get_variant_edr_sd(variant, ind_nr)

        trid = variant.INFO.get('TRID')
        if trid not in chrom:
            chrom[trid] = variant.CHROM

        append_safe(edit_ratios, trid, edit_ratio)
        append_safe(training_data, trid, sd)

    return True


def append_safe(obj, obj_index, elem):
    """Append `elem` to list in `obj` at `obj_index`.
    If no list exists `elem` will be first element catching
    the KeyError raised."""
    try:
        obj[obj_index].append(elem)
    except KeyError:
        obj[obj_index] = [elem]


def write_training_data(training_set: Path, data: List[dict]):
    """Write training set dictionaries to json file """
    with open(training_set, "w") as f:
        json.dump(data, f)


def read_training_data(training_set: Path):
    """Read training set dictionaries to json file"""
    with open(training_set) as f:
        data = json.load(f)
    return data


def parse_training_data(training_set)-> Tuple[dict, dict]:

    if not os.path.isfile(training_set):
        training_data = {}
        training_edit_ratio = {}
        n_training_cases = 0
        for training_file in os.listdir(training_set):
            tf = os.path.join(training_set, training_file)
            if parse_sds(tf, training_data, training_edit_ratio):
                n_training_cases = n_training_cases + 1
        for trid in training_data.keys():
            training_data[trid] = sorted(training_data[trid])
            training_edit_ratio[trid] = sorted(training_edit_ratio[trid])
    else:
        training_data, training_edit_ratio = read_training_data(training_set)

    return (training_data, training_edit_ratio)


def call_test_file(input_file: Path, xy: bool, training_data:dict, alpha, edit, fraction) -> dict:

    annotation = {}
    test_data = {}
    test_edit_ratio = {}
    test_chrom = {}
    parse_sds(input_file, test_data, test_edit_ratio, test_chrom)

    p_threshold = alpha / len(test_data.keys())

    case_total = 0
    for trid in test_data.keys():
        td = pd.Series(sorted(training_data[trid]))
        count_value = (td[td < test_data[trid][0]]).sum()
        total_value = td.sum()
        case_total += test_data[trid][0]
        p = count_value / total_value if total_value > 0 else 0
        annotation[trid] = {"p": p, "edit_ratio": test_edit_ratio[trid][0]}

    case_total_n_trids = len(test_data.keys())
    case_average_depth = case_total / case_total_n_trids
    logger.info(f"Case average depth {case_average_depth}")

    for trid in test_data.keys():
        locus_depth = test_data[trid][0] / case_average_depth
        annotation[trid]["depth_ratio"] = locus_depth

        if (annotation[trid]["p"] < p_threshold) and (test_edit_ratio[trid][0] > edit):
            logger.info(f"{trid} locus overall low with {test_data[trid][0]} (P={p}) and ratio is less over edit distance cutoff {test_edit_ratio[trid]}.")
            annotation[trid]["coverage_warning"] = True

        fraction_cutoff = fraction
        if xy and "X" in test_chrom[trid] or "Y" in test_chrom[trid]:
            fraction_cutoff = fraction - 0.5 if fraction - 0.5 > 0 else 0.05

        if locus_depth < fraction_cutoff and trid in test_edit_ratio and test_edit_ratio[trid][0] > edit:
            logger.info(f"{trid} locus coverage low with {test_data[trid][0]}, below {fraction} of case average and edit distance ratio is over cutoff {test_edit_ratio[trid]}.")
            annotation[trid]["coverage_warning"] = True

        if locus_depth < fraction_cutoff and (annotation[trid]["p"] < p_threshold) and trid in test_edit_ratio and test_edit_ratio[trid][0] > edit:
            logger.warning(f"Calling coverage drop for {trid}")
            annotation[trid]["coverage_drop"] = True

    return annotation