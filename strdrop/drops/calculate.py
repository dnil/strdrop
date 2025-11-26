import os
import Levenshtein
import pandas as pd

from pathlib import Path
from typing import Tuple

import logging
logger = logging.getLogger(__name__)

from cyvcf2 import VCF

def parse_sds(file: Path, training_data:dict = {}, edit_ratios:dict={}, chrom:dict={}) -> bool:
    """Parse SDs from VCF. Return False if file was not found."""
    if os.path.isfile(file):
        training_vcf = VCF(file)
        for variant in training_vcf:

            trid = variant.INFO.get('TRID')

            a1 = variant.REF
            a2 = variant.REF
            if len(variant.ALT) == 1:
                a2 = variant.ALT[0]
            elif len(variant.ALT) > 1:
                a2 = variant.ALT[1]
            edit_ratio = Levenshtein.ratio(a1, a2)
            if trid not in chrom:
                chrom[trid] = variant.CHROM

            if trid in edit_ratios:
                edit_ratios[trid].append(edit_ratio)
            else:
                edit_ratios[trid] = [edit_ratio]

            ref_sd = 0
            alt_sd = 0
            sd_values = variant.format('SD')[0]
            if len(sd_values) == 1:
                alt_sd = int(sd_values[0])
            if len(sd_values) == 2:
                alt_sd = int(sd_values[1])
                ref_sd = int(sd_values[0])

            if trid in training_data:
                # training_data[trid].append(ref_value)
                training_data[trid].append(ref_sd + alt_sd)
                continue
            training_data[trid] = [ref_sd + alt_sd]
        return True
    return False

def parse_training_data(training_set)-> Tuple[dict, dict]:
    training_data = {}
    training_edit_ratio = {}
    n_training_cases = 0
    for training_file in os.listdir(training_set):
        tf = os.path.join(training_set, training_file)
        if parse_sds(tf, training_data, training_edit_ratio):
            n_training_cases = n_training_cases + 1

    return (training_data, training_edit_ratio)

def call_test_file(input_file: Path, xy: bool, training_data:dict, alpha, edit, fraction) -> dict:

    annotation = {}
    test_data = {}
    test_edit_ratio = {}
    parse_sds(input_file, test_data, test_edit_ratio)

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

        if (annotation[trid]["p"] < p_threshold) and test_edit_ratio[trid][0] > edit:
            logger.info(f"{trid} locus overall low with {test_data[trid][0]} (P={p}) and ratio is less over edit distance cutoff {test_edit_ratio[trid]}.")
            annotation[trid]["coverage_warning"] = True

        if locus_depth < fraction and trid in test_edit_ratio and test_edit_ratio[trid][0] > edit:
            logger.info(f"{trid} locus coverage low with {test_data[trid][0]}, below 0.5 of case average and edit distance ratio is over cutoff {test_edit_ratio[trid]}.")
            annotation[trid]["coverage_warning"] = True

        if locus_depth < fraction and (annotation[trid]["p"] < p_threshold) and trid in test_edit_ratio and test_edit_ratio[trid][0] > edit:
            logger.warning(f"Calling coverage drop for {trid}")
            annotation[trid]["coverage_drop"] = True

    return annotation