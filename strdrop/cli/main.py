import coloredlogs
import logging
import typer

import os

import Levenshtein
import pandas as pd

from cyvcf2 import VCF

from pathlib import Path
from typing import Optional
from typing_extensions import Annotated

from strdrop import __version__

coloredlogs.install(level="DEBUG")
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

ascii_logo = r"""
 ░▒▓███████▓▒░▒▓████████▓▒░▒▓███████▓▒░░▒▓███████▓▒░░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓███████▓▒░  
░▒▓█▓▒░         ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
░▒▓█▓▒░         ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░ 
 ░▒▓██████▓▒░   ░▒▓█▓▒░   ░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░  
       ░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
       ░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░        
░▒▓███████▓▒░   ░▒▓█▓▒░   ░▒▓█▓▒░░▒▓█▓▒░▒▓███████▓▒░░▒▓█▓▒░░▒▓█▓▒░░▒▓██████▓▒░░▒▓█▓▒░        
"""

app = typer.Typer(
    rich_markup_mode="rich",
    invoke_without_command=True,
    pretty_exceptions_show_locals=False,
    add_completion=False,
    help="Call coverage drops over alleles in STR VCFs"
)

EDIT_DISTANCE_CUTOFF = 0.9

def parse_sds(file: Path, training_data:dict = {}, edit_ratios:dict={}) -> bool:
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


@app.callback()
def main(version: Annotated[Optional[bool], typer.Option("--version", "-v", is_flag=True, help="Show version and exit.")]=False):
    if version:
        typer.echo(f"Strdrop version {__version__}")
        raise typer.Exit()

@app.command(
    help="[bold]STRdrop[/bold]: Detect drops in STR coverage 🧬",
    no_args_is_help=True,
)
def call(training_set: Annotated[Path,
                                    typer.Option(exists=True, dir_okay=True, help="Input directory with reference data")],
            input_file: Annotated[Path, typer.Argument(help="Input STR call VCF file")],
    ) -> None:
    """
    Call drops, given reference files dir

    (Dev: Maybe add a manifest, possibly including karyotype or just dir)
    """

    logger.info(ascii_logo)

    training_data = {}
    training_edit_ratio = {}
    n_training_cases = 0
    for training_file in os.listdir(training_set):
        tf = os.path.join(training_set, training_file)
        if parse_sds(tf, training_data, training_edit_ratio):
            n_training_cases = n_training_cases + 1

    test_data = {}
    test_edit_ratio = {}
    parse_sds(input_file, test_data, test_edit_ratio)

    alpha = 0.05
    p_threshold = alpha / len(test_data.keys())

    # plain per allele counts
    case_total=0
    for trid in test_data.keys():
        td = pd.Series(sorted(training_data[trid]))
        count_value = (td[td<test_data[trid][0]]).sum()
        total_value = td.sum()
        case_total += test_data[trid][0]
        p = count_value / total_value if total_value > 0 else 0
        if(p < p_threshold) and trid in test_edit_ratio and test_edit_ratio[trid][0] > EDIT_DISTANCE_CUTOFF:
            result = f"{trid} locus overall low with {test_data[trid][0]} (P={p}) and ratio is less than cutoff {test_edit_ratio[trid]}."
            print(result)

    case_total_n_trids = len(test_data.keys())
    case_average_depth = case_total / case_total_n_trids
    print(f"case average depth {case_average_depth}")
    for trid in test_data.keys():
        locus_depth = test_data[trid][0] / case_average_depth
        if locus_depth < 0.5 and trid in test_edit_ratio and test_edit_ratio[trid][0] > EDIT_DISTANCE_CUTOFF:
            result = f"{trid} locus coverage low with {test_data[trid][0]}, below 0.5 of case average and ratio is less than cutoff {test_edit_ratio[trid]}."
            print(result)

def run():
    app()