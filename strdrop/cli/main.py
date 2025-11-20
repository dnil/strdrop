import coloredlogs
import logging
import typer

import os

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

def parse_sds(file: Path, training_data:dict = {}):

    if os.path.isfile(file):
        training_vcf = VCF(file)
        for variant in training_vcf:
            trid = variant.INFO.get('TRID')
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
            training_data[trid] = [alt_sd]



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

    for training_file in os.listdir(training_set):
        tf = os.path.join(training_set, training_file)
        parse_sds(tf, training_data)


    test_data = {}
    parse_sds(input_file, test_data)

    for trid in test_data.keys():
        td = pd.Series(sorted(training_data[trid]))
        count_value = (td[td<test_data[trid][0]]).sum()
        total_value = td.sum()
        p = count_value / total_value if total_value > 0 else 0
        if(p < 0.05):
            print(f"{trid} low with {test_data[trid][0]} (P={p})")

def run():
    app()