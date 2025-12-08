import coloredlogs
import logging
import typer

from pathlib import Path
from typing import List, Optional
from typing_extensions import Annotated

from strdrop import __version__
from strdrop.drops import call_test_file, parse_training_data, write_output, write_training_data

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
    help="Call coverage drops over alleles in STR VCFs",
)

EDIT_RATIO_CUTOFF = 0.9
CASE_COVERAGE_RATIO_CUTOFF = 0.55
ALPHA = 0.05


@app.callback()
def main(
    version: Annotated[
        Optional[bool], typer.Option("--version", "-v", is_flag=True, help="Show version and exit.")
    ] = False,
):
    if version:
        typer.echo(f"Strdrop version {__version__}")
        raise typer.Exit()


@app.command(
    help="[bold]STRdrop[/bold]: Detect drops in STR coverage 🧬",
    no_args_is_help=True,
)
def call(
    training_set: Annotated[
        Path,
        typer.Option(
            exists=True, dir_okay=True, help="Training VCF directory or json with reference data"
        ),
    ],
    input_file: Annotated[Path, typer.Argument(help="Input STR call VCF file")],
    output_file: Annotated[Path, typer.Argument(help="Output annotated VCF file")],
    xy: Annotated[
        Optional[List[str]],
        typer.Option(
            help="Sample in VCF to treat as karyotype XY. May be specified multiple times.",
            metavar="SAMPLE ID",
        ),
    ] = None,
    alpha: Annotated[
        float, typer.Option(help="Unadjusted probability confidence level for coverage test")
    ] = ALPHA,
    fraction: Annotated[
        float, typer.Option(help="Case average adjusted sequencing depth ratio cutoff")
    ] = CASE_COVERAGE_RATIO_CUTOFF,
    edit: Annotated[
        float, typer.Option(help="Allele similarity Levenshtein edit distance ratio cutoff")
    ] = EDIT_RATIO_CUTOFF,
) -> None:
    """
    Call drops, given reference files dir

    (Dev: Maybe add a manifest, possibly including karyotype)
    """

    logger.info(ascii_logo)

    (training_data, training_edit_ratio) = parse_training_data(training_set)

    annotation = call_test_file(input_file, xy, training_data, alpha, edit, fraction)

    write_output(input_file, annotation, output_file)


@app.command(
    help="[bold]STRdrop[/bold]: Build reference json from sequencing coverage in STR VCFs",
    no_args_is_help=True,
)
def build(
    training_set: Annotated[
        Path, typer.Option(exists=True, dir_okay=True, help="Input directory with reference data")
    ],
    reference_file: Annotated[Path, typer.Argument(help="Output reference archive")],
) -> None:
    (training_data, training_edit_ratio) = parse_training_data(training_set)
    write_training_data(reference_file, [training_data, training_edit_ratio])


def run():
    app()
