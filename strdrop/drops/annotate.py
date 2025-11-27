import logging

from cyvcf2 import VCF, Writer

from pathlib import Path

logger = logging.getLogger(__name__)

def write_output(input_file: Path, annotation: dict, output_file: Path):
    """Write out VCF"""

    vcf = VCF(input_file)

    vcf.add_info_to_header({'ID': 'STRDROP_P', 'Description': 'Strdrop coverage sequencing depth level probability',
                            'Type': 'Float', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'STRDROP_EDR', 'Description': 'Strdrop allele similarity Levenshtein edit distance ratio',
                            'Type': 'Float', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'STRDROP_SDR', 'Description': 'Strdrop case average adjusted sequencing depth ratio',
                            'Type': 'Float', 'Number': '1'})

    vcf.add_info_to_header({'ID': 'STRDROP', 'Description': 'Strdrop coverage drop detected',
                            'Type': 'Flag', 'Number': '0'})

    vcf.add_filter_to_header({'ID': 'LowDepth', 'Description': 'Strdrop coverage drop detected'})
    # Add LowDepth flag?

    w = Writer(output_file, vcf)

    for v in vcf:
        trid = v.INFO.get('TRID')
        if trid in annotation:
            v.INFO['STRDROP_P'] = annotation[trid]["p"]
            v.INFO['STRDROP_EDR'] = annotation[trid]["edit_ratio"]
            v.INFO['STRDROP_SDR'] = annotation[trid]["depth_ratio"]
            if "coverage_drop" in annotation[trid]:
                v.INFO['STRDROP'] = annotation[trid]["coverage_drop"]
                filter = v.FILTER
                if not filter:
                    filter = "LowDepth"
                else:
                    filter += ";LowDepth"
                v.FILTER = filter
        w.write_record(v)

    w.close()
    vcf.close()