import argparse
from .tools import analysis_genomic_element, analysis_genomic_element_from_gz
from .utils import is_gz_file


def main():
    parser = argparse.ArgumentParser(description='Tools for processing methylation data')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')
    # =====================================================================
    parser_signal = subparsers.add_parser('signal',
                                          help='Extract the signal value of a region in the genome')
    parser_signal.add_argument('-i', '--input', required=True, help='The input bed ratio file with header')
    parser_signal.add_argument('-e', '--element', choices=['tss', 'cgi'], default='tss', help='tss or cgi')
    parser_signal.add_argument('-u', '--up', default=10000, type=int, help='tss upstream')
    parser_signal.add_argument('-d', '--down', default=10000, type=int, help='tss downstream')
    parser_signal.add_argument('-o', '--out', help='The output signal file')
    parser_signal.add_argument('-v', '--version', action='version', version='0.1')
    # =====================================================================
    parser_region_methy = subparsers.add_parser('region-methy',
                                                help='Extract the signal value of a region in the genome')
    parser_region_methy.add_argument('-i', '--input', required=True, help='The input bed ratio file with header')
    parser_region_methy.add_argument('-e', '--element', choices=['tss', 'cgi'], default='tss', help='tss or cgi')
    parser_region_methy.add_argument('-u', '--up', default=-50, type=int,
                                     help='the left bounder to the center, center is 0,left is negative, right is positive')
    parser_region_methy.add_argument('-d', '--down', default=50, type=int,
                                     help='the right bounder to the center, center is 0,left is negative, right is positive')
    parser_region_methy.add_argument('-o', '--out', help='The output signal file')
    parser_region_methy.add_argument('-v', '--version', action='version', version='0.1')

    args = parser.parse_args()
    if args.sub == 'signal':
        analysis_genomic_element(args.input, args.element, args.out, args.up, args.down)
    elif args.sub == 'region-methy':
        if is_gz_file(args.input):
            analysis_genomic_element_from_gz(args.input, args.element, args.out, args.up, args.down)
        else:
            raise Exception("You should input a bgziped and indexed bed file")
    else:
        parser.print_help()
