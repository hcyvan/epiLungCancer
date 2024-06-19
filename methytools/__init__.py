import argparse
from .tools import analysis_genomic_element


def main():
    parser = argparse.ArgumentParser(description='meth matrix filter and select')
    parser.add_argument('-s', '--input', required=True, help='The input bed ratio file with header')
    parser.add_argument('-e', '--element', choices=['tss', 'cgi'], default='tss', help='tss or cgi')
    parser.add_argument('-u', '--up', default=10000, type=int, help='tss upstream')
    parser.add_argument('-d', '--down', default=10000, type=int, help='tss downstream')
    parser.add_argument('-o', '--out', help='The output signal file')
    parser.add_argument('-v', '--version', action='version', version='0.1')
    args = parser.parse_args()
    analysis_genomic_element(args.input, args.element, args.out, args.up, args.down)
