from argparse import ArgumentError, ArgumentParser, FileType
from magstruct import plane, transformations
from numpy import float as npfloat

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--input', type=FileType('r'),
        help='file to read RA/Dec/Distance values from')
    parser.add_argument('-o', '--output', type=str,
        help='directory to output plots to')
    parser.add_argument('-c', '--center', type=npfloat, nargs=3,
        metavar='d RA Dec',
        help='coordinates to center of target galaxy')

    args = parser.parse_args()

    return args

def main():
    print("Hello")


if __name__ == '__main__':
    exit(main())
