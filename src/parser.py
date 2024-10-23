from argparse import ArgumentParser
from Bio import SeqIO


def args():
    parser = ArgumentParser()
    parser.add_argument("query", type=str)
    parser.add_argument("db", type=str)
    parser.add_argument("-v", "--verbose", help="verbose",
                        type=bool, default=False)
    parser.add_argument("-m", "--mode", help="mode of the algorihtm", type=str)
    return parser.parse_args()


def open_fasta(path: str):
    return SeqIO.parse(path, "fasta")
