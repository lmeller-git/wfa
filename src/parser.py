from argparse import ArgumentParser
from Bio import SeqIO


def args():
    parser = ArgumentParser()
    parser.add_argument("query", type=str)
    parser.add_argument("db", type=str)
    parser.add_argument("-v", "--verbose", help="verbose",
                        default=False, action="store_true")
    parser.add_argument("-m", "--mode", help="mode of the algorihtm",
                        type=str, choices=["global", "local", "semi_global"], default="global")
    parser.add_argument("-p", "--plot", default=False,
                        action="store_true", help="create an animation of the algorithm")
    parser.add_argument("-o", "--out", help="path to save results",
                        default="./temp/animation", type=str)
    parser.add_argument("-c", "--compare", default=False,
                        action="store_true", help="compare needleman wunsch with wfa")
    return parser.parse_args()


def open_fasta(path: str):
    return SeqIO.parse(path, "fasta")
