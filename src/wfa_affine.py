from Bio import SeqIO, SeqRecord, Seq
from argparse import Namespace
from dataclasses import dataclass
from typing import Self
from src.utils import timeit
from enum import StrEnum
from collections import defaultdict

"""
diagonal and offsets impkementation for gap affine scoring 
"""

type Records = SeqIO.FastIO.FastaIterator
type Record = SeqRecord.SeqRecord
type Sequence = Seq.Seq


@timeit
def align(query: Records, db: Records, args: Namespace) -> None:
    global scheme
    scheme = ScoringScheme(0, 2, 2, 6)
    for q in query:
        for d in db:
            ocean = Ocean(q, d, args)
            if ocean.align():
                print("i guess")


@dataclass
class Point:
    x: int
    y: int


@dataclass
class ScoringScheme:
    match_: int
    mismatch: int
    gap_opening: int
    gap_extension: int


@dataclass
class WaveFront:
    M: int
    I: int
    D: int
    s: int
    traceback: dict


def extend_wf(wf: WaveFront, q: Sequence, db: Sequence, k):
    while q[wf.M] == db[wf.M - k - 1]:
        wf.M += 1


def new_wf(*args):
    return WaveFront(0, 0, 0, 0)


def new_d(*args):
    return defaultdict(new_wf)


class Ocean:
    def __init__(self, query: Record, db: Record, args: Namespace) -> None:
        self.query = query
        self.db = db
        self.mode = args.mode
        self.verbose = args.verbose
        # starting condition
        self.wavefronts = defaultdict(new_d)

    def backtrace(self) -> None:
        pass

    def align(self) -> bool:

        X = 0
        diagonal_offset = len(self.query) - len(self.db)
        max_offset = len(self.query)

        current_score = 0
        while True:
            for k, wf in self.wavefronts[current_score].items():
                extend_wf(wf, self.query, self.db, k)
                if wf.M >= max_offset and k == diagonal_offset:
                    # Pass the final wavefront for backtracking
                    self.backtrace(wf)
                    return True

            current_score += min(scheme.mismatch, scheme.gap_extension)
            self.get_next(current_score)
        return False
        """
        # starting condition
        X = 0
        diagonal_offset = len(self.query) - len(self.db)
        max_offset = len(self.query)
        current_score = 0
        while True:
            for k, wf in self.wavefronts[current_score].items():
                extend_wf(wf, self.query, self.db, k)
                if wf.M >= max_offset and k == diagonal_offset:
                    self.backtrace()
                    return True
            current_score += min(scheme.mismatch, scheme.gap_extension)
            self.get_next(current_score)
        return False
        """

    def get_next(self, s):

        new_wfs = defaultdict(new_wf)
        M_prev = self.wavefronts[s - scheme.mismatch]
        I_prev = self.wavefronts[s - scheme.gap_extension]
        D_prev = self.wavefronts[s - scheme.gap_extension]

        low = min(min(M_prev.keys()), min(I_prev.keys()), min(D_prev.keys()))
        high = max(max(M_prev.keys()), max(I_prev.keys()), max(D_prev.keys()))

        for k in range(low, high + 1):
            new_wfs[k] = WaveFront(
                M=max(M_prev[k].M + scheme.match_, I_prev[k].I, D_prev[k].D),
                I=max(M_prev[k - 1].I + scheme.gap_extension, I_prev[k - 1].I),
                D=max(M_prev[k + 1].D + scheme.gap_extension, D_prev[k + 1].D),
                traceback={'M': M_prev[k], 'I': I_prev[k], 'D': D_prev[k]}
            )

        self.wavefronts[s] = new_wfs
        """
           new_wfs = defaultdict(new_wf)
           M1 = self.wavefronts[s - scheme.mismatch]
           M2 = self.wavefronts[s - scheme.gap_opening - scheme.gap_extension]
           I_D = self.wavefronts[s - scheme.gap_extension]
           low = min(M1[min(M1.keys())], M2[min(M2.keys())], I_D[min(I_D.keys())])
           high = max(M1[max(M1.keys())], M2[max(M2.keys())],
                      I_D[max(I_D.keys())])
           for k in range(low, high):
               new_wfs[k] = WaveFront(k, max(
                   M1[k] + 1, I_D[k]), max(M2[k - 1], I_D[k - 1]) + 1, max(M2[k + 1], I_D[k + 1]))
           self.wavefronts[s] = new_wfs
           """

    def prune(self) -> None:
        pass
