from enum import StrEnum
from src.utils import timeit
from typing import Self
from dataclasses import dataclass
from argparse import Namespace
from Bio import SeqIO, SeqRecord, Seq


type Records = SeqIO.FastIO.FastaIterator
type Record = SeqRecord.SeqRecord
type Sequence = Seq.Seq


@timeit
def align(query: Records, db: Records, args: Namespace) -> None:
    global scheme
    scheme = ScoringScheme(0, 2, 0, 4)
    for q in query:
        for d in db:
            # print(q.seq)
            # print(d.seq)
            ocean = Ocean(q, d, args)
            while not ocean.step():
                ocean.prune()


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


class WaveFront:
    def __init__(self, diag: int, offset: int, q: Sequence, db: Sequence) -> None:
        self.offset = offset
        self.diag = diag
        self.q = q
        self.db = db

    def expand(self) -> bool:
        if self.is_converged():
            return True
        if self.is_at_end():
            return False
        # print(self.offset, self.diag, "huhu")
        while self.q[self.offset + max(self.diag, 0)] == self.db[self.offset - self.diag]:
            self.offset += 1
            # print(self.offset, self.diag)
            if self.is_converged():
                return True
            if self.is_at_end():
                return False
        return False

    def is_converged(self) -> bool:
        return self.diag == len(self.q) - len(self.db) and self.offset >= len(self.q)

    def is_at_end(self):
        return self.offset + max(0, self.diag) >= len(self.q) or self.offset - self.diag >= len(self.db) or self.offset - self.diag < 0


def fill_diff(wave1: WaveFront, wave2: WaveFront, q: Sequence, db: Sequence, q_str: str, db_str: str):
    delta_offset = wave1.offset - wave2.offset
    # print("vars:")
    # print(wave1.diag, wave1.offset)
    # print(wave2.diag, wave2.offset)
    # print(delta_offset)
    # print("prior:")
    # print(q_str)
    # print(db_str)
    # print("after:")
    q_str_ = q[wave1.offset -
               delta_offset:wave1.offset] + q_str
    db_str_ = db[wave1.offset - wave1.diag -
                 delta_offset:wave1.offset - wave1.diag] + db_str

    # print(q_str_)
    # print(db_str_)
    # print("\n")
    return q_str_, db_str_


class Ocean:
    def __init__(self, query: Record, db: Record, args: Namespace) -> None:
        self.query = query
        self.db = db
        self.mode = args.mode
        self.verbose = args.verbose
        self.wavefronts = {0: [WaveFront(0, 0, self.query.seq, self.db.seq)]}
        self.current_score = 0
        self.max_wf_len = 100

    def backtrace(self, wave: WaveFront, q: str, db: str, score: int) -> None:
        if score == 0:
            q_, db_ = fill_diff(wave, WaveFront(
                0, 0, self.query.seq, self.db.seq), self.query.seq, self.db.seq, q, db)
            print(f"found alignemnt with score {self.current_score}")
            print(f"query: {q_}")
            print("       ", end='')
            [print(' ' if q__ != d__ else '|', end='')
             for q__, d__ in zip(q_, db_)]
            print(f"\ndb:    {db_}")
            print("\n")
            return
        else:
            # print(f"Score: {score}, Offset: {
            #     wave.offset}, Diagonal: {wave.diag}")
            # print(f"Query alignment so far: {q}")
            # print(f"DB alignment so far: {db}")
            gscore = score - scheme.gap_extension
            mscore = score - scheme.mismatch
            if gscore in self.wavefronts:
                for wf in self.wavefronts[gscore]:
                    # since I keep only the furthest reaching point per diagonal, i know that this is the true predecessor
                    if wf.offset > wave.offset:
                        continue
                    if wf.diag == wave.diag - 1:
                        q_, db_ = fill_diff(
                            wave, wf, self.query.seq, self.db.seq, q, db)
                        self.backtrace(
                            wf, self.query.seq[wf.offset - 1] + q_, '-' + db_, gscore)
                        # return
                    if wf.diag == wave.diag + 1:
                        q_, db_ = fill_diff(
                            wave, wf, self.query.seq, self.db.seq, q, db)
                        self.backtrace(
                            wf, '-' + q_, self.db.seq[wf.offset - wf.diag] + db_, gscore)
                        # return

            if mscore in self.wavefronts:
                for wf in self.wavefronts[mscore]:
                    if wf.diag == wave.diag:
                        q_, db_ = fill_diff(
                            wave, wf, self.query.seq, self.db.seq, q, db)
                        self.backtrace(wf, q_, db_, mscore)
                        # return

    def step(self) -> bool:
        # print(self.current_score, len(
        # self.wavefronts[self.current_score]), "wtf")
        for wf in self.wavefronts[self.current_score]:
            if wf.expand():
                self.backtrace(wf, "", "", self.current_score)
                return True
        self.current_score += min(scheme.mismatch, scheme.gap_extension)
        self.next_wave()
        return False

    def next_wave(self) -> None:
        gscore = self.current_score - scheme.gap_extension
        mscore = self.current_score - scheme.mismatch
        new_wfs = []
        if gscore in self.wavefronts:
            for wave in self.wavefronts[gscore]:
                new_wfs += [WaveFront(wave.diag + 1, wave.offset, wave.q, wave.db),
                            WaveFront(wave.diag - 1, wave.offset, wave.q, wave.db)]

        if mscore in self.wavefronts:
            for wave in self.wavefronts[mscore]:
                new_wfs += [WaveFront(wave.diag,
                                      wave.offset + 1, wave.q, wave.db)]
        self.wavefronts[self.current_score] = new_wfs

    def is_prunable(self, wave: WaveFront, max_diag: int, max_offset: int, min_diag: int,  min_d: int, idx):
        if wave.diag < min_diag or wave.diag > max_diag or wave.offset > max_offset:
            return True
        for i, wv in enumerate(self.wavefronts[self.current_score]):
            if i == idx:
                continue
            if (wv.diag == wave.diag and wv.offset == wave.offset) or (wv.diag == wave.diag and wv.offset > wv.offset):
                return True
        return False
        if wave.diag > 0 and len(self.query.seq) > len(self.db.seq) or wave.diag < 0 and len(self.query.seq) < len(self.db.seq):
            return True
        return False

    def get_threshold(self):
        d_min = 0
        for wave in self.wavefronts[self.current_score]:
            d_q = len(self.query.seq) - wave.offset + wave.diag
            d_db = len(self.db.seq) - wave.offset
            d = max(d_q, d_db)
            d_min = min(d, d_min)
        return 2 * d_min

    def prune(self) -> None:
        max_diag = len(self.query.seq)
        min_diag = -len(self.db.seq)
        max_offset = len(self.query.seq)
        # diag_offset = len(self.query.seq) - len(self.db.seq)
        # if len(self.wavefronts[self.current_score]) > self.max_wf_len:
        # threshold = self.get_threshold()
        # else:
        threshold = len(self.query.seq)
        self.wavefronts[self.current_score] = [
            wave for i, wave in enumerate(self.wavefronts[self.current_score])
            if not self.is_prunable(wave, max_diag, max_offset, min_diag, threshold, i)
        ]
