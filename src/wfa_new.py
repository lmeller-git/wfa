from Bio import SeqIO, SeqRecord, Seq
from argparse import Namespace
from dataclasses import dataclass
from typing import Self
from src.utils import timeit
from enum import StrEnum

type Records = SeqIO.FastIO.FastaIterator
type Record = SeqRecord.SeqRecord
type Sequence = Seq.Seq


@timeit
def align(query: Records, db: Records, args: Namespace) -> None:
    global scheme
    scheme = ScoringScheme(0, 2, 0, 4)
    for q in query:
        for d in db:
            print(q.seq)
            print(d.seq)
            ocean = Ocean(q, d, args)
            while not ocean.step():
                ocean.prune()
                pass
            ocean.backtrace()


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


class Move(StrEnum):
    Insertion = "insertion"
    Deletion = "deletion"
    MatchMismatch = "match_mismatch"
    NoMove = "none"


@dataclass
class MoveInfo:
    move: Move
    dscore: int
    dx: int
    dy: int


class MoveIterator:

    def __init__(self):
        self.c = 0

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> MoveInfo:
        match self.c:
            case 0:
                self.c += 1
                return MoveInfo(Move.Deletion, scheme.gap_extension, 1, 0)
            case 1:
                self.c += 1
                return MoveInfo(Move.MatchMismatch, scheme.mismatch, 1, 1)
            case 2:
                self.c += 1
                return MoveInfo(Move.Insertion, scheme.gap_extension, 0, 1)
            case _:
                raise StopIteration


class WaveFront:
    def __init__(self, x: int, y: int, score: int, query: Sequence, db: Sequence, all_moves: list[Move]) -> None:
        self.x = x
        self.y = y
        self.score = score
        self.db = db
        self.query = query
        self.all_moves = all_moves

    def expand(self) -> bool:
        while (self.x < len(self.db) and self.y < len(self.query)) and (self.db[self.x] == self.query[self.y]):
            # print("expand")
            self.x += 1
            self.y += 1
            self.score += scheme.match_
            self.all_moves.append(Move.MatchMismatch)
            if self.is_converged():
                return True
        return False

    def is_converged(self) -> bool:
        return self.x == len(self.db) and self.y == len(self.query)

    def pprint(self):
        db = ""
        q = ""
        in_db = 1
        in_q = 1
        moves = iter(reversed(self.all_moves))
        while (move := next(moves)) != Move.NoMove:
            match move:
                case Move.Deletion:
                    db = self.db[len(self.db) - in_db] + db
                    q = '-' + q
                    in_db += 1
                    pass
                case Move.Insertion:
                    db = '-' + db
                    q = self.query[len(self.query) - in_q] + q
                    in_q += 1
                    pass
                case Move.MatchMismatch:
                    db = self.db[len(self.db) - in_db] + db
                    q = self.query[len(self.query) - in_q] + q
                    in_db += 1
                    in_q += 1
                    pass
                case _:
                    pass

        def item(x, y): return ' ' if x != y else '|'
        print(db)
        [print(item(q_, d_), end='') for q_, d_ in zip(q, db)]
        print()
        print(q)
        print(self.score)

    def get_fr(self) -> tuple[int, int]:
        diag = self.x - self.y
        if diag <= 0:
            offset = self.y
        else:
            offset = self.x
        return diag, offset


class Ocean:
    def __init__(self, query: Record, db: Record, args: Namespace) -> None:
        self.query = query
        self.db = db
        self.mode = args.mode
        self.verbose = args.verbose
        self.current_score = 0
        self.wavefronts = {0: [
            WaveFront(0, 0, 0, query.seq, db.seq, [Move.NoMove])]}
        # print("i:", *[i for i in self.wavefronts])

    def backtrace(self) -> None:
        pass

    def step(self) -> bool:
        if self.current_score in self.wavefronts.keys():
            # print(len(self.wavefronts[self.current_score]))
            for wf in self.wavefronts[self.current_score]:
                if wf.expand() or wf.is_converged():
                    wf.pprint()
                    return True
        self.new_wave()
        # print(self.current_score)
        return False

    def new_wave(self) -> None:
        wfs = []
        self.current_score += 1
        mscore = self.current_score - scheme.mismatch
        gscore = self.current_score - scheme.gap_extension
        if mscore in self.wavefronts:
            # print("y", len(self.wavefronts[mscore]))
            for wf in self.wavefronts[mscore]:
                wfs.append(WaveFront(wf.x + 1, wf.y + 1, self.current_score,
                           self.query.seq, self.db.seq, wf.all_moves + [Move.MatchMismatch]))
        if gscore in self.wavefronts:
            # print("yu", len(self.wavefronts[gscore]))
            for wf in self.wavefronts[gscore]:
                wfs += [WaveFront(wf.x + 1, wf.y, self.current_score, self.query.seq, self.db.seq, wf.all_moves + [Move.Deletion]),
                        WaveFront(wf.x, wf.y + 1, self.current_score, self.query.seq, self.db.seq, wf.all_moves + [Move.Insertion])]
        if len(wfs) > 0:
            self.wavefronts[self.current_score] = wfs

    def can_prune(self, wf: WaveFront, i: int) -> bool:
        if wf.x > len(self.db.seq) or wf.y > len(self.query.seq):
            return True
        return False
        d, o = wf.get_fr()
        for j, wave in enumerate(self.wavefronts[self.current_score]):
            d_, o_ = wave.get_fr()
            if d_ == d and o_ >= o:
                return True
        return False

    def prune(self) -> None:
        if self.current_score not in self.wavefronts:
            return
        self.wavefronts[self.current_score] = [
            wf for i, wf in enumerate(self.wavefronts[self.current_score]) if not self.can_prune(wf, i)]
        pass
