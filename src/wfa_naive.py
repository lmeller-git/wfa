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
            ocean = Ocean(q, d, args)
            while not ocean.step():
                ocean.prune()
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

    def __init__(self, last_move: Move):
        self.c = 0
        self.last_move = last_move

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> MoveInfo:
        match self.c:
            case 0:
                self.c += 1
                if self.last_move == Move.Deletion:
                    return MoveInfo(Move.Deletion, scheme.gap_extension, 1, 0)
                else:
                    return MoveInfo(Move.Deletion, scheme.gap_opening + scheme.gap_extension, 1, 0)
            case 1:
                self.c += 1
                return MoveInfo(Move.MatchMismatch, scheme.mismatch, 1, 1)
            case 2:
                self.c += 1
                if self.last_move == Move.Insertion:
                    return MoveInfo(Move.Insertion, scheme.gap_extension, 0, 1)
                else:
                    return MoveInfo(Move.Insertion, scheme.gap_opening + scheme.gap_extension, 0, 1)
            case _:
                self.c = 0
                raise StopIteration


class WaveFront:
    def __init__(self, x: int, y: int, score: int, query: Sequence, db: Sequence, last_move: Move, all_moves: list[Move]) -> None:
        self.x = x
        self.y = y
        self.score = score
        self.db = db
        self.query = query
        self.last_move = last_move
        self.all_moves = all_moves

    def expand(self) -> list[Self]:
        while (self.x < len(self.db) and self.y < len(self.query)) and (self.db[self.x] == self.query[self.y]):
            self.x += 1
            self.y += 1
            self.score += scheme.match_
            self.last_move = Move.MatchMismatch
            self.all_moves.append(self.last_move)
        if self.x == len(self.db):
            self.expand_insert()
        elif self.y == len(self.query):
            self.expand_delet()
        wfs = []
        for move_info in MoveIterator(self.last_move):
            if self.x + move_info.dx > len(self.db) or self.y + move_info.dy > len(self.query):
                continue
            wfs.append(WaveFront(self.x + move_info.dx, self.y + move_info.dy,
                       self.score + move_info.dscore, self.query, self.db, move_info.move, self.all_moves + [move_info.move]))
        return wfs

    def is_converged(self) -> bool:
        return self.x == len(self.db) and self.y == len(self.query)

    def expand_insert(self):
        while self.y < len(self.query):
            self.last_move = Move.Insertion
            self.score += scheme.gap_extension
            self.y += 1
            self.all_moves.append(self.last_move)

    def expand_delet(self):
        while self.x < len(self.db):
            self.last_move = Move.Deletion
            self.score += scheme.gap_extension
            self.x += 1
            self.all_moves.append(self.last_move)

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


class Ocean:
    def __init__(self, query: Record, db: Record, args: Namespace) -> None:
        self.query = query
        self.db = db
        self.mode = args.mode
        self.verbose = args.verbose
        self.wavefronts = [
            WaveFront(0, 0, 0, query.seq, db.seq, Move.NoMove, [Move.NoMove])]

    def backtrace(self) -> None:
        pass

    def step(self) -> bool:
        new_wfs = []
        for wf in self.wavefronts:
            new_wfs += wf.expand()
            if wf.is_converged():
                wf.pprint()
                return True

        self.wavefronts += new_wfs
        print("false")
        return False

    def prune(self) -> None:
        pass
