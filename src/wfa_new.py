from Bio import SeqIO, SeqRecord, Seq
from argparse import Namespace
from dataclasses import dataclass
from typing import Self
from src.utils import timeit
from enum import StrEnum

"""
implementation based on x and y coordinates with move vectors using edit distance.
TODO maybe remove x, y and moves, optimize
rewrite this mess to allow for better local and semi global implementations
OR: abandon and reerite in rust
"""
type Records = SeqIO.FastIO.FastaIterator
type Record = SeqRecord.SeqRecord
type Sequence = Seq.Seq


@timeit
def align(query: Records, db: Records, args: Namespace) -> None:
    global scheme
    scheme = ScoringScheme(0, 4, 0, 4)
    for q in query:
        for d in db:
            print(q.seq)
            print(d.seq)
            ocean = Ocean(q, d, args)
            while not ocean.step():
                ocean.prune()
                pass


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
    def __init__(
        self,
        x: int,
        y: int,
        score: int,
        query: Sequence,
        db: Sequence = None,
        all_moves: list[Move] = [],
        diag: int = 0,
        offset: int = 0,
        *args,
        **kwargs,
    ) -> None:
        self.x_ = x
        self.y_ = y
        self.score = score
        self.db = db
        self.query = query
        self.all_moves = all_moves
        self.diag = diag
        self.offset = offset
        self.local_in_q = False
        assert self.get_fr() == (
            self.diag,
            self.offset,
        ), f"{
            self.get_fr()}, {(self.diag, self.offset)}, {self.x_}, {self.y_}"

    def x(self) -> int:
        return self.offset - min(self.diag, 0)

    def y(self) -> int:
        return self.offset + max(self.diag, 0)

    def expand(self) -> bool:
        while (self.x_ < len(self.db) and self.y_ < len(self.query)) and (
            self.db[self.x_] == self.query[self.y_]
        ):
            assert self.get_fr() == (
                self.diag,
                self.offset,
            ), f"{
                self.get_fr()}, {(self.diag, self.offset)}, {self.x_}, {self.y_}"
            assert self.x_ == self.offset - min(
                0, self.diag
            ), f"{self.x_}, {
                    self.offset - min(0, self.diag - 1)}, {self.diag}, {self.offset}"
            assert self.y_ == self.offset + max(
                0, self.diag
            ), f"{self.y_}, {
                    self.offset + max(0, self.diag)}, {self.diag}, {self.offset}"
            self.x_ += 1
            self.y_ += 1
            self.offset += 1
            self.score += scheme.match_
            self.all_moves.append(Move.MatchMismatch)
            if self.is_converged():
                return True
        return False

    def is_converged(self) -> bool:
        return self.x_ == len(self.db) and self.y_ == len(self.query)

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
                    q = "-" + q
                    in_db += 1
                    pass
                case Move.Insertion:
                    db = "-" + db
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

        def item(x, y):
            return " " if x != y else "|"

        print(q)
        [print(item(q_, d_), end="") for q_, d_ in zip(q, db)]
        print()
        print(db)
        print(self.score)
        print("actual score: ", get_score(q, db))

    def get_fr(self) -> tuple[int, int]:
        diag = self.y_ - self.x_
        offset = min(self.x_, self.y_)
        return diag, offset

    def __repr__(self) -> str:
        return f"diag: {self.diag} offset: {self.offset} x: {self.x_} y: {self.y_}"


class SemiWaveFront(WaveFront):
    def __init__(
        self,
        x: int,
        y: int,
        score: int,
        query: Sequence,
        db: Sequence,
        moves: list[Move],
        diag: int = 0,
        offset: int = 0,
        local_in_q: bool = True,
    ) -> None:
        super().__init__(x, y, score, query, db, moves, diag, offset)
        self.local_in_q = local_in_q

    def expand(self) -> bool:
        while (self.x_ < len(self.db) and self.y_ < len(self.query)) and (
            self.db[self.x_] == self.query[self.y_]
        ):
            self.x_ += 1
            self.y_ += 1
            self.offset += 1
            self.score += scheme.match_
            self.all_moves.append(Move.MatchMismatch)
            if self.is_converged():
                return True
            # continue
            if self.y() == len(self.query) and self.local_in_q:
                self.fill_last_row()
                return True
            if self.x() == len(self.db) and not self.local_in_q:
                self.fill_last_col()
                return True
        return False

    def fill_last_row(self) -> None:
        self.all_moves += [Move.Deletion for i in range(len(self.db) - self.x())]
        self.x_ = len(self.db)
        self.diag, self.offset = self.get_fr()

    def fill_last_col(self) -> None:
        self.all_moves += [Move.Insertion for i in range(len(self.query) - self.y())]
        self.y_ = len(self.query)
        self.diag, self.offset = self.get_fr()

    def pprint(self) -> None:
        super().pprint()


class LocalWaveFront(WaveFront):
    def __init__(
        self,
        x: int,
        y: int,
        query: Sequence,
        db: Sequence,
        diag: int,
        offset: int,
        score: int,
        moves: list[Move],
    ) -> None:
        super.__init__(x, y, score, query, db, moves, diag, offset)


def fill_diff(
    parent: WaveFront,
    child: WaveFront,
    q: Sequence,
    db: Sequence,
    q_str: str,
    db_str: str,
):
    q_ = (
        q[parent.offset + max(parent.diag, 0) : child.offset + max(child.diag, 0)]
        + q_str
    )
    db_ = (
        db[parent.offset - min(parent.diag, 0) : child.offset - min(child.diag, 0)]
        + db_str
    )
    return q_, db_


def is_invalid_insertion(parent: WaveFront, child: WaveFront) -> bool:
    if child.diag <= 0:
        return parent.offset >= child.offset
    if parent.offset > child.offset:
        return True
    if parent.db[parent.x() : child.x()] != parent.query[parent.y() + 1 : child.y()]:
        return True
    return False


def is_invalid_deletion(parent: WaveFront, child: WaveFront) -> bool:
    if child.diag >= 0:
        return parent.offset >= child.offset
    if parent.offset > child.offset:
        return True
    if parent.db[parent.x() + 1 : child.x()] != parent.query[parent.y() : child.y()]:
        return True
    return False


def is_invalid_diagonal(parent: WaveFront, child: WaveFront) -> bool:
    return (
        parent.offset >= child.offset
        or parent.diag != child.diag
        or parent.db[parent.x() + 1 : child.x()]
        != parent.query[parent.y() + 1 : child.y()]
    )


def get_score(seq1: Sequence, seq2: Sequence) -> int:
    matches = sum([1 for q, d in zip(seq1, seq2) if seq1 == seq2 and seq1 != "-"])
    mismatches = sum(
        [1 for q, d in zip(seq1, seq2) if q != d and q != "-" and d != "-"]
    )
    gaps = seq1.count("-") + seq2.count("-")
    return (
        scheme.match_ * matches
        + scheme.gap_extension * gaps
        + mismatches * scheme.mismatch
    )


class Ocean:
    def __init__(self, query: Record, db: Record, args: Namespace) -> None:
        self.query = query
        self.db = db
        self.mode = args.mode
        self.verbose = args.verbose
        self.current_score = 0
        match args.mode:
            case "global":
                self.wavetype = WaveFront
            case "semi_global":
                self.wavetype = SemiWaveFront
            case "local":
                self.wavetype = LocalWaveFront
            case _:
                raise Exception("unknown modus")
        self.wavefronts = {
            0: [self.wavetype(0, 0, 0, query.seq, db.seq, [Move.NoMove])],
            scheme.gap_extension: [
                self.wavetype(
                    0,
                    1,
                    scheme.gap_extension,
                    query.seq,
                    db.seq,
                    [Move.NoMove, Move.Insertion],
                    1,
                ),
                self.wavetype(
                    1,
                    0,
                    scheme.gap_extension,
                    query.seq,
                    db.seq,
                    [Move.NoMove, Move.Deletion],
                    -1,
                ),
            ],
        }
        if args.mode == "semi_global":
            self.wavefronts[0] += [
                self.wavetype(
                    j,
                    0,
                    0,
                    query.seq,
                    db.seq,
                    [Move.NoMove] + [Move.Deletion for i in range(j)],
                    -j,
                    local_in_q=True,
                )
                for j in range(1, len(self.db.seq))
            ]
            self.wavefronts[0] += [
                self.wavetype(
                    0,
                    j,
                    0,
                    query.seq,
                    db.seq,
                    [Move.NoMove] + [Move.Insertion for i in range(j)],
                    j,
                    local_in_q=False,
                )
                for j in range(1, len(self.db.seq))
            ]

    def step(self) -> bool:
        found = False
        if self.current_score in self.wavefronts.keys():
            for wf in self.wavefronts[self.current_score]:
                if wf.expand() or wf.is_converged():
                    wf.pprint()
                    print(len(self.wavefronts[self.current_score]))
                    # the next to lines are this way for comparison to n_w
                    # self.backtrace(wf, "", "", self.current_score)
                    return True
                    found = True
        if found:
            return True
        self.new_wave()
        return False

    def new_wave(self) -> None:
        wfs = []
        self.current_score += 1
        mscore = self.current_score - scheme.mismatch
        gscore = self.current_score - scheme.gap_extension
        if mscore in self.wavefronts:
            for wf in self.wavefronts[mscore]:
                wfs.append(
                    self.wavetype(
                        wf.x_ + 1,
                        wf.y_ + 1,
                        self.current_score,
                        self.query.seq,
                        self.db.seq,
                        wf.all_moves + [Move.MatchMismatch],
                        wf.diag,
                        wf.offset + 1,
                    )
                )
        if gscore in self.wavefronts:
            for wf in self.wavefronts[gscore]:
                wfs += [
                    self.wavetype(
                        wf.x_ + 1,
                        wf.y_,
                        self.current_score,
                        self.query.seq,
                        self.db.seq,
                        wf.all_moves + [Move.Deletion],
                        wf.diag - 1,
                        wf.offset + 1 if wf.diag > 0 else wf.offset,
                        local_in_q=wf.local_in_q,
                    ),
                    self.wavetype(
                        wf.x_,
                        wf.y_ + 1,
                        self.current_score,
                        self.query.seq,
                        self.db.seq,
                        wf.all_moves + [Move.Insertion],
                        wf.diag + 1,
                        wf.offset + 1 if wf.diag < 0 else wf.offset,
                        local_in_q=wf.local_in_q,
                    ),
                ]
        if len(wfs) > 0:
            if self.current_score in self.wavefronts:
                self.wavefronts[self.current_score] += wfs
            else:
                self.wavefronts[self.current_score] = wfs

    def can_prune(self, wf: WaveFront, i: int) -> bool:
        if wf.x_ > len(self.db.seq) or wf.y_ > len(self.query.seq):
            return True
        for j, wave in enumerate(self.wavefronts[self.current_score]):
            if j == i:
                continue
            if wave.diag == wf.diag and wave.offset > wf.offset:
                return True
        return False

    def prune(self) -> None:
        if self.current_score not in self.wavefronts:
            return
        self.wavefronts[self.current_score] = [
            wf
            for i, wf in enumerate(self.wavefronts[self.current_score])
            if not self.can_prune(wf, i)
        ]
        pass

    def backtrace(self, wave: WaveFront, q: str, db: str, score: int) -> None:
        # TODO there seems to be the possibility for an eternal loop. also in some conditions,
        # false parents are chosen or the correct parents are not chosen
        if score == 0:
            q_, db_ = fill_diff(
                self.wavetype(
                    0, 0, 0, self.query.seq, self.db.seq, [Move.NoMove], 0, 0
                ),
                wave,
                self.query.seq,
                self.db.seq,
                q,
                db,
            )
            if get_score(q_, db_) != self.current_score:
                return
            print(f"found alignemnt with score {self.current_score}")
            print(f"query: {q_}")
            print("       ", end="")
            [print(" " if q__ != d__ else "|", end="") for q__, d__ in zip(q_, db_)]
            print(f"\ndb:    {db_}")
            print(f"actual score: {get_score(q_, db_)}")
            print("\n")
            return

        gscore = score - scheme.gap_extension
        mscore = score - scheme.mismatch
        if gscore in self.wavefronts:
            for parent in self.wavefronts[gscore]:
                if parent.diag == wave.diag - 1:
                    # Insertion
                    if is_invalid_insertion(parent, wave):
                        continue
                    q_, db_ = fill_diff(
                        parent, wave, self.query.seq, self.db.seq, q, db
                    )
                    self.backtrace(parent, q_, "-" + db_, gscore)

                if parent.diag == wave.diag + 1:
                    # Deletion
                    if is_invalid_deletion(parent, wave):
                        continue
                    q_, db_ = fill_diff(
                        parent, wave, self.query.seq, self.db.seq, q, db
                    )
                    self.backtrace(parent, "-" + q_, db_, gscore)

        if mscore in self.wavefronts:
            for parent in self.wavefronts[mscore]:
                if is_invalid_diagonal(parent, wave):
                    continue
                q_, db_ = fill_diff(parent, wave, self.query.seq, self.db.seq, q, db)
                self.backtrace(parent, q_, db_, mscore)
