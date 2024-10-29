from src.wfa_new import align
from Bio import SeqIO, SeqRecord, Seq
from argparse import Namespace
import numpy as np
from src.utils import timeit

type Records = SeqIO.FastIO.FastaIterator
type Record = SeqRecord.SeqRecord
type Sequence = Seq.Seq


def compare(query: Records, db: Records, args: Namespace) -> None:
    query = list(query)
    db = list(db)
    run_n_w(query, db)
    align(query, db, args)


@timeit
def run_n_w(query: Records, db: Records) -> None:
    for q in query:
        for d in db:
            n_w(q.seq, d.seq)


def n_w(query: Sequence, db: Sequence) -> None:

    text = db
    pattern = query

    if len(text) < len(pattern):
        p, t = text, pattern
        print(
            "NOTE: Text and pattern were exchanged such that text is larger than pattern"
        )
    else:
        t = text
        p = pattern

    matrix = np.zeros([len(t) + 1, len(p) + 1], dtype=np.int16)

    def print_matrix(t, p, m):
        return
        """
        This function prints a graphic representation of the score matrix m for the global alignment
        of a text t and patter p.
        t, p ... string
        m    ... numpy array
        Note: Matrix m is a 2-dimensional array whereby the first dimension of matrix m is required
        to equal the length text t plus 1, and the second dimension of matrix m is required to equal
        the length of pattern p plus 1.
        """
        print("  | %3s" % " ", end="")
        [print(" %3s" % c, end="") for c in t]
        print()
        print("-" * ((len(t) + 1) * 4 + 3))
        [
            [print("%s |" % [c for c in " " + p][j], end="")]
            + [print(" %3d" % e, end="") for e in m[:, j]]
            + [print()]
            for j in range(len(p) + 1)
        ]
        print()

    # """ <-----(remove the # to use the inverse scoring below
    # calculate scores counting edit operations
    MTCH = 0  # match
    MM = 1  # mismatch
    INS = 1  # INSERTION into referece = vertical progression: MOVE DOWN
    DEL = 1  # DELETION from reference = horizontal progression: MOVE RIGHT
    BEST_SCORE_F = min  # BEST_SCORE_F used as function pointer use () to call min
    """
    # calculate scores with penalties for operations
    MTCH = 0  # match
    MM   = -1  # mismatch
    INS  = -1  # INSERTION into referece = vertical progression: MOVE DOWN
    DEL  = -1  # DELETION from reference = horizontal progression: MOVE RIGHT
    BEST_SCORE_F = max # BEST_SCORE_F() will call the max function
    #"""

    # fill in initial values
    for i in range(len(t) + 1):
        matrix[i, 0] = i * DEL  # initialize first row
    for j in range(1, len(p) + 1):
        matrix[0, j] = j * INS  # initialize first column

    # print_matrix(t,p,matrix)

    # CIGAR (Compact Idiosyncratic Gapped Alignment Report)
    cigar_matrix = np.ndarray([len(t) + 1, len(p) + 1], dtype=object)
    for i in range(len(t) + 1):
        cigar_matrix[i, 0] = [
            "D" * i,
        ]  # initialize first row
    for j in range(1, len(p) + 1):
        cigar_matrix[0, j] = [
            "I" * j,
        ]

    # print(cigar_matrix)

    for j_p in range(len(p)):
        j = j_p + 1  # index for pattern and matrix differ by initialization colum
        for i_t in range(len(t)):
            i = i_t + 1  # index for text and matrix differ by initialization row
            diag_sc = MTCH if t[i_t] == p[j_p] else MM
            diag_sc += matrix[i - 1, j - 1]  # from diagonally left upper cell
            ins_sc = matrix[i, j - 1] + INS  # from cell above
            del_sc = matrix[i - 1, j] + DEL  # from left cell
            score = BEST_SCORE_F(diag_sc, ins_sc, del_sc)
            matrix[i, j] = score  # store the best score

            # calculate list of CIGAR strings for best score(s)
            cigar_list = []
            # """ <-----(remove the # to use the inverse scoring below
            if score == diag_sc:
                mtch_char = "M" if t[i_t] == p[j_p] else "X"
                cigar_list += [cig + mtch_char for cig in cigar_matrix[i - 1, j - 1]]
            if score == ins_sc:
                cigar_list += [cig + "I" for cig in cigar_matrix[i, j - 1]]
            if score == del_sc:
                cigar_list += [cig + "D" for cig in cigar_matrix[i - 1, j]]
            """
            # collect al possible paths
            cigar_list = []
            mtch_char = 'M' if t[i_t]==p[j_p] else 'X'
            cigar_list += [cig+mtch_char for cig in cigar_matrix[i-1,j-1]]
            cigar_list += [cig+'I' for cig in cigar_matrix[i,j-1]]
            cigar_list += [cig+'D' for cig in cigar_matrix[i-1,j]]
            #"""
            cigar_matrix[i, j] = cigar_list

    print_matrix(t, p, matrix)
    # print("CIGAR strings of best alignment(s):\n\t" +
    #     "\n\t".join(cigar_matrix[-1, -1]) + "\n")

    def print_backtrace(t, p, cigar):
        """
        This function prints a graphic representation of gap alignments of a text t and patter p from
        a cigar string.
        t, p, cigar ... string
        Note: cigar needs to hold a valid alignment path which requires that the sum of occurances of
        "M", "X", and "D" equals the length of text t, and the sum of occurances of "M", "X", and "I"
        equals the length of pattern p.
        """
        t_a = ""
        t_i = 0
        p_a = ""
        p_i = 0
        mstr = ""
        for c in cigar:
            if (
                c == "M" or c == "X" or c == "D"
            ):  # white to text alignment string for diagonal and deletion moves
                t_a += t[t_i]
                t_i += 1
            if (
                c == "M" or c == "X" or c == "I"
            ):  # write to pattern alignment string for diagonal and insertion moves
                p_a += p[p_i]
                p_i += 1
            if c == "D":  # for deletion move add a gap to text alignment string
                p_a += "-"
            elif c == "I":  # for insertion move add a gap to pattern alignment string
                t_a += "-"
            # add | to alignment match string for a match else .
            if c == "M":
                mstr += "|"
            elif c == "X" or c == "D" or c == "I":
                mstr += "."
        print("reference: %s" % t_a)
        print("           %s" % mstr)
        print("query    : %s" % p_a)
        print()

    print_backtrace(t, p, cigar_matrix[-1, -1][0])
    # _ = [print_backtrace(t, p, cig) for cig in cigar_matrix[-1, -1]]
    # _ = [print_backtrace(t, p, cig) for cig in cigar_matrix[-1, -1]]
