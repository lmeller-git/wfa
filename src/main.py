from src import align, open_fasta, args


def main(args):
    query = open_fasta(args.query)
    db = open_fasta(args.db)
    align(query, db, args)


if __name__ == "__main__":
    args_ = args()
    main(args_)
