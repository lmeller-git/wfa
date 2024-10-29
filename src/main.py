from src import align, open_fasta, args, plot_align, compare


def main(args):
    query = open_fasta(args.query)
    db = open_fasta(args.db)
    if args.compare:
        compare(query, db, args)
        return
    if args.plot:
        plot_align(query, db, args)
        return
    align(query, db, args)


if __name__ == "__main__":
    args_ = args()
    main(args_)
