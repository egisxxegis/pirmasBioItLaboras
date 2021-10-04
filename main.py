from Bio import SeqIO, Seq


def split_into_reading_frames(seq_obj, offset=0):
    return 420


def find_start_end(something):
    return 6, 9


if __name__ == '__main__':
    bacterials = [f'sources\\data\\bacterial{x+1}.fasta' for x in range(4)]
    mamalians = [f'sources\\data\\mamalian{x+1}.fasta' for x in range(4)]

    for seq_record in SeqIO.parse(bacterials[2], "fasta"):
        print(seq_record.id)
        print(f'original:   {repr(seq_record.seq)}')
        print(f'translated: {repr(seq_record.seq.translate())}')
        print(f'reversed:   {repr(seq_record.seq.reverse_complement())}')
        print(f'symbols:    {len(seq_record)}')

else:
    print(f'Execution cancelled, not the main.py called')
