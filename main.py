from Bio import SeqIO
from Bio.Seq import Seq


def split_into_reading_frames(seq_obj):
    # we get Seq(AAAGGGTTT)
    # we return [Seq(AAAGGGTTT), Seq(AAGGGT), Seq(AGGGTT)]
    the_split = [seq_obj[x:(len(seq_obj)-x)//3*3+x] for x in range(3)]
    return the_split


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
        main_frames = split_into_reading_frames(seq_record)
        reverse_complement_frames = split_into_reading_frames(seq_record.reverse_complement())

else:
    print(f'Execution cancelled, not the main.py called')
