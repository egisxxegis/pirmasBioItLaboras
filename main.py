from Bio import SeqIO
from Bio.Seq import Seq


def split_into_reading_frames(seq_obj):
    # we get Seq(AAAGGGTTT)
    # we return [Seq(AAAGGGTTT), Seq(AAGGGT), Seq(AGGGTT)]
    the_split = [seq_obj[x:(len(seq_obj)-x)//3*3+x] for x in range(3)]
    return the_split


def find_start_ends(translated_seq_obj):
    # we get Seq(Q*MQMM*M**Q)
    # we return [(2,6), (4,6), (5,6), (7,8)] - pairs of M and * pos
    the_start_ends = []
    the_starts = []
    the_ends = []
    if len(translated_seq_obj) < 1:
        return the_start_ends

    # find all starts ant ends
    for position in range(len(translated_seq_obj)):
        if translated_seq_obj[position] == 'M':
            the_starts.append(position)
        elif translated_seq_obj[position] == '*':
            the_ends.append(position)
    if len(the_starts) < 1 or len(the_ends) < 1:
        return the_start_ends

    # make pairs
    the_index_to_check_start = 0
    the_index_to_check_end = 0
    while True:
        start = the_starts[the_index_to_check_start]
        end = the_ends[the_index_to_check_end]
        if end < start:
            # take another end
            the_index_to_check_end += 1
            if the_index_to_check_end >= len(the_ends):
                break
            continue

        the_start_ends.append((start, end))
        the_index_to_check_start += 1
        if the_index_to_check_start >= len(the_starts):
            break
        continue

    return the_start_ends


def filter_start_end_to_longest_pairs(start_stops):
    # we get [(2,6), (4,6), (5,6), (7,8)]
    # we return [(2,6), (7, 8)]
    the_last_checked_end = 0
    the_start_ends = []
    for start, stop in start_stops:
        if stop == the_last_checked_end:
            continue
        the_start_ends.append((start, stop))
        the_last_checked_end = stop
    return the_start_ends


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
        start_ends = [find_start_ends(frame.translate()) for frame in main_frames]
        reverse_start_ends = [find_start_ends(frame.translate()) for frame in reverse_complement_frames]
        # task 1 done.
        start_ends_no_overlap = [filter_start_end_to_longest_pairs(pairs) for pairs in start_ends]
        reverse_start_ends_no_overlap = [filter_start_end_to_longest_pairs(pairs) for pairs in reverse_start_ends]
        # task 2 done.


else:
    print(f'Execution cancelled, not the main.py called')
