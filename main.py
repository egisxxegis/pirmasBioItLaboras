from Bio import SeqIO
from Bio.Seq import Seq
from math import ceil
import numpy as np


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


def extract_fragments(sequence_object, start_stops, minimum_length=1):
    # we get Seq(Q*MQMM*M**Q), [(2, 6), (7, 8)], 2
    # we return [Seq(MQMM), Seq(M)]
    the_fragments = []
    for start, stop in start_stops:
        if stop - start < minimum_length:
            continue
        the_fragment = sequence_object[start:stop]
        the_fragments.append(the_fragment)
    return the_fragments


def occurrences(the_haystack, the_needle):
    the_count = 0
    the_start = 0
    while True:
        the_start = the_haystack.find(the_needle, the_start) + 1
        if the_start > 0:
            the_count += 1
        else:
            return the_count


def count_frequency(arr_where, arr_what):
    the_total_frequency = np.zeros(len(arr_what))
    for seq in arr_where:
        the_length_of_seq = len(seq) - len(arr_what[0]) + 1
        for i in range(len(arr_what)):
            the_total_frequency[i] += occurrences(seq, arr_what[i]) / the_length_of_seq
    return the_total_frequency / len(arr_where)


if __name__ == '__main__':
    bacterials = [f'sources\\data\\bacterial{x+1}.fasta' for x in range(4)]
    mamalians = [f'sources\\data\\mamalian{x+1}.fasta' for x in range(4)]

    for seq_record in SeqIO.parse(bacterials[2], "fasta"):
        print(seq_record.id)
        print(f'original:   {repr(seq_record.seq)}')
        print(f'translated: {repr(seq_record.seq.translate())}')
        print(f'reversed:   {repr(seq_record.seq.reverse_complement())}')
        print(f'symbols:    {len(seq_record)}')
        main_frames = split_into_reading_frames(seq_record.seq)
        reverse_complement_frames = split_into_reading_frames(seq_record.seq.reverse_complement())
        all_frames = main_frames + reverse_complement_frames
        all_frames_start_ends = [find_start_ends(frame.translate()) for frame in all_frames]
        # task 1 done above.
        all_frames_start_ends_no_overlap = [filter_start_end_to_longest_pairs(pairs) for pairs in all_frames_start_ends]
        # task 2 done above.
        minimum_bp = 100  # three bp codes one amino acid
        fragments = [extract_fragments(all_frames[i], all_frames_start_ends_no_overlap[i], ceil(minimum_bp / 3))
                     for i in range(len(all_frames))]
        # task 3 done above.
        translated_nucleotides = "ARNDCEQGHILKMFPSTWYV"
        translated_dicodones = [codone + codone2
                                for codone in translated_nucleotides
                                for codone2 in translated_nucleotides]
        frequency_of_codones = count_frequency(all_frames, translated_nucleotides)
        frequency_of_dicodones = count_frequency(all_frames, translated_dicodones)


else:
    print(f'Execution cancelled, not the main.py called')
