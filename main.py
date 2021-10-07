from Bio import SeqIO
from Bio.Seq import Seq
from math import ceil
import numpy as np
from tabulate import tabulate


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


def print_data_in_table(data, top_down_headers=None, left_right_headers=None, dimensions=1):
    if left_right_headers is None:
        left_right_headers = ["data"]
    if dimensions < 1 or dimensions > 2:
        print(f'{dimensions} dimensions are not supported in print_data_in_table.')
        return
    if dimensions == 1:
        if top_down_headers is None:
            print(tabulate(data, left_right_headers))
        else:
            new_data = []
            for index in range(len(data)):
                new_data.append([top_down_headers[index], data[index]])
            print(tabulate(new_data, ["name"] + left_right_headers))
        return
    elif dimensions == 2:
        # we substract one because top_down_headers will be added
        if len(data) % (len(left_right_headers)-1) > 0:
            print(f'Given data cannot be converted to {len(left_right_headers)} column table.')
            return
        elif len(left_right_headers) == 1:
            print(f'1 column header can not form 2 dimensional table.')
            return
        if top_down_headers is None:
            print(f'top down headers are must for 2 dimensional table. Else use 1 dimensional')
            return
        new_data = []
        the_index = 0
        the_index_of_header = 0
        the_length = len(left_right_headers) - 1
        while the_index + the_length < len(data)+1:
            new_data.append(
                [top_down_headers[the_index_of_header]]
                + data[the_index: the_index + the_length])
            the_index += the_length
            the_index_of_header += 1
        print(tabulate(new_data, left_right_headers))
        return


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
        all_frames_translated = [frame.translate() for frame in all_frames]
        all_frames_start_ends = [find_start_ends(translated_frame) for translated_frame in all_frames_translated]
        # task 1 done above.
        all_frames_start_ends_no_overlap = [filter_start_end_to_longest_pairs(pairs) for pairs in all_frames_start_ends]
        # task 2 done above.
        minimum_bp = 100  # three bp codes one amino acid
        fragments = [extract_fragments(all_frames_translated[i],
                                       all_frames_start_ends_no_overlap[i],
                                       ceil(minimum_bp / 3))
                     for i in range(len(all_frames_translated))]
        # task 3 done above.
        translated_nucleotides = "ARNDCEQGHILKMFPSTWYV"
        translated_dicodones = [codone + codone2
                                for codone in translated_nucleotides
                                for codone2 in translated_nucleotides]
        frequency_of_codones = count_frequency(all_frames_translated, translated_nucleotides)
        frequency_of_dicodones = count_frequency(all_frames_translated, translated_dicodones)
        print_data_in_table(frequency_of_codones.tolist(),
                            top_down_headers=translated_nucleotides,
                            left_right_headers=["Codone", "Frequency"],
                            dimensions=2)
        print_data_in_table(frequency_of_dicodones.tolist(),
                            top_down_headers=translated_nucleotides,
                            left_right_headers=["name"] + [*translated_nucleotides],
                            dimensions=2)
        # task 4 done above.


else:
    print(f'Execution cancelled, not the main.py called')
