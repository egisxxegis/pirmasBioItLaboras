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


def format_numpy_array_digits(numpy_array,
                              decimal_spaces,
                              should_near_zero_be_adjusted_to_low_value=True,
                              do_the_round=False):
    the_divisor = 10 ** decimal_spaces
    the_low_limit = 10 ** -decimal_spaces
    the_return = numpy_array * the_divisor // 1 / the_divisor
    if do_the_round:
        the_divisor = the_divisor * 10
        the_return = ((numpy_array * the_divisor % 10) >= 5) \
            * the_low_limit \
            + the_return
    if should_near_zero_be_adjusted_to_low_value:
        the_low_limit_compare = the_low_limit if not do_the_round else the_low_limit / 2
        the_return = (numpy_array < the_low_limit_compare) * (numpy_array > 0) \
            * the_low_limit \
            + the_return
    return the_return


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
        if len(data) % (len(left_right_headers) - 1) > 0:
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


def calculate_difference_frequencies(numpy_array_frequencies1, numpy_array_frequencies2):
    the_freq1 = numpy_array_frequencies1.copy()
    the_freq2 = numpy_array_frequencies2.copy()

    # substitute zero entries
    the_freq1 = (the_freq1 <= 0) * 1e-9 + the_freq1
    the_freq2 = (the_freq2 <= 0) * 1e-9 + the_freq2

    # f1 / f2 only with entries of >1 values
    the_mask = the_freq1 / the_freq2 >= 1
    the_ratio1 = the_mask * the_freq1 / the_freq2
    the_ratio2 = np.invert(the_mask) * the_freq2 / the_freq1

    # rate
    the_ratio = the_ratio1 + the_ratio2
    # limit max
    the_mask = (the_ratio > 4)
    the_ratio = np.invert(the_mask) * the_ratio + the_mask * 4
    # rate values
    # the_thresholds = [3.5, 3.0, 2.6, 2.2, 2.0, 1.8, 1.6, 1.4, 1.3, 1.0]
    the_thresholds = [4.0, 3.0, 2.0, 1.3, 1.0]
    the_max_weight = -1 * len(the_thresholds) - 1
    for the_index in range(len(the_thresholds)):
        the_mask = (the_ratio >= the_thresholds[the_index])
        the_ratio = the_mask * the_ratio * (the_max_weight + the_index) + the_ratio
    the_ratio = the_ratio * -1

    # calculate final rate
    return the_ratio.mean()-1


def print_in_phylip(phylip_array):
    the_row_count = phylip_array[0]
    print(the_row_count)
    the_start = True
    for the_row in phylip_array:
        if the_start:
            the_start = False
            continue
        print(' '.join(the_row))


if __name__ == '__main__':
    bacterials = [f'sources\\data\\bacterial{x+1}.fasta' for x in range(4)]
    mamalians = [f'sources\\data\\mamalian{x+1}.fasta' for x in range(4)]

    all_viruses = bacterials + mamalians
    all_viruses_codones_frequencies = []
    all_viruses_dicodones_frequencies = []

    for one_file_name in all_viruses:
        for seq_record in SeqIO.parse(one_file_name, "fasta"):
            print(seq_record.id)
            print(f'original:   {repr(seq_record.seq)}')
            print(f'symbols:    {len(seq_record)}')
            main_frames = split_into_reading_frames(seq_record.seq)
            reverse_complement_frames = split_into_reading_frames(seq_record.seq.reverse_complement())
            all_frames = main_frames + reverse_complement_frames
            all_frames_translated = [frame.translate() for frame in all_frames]
            all_frames_start_ends = [find_start_ends(translated_frame) for translated_frame in all_frames_translated]
            # task 1 done above.
            all_frames_start_ends_no_overlap = [filter_start_end_to_longest_pairs(pairs)
                                                for pairs in all_frames_start_ends]
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
            print(f'\nthe table of codon frequencies in {seq_record.id}\n')
            print_data_in_table(frequency_of_codones.tolist(),
                                top_down_headers=translated_nucleotides,
                                left_right_headers=["Codone", "Frequency"],
                                dimensions=2)
            print(f'\nthe table of dicodon frequencies in {seq_record.id}\n')
            print_data_in_table(format_numpy_array_digits(frequency_of_dicodones,
                                                          decimal_spaces=4,
                                                          should_near_zero_be_adjusted_to_low_value=True,
                                                          do_the_round=True).tolist(),
                                top_down_headers=translated_nucleotides,
                                left_right_headers=["name"] + [*translated_nucleotides],
                                dimensions=2,)
            # task 4 done above.
            all_viruses_codones_frequencies.append((seq_record.id, frequency_of_codones))
            all_viruses_dicodones_frequencies.append((seq_record.id, frequency_of_dicodones))

    # end for - all records have been read
    # now it is time to compare
    distance_codone_in_phylip = [len(all_viruses_codones_frequencies)]
    for frequency in all_viruses_codones_frequencies:
        frequency_row = [frequency[0]]  # name of sample
        for i in range(len(all_viruses_codones_frequencies)):
            distance = calculate_difference_frequencies(frequency[1],
                                                        all_viruses_codones_frequencies[i][1])
            frequency_row.append(f'{distance:2.5f}')
        distance_codone_in_phylip.append(frequency_row)

    distance_dicodone_in_phylip = [len(all_viruses_dicodones_frequencies)]
    for frequency in all_viruses_dicodones_frequencies:
        frequency_row = [frequency[0]]  # name of sample
        for i in range(len(all_viruses_dicodones_frequencies)):
            distance = calculate_difference_frequencies(frequency[1],
                                                        all_viruses_dicodones_frequencies[i][1])
            frequency_row.append(f'{distance:2.5f}')
        distance_dicodone_in_phylip.append(frequency_row)

    print("\nCodon distance matrix in Phylip")
    print_in_phylip(distance_codone_in_phylip)

    print("\nDicodon distance matrix in Phylip")
    print_in_phylip(distance_dicodone_in_phylip)
    # task 5 done above.

else:
    print(f'Execution cancelled, not the main.py called')
