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
