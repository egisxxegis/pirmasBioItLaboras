import numpy as np
from main import \
    split_into_reading_frames, \
    find_start_ends, \
    filter_start_end_to_longest_pairs, \
    extract_fragments, \
    occurrences, \
    count_frequency, \
    format_numpy_array_digits, \
    calculate_difference_frequencies

the_split = split_into_reading_frames("AAAGGGTTT")
if the_split != ["AAAGGGTTT", "AAGGGT", "AGGGTT"]:
    print("---split into reading frames failed.")

the_start_ends = find_start_ends("Q*MQMM*M**Q")
if the_start_ends != [(2, 6), (4, 6), (5, 6), (7, 8)]:
    print("---find start ends failed.")

the_start_ends_filtered = filter_start_end_to_longest_pairs([(2, 6), (4, 6), (5, 6), (7, 8)])
if the_start_ends_filtered != [(2, 6), (7, 8)]:
    print("---filter start end to longest pairs failed.")

the_fragments = extract_fragments("Q*MQMM*M**Q", [(2, 6), (7, 8)], 2)
if the_fragments != ["MQMM"]:
    print("---extract fragments failed.")

the_occurences = occurrences("AAABAAA", "AA")
if the_occurences != 4:
    print("---The occurences failed.")

the_where1 = ["Q*MQMM*M**Q", "MAAQQA*QQ", "QQ*TA*AAMM*"]
the_what1 = "QMAR"
the_where2 = ["Q*MQMM*M**Q", "MAAAAA*AM"]
the_what2 = ["QQ", "QM", "QA", "MQ", "MM", "MA", "AQ", "AM", "AA"]
the_frequency = count_frequency(the_where1, the_what1)
the_answer = np.array([0.2996632996632997, 0.21885521885521886, 0.202020202020202, 0.0])
if not (the_frequency == the_answer).all():
    print("---frequency with one letters failed.")
the_frequency = count_frequency(the_where2, the_what2)
the_answer = np.array([0., 0.05, 0., 0.05, 0.05, 0.0625, 0., 0.0625, 0.25])
if not (the_frequency == the_answer).all():
    print("---frequency with two letters failed.")

the_numpy_array = np.array([0.0, 12.5, 0.00231, 0.00235, 0.12345, 0.00001, 0.00005])
the_formatted_numpy_array = format_numpy_array_digits(the_numpy_array, 4, True, True)
the_answer = np.array([0.0, 12.5, 0.0023, 0.0024, 0.1235, 0.0001, 0.0001])
if not (the_formatted_numpy_array == the_answer).all():
    print("---numpy array formatting with near zeros and rounding failed.")
the_formatted_numpy_array = format_numpy_array_digits(the_numpy_array, 4, True, False)
the_answer = np.array([0.0, 12.5, 0.0023, 0.0023, 0.1234, 0.0001, 0.0001])
if not (the_formatted_numpy_array == the_answer).all():
    print("---numpy array formatting with near zeros and NO rounding failed.")
the_formatted_numpy_array = format_numpy_array_digits(the_numpy_array, 4, False, True)
the_answer = np.array([0.0, 12.5, 0.0023, 0.0024, 0.1235, 0, 0.0001])
if not (the_formatted_numpy_array == the_answer).all():
    print("---numpy array formatting with NO near zeros and rounding failed.")

the_score = calculate_difference_frequencies(np.array([0, 1, 2]), np.array([0, 1, 2]))
if the_score != 0:
    print("---calculate difference failed calculating difference of two same arguments.")

print("tests done")
