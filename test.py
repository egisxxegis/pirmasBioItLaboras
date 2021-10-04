from main import split_into_reading_frames, find_start_ends, filter_start_end_to_longest_pairs

the_split = split_into_reading_frames("AAAGGGTTT")
if the_split != ["AAAGGGTTT", "AAGGGT", "AGGGTT"]:
    print("---split into reading frames failed.")

the_start_ends = find_start_ends("Q*MQMM*M**Q")
if the_start_ends != [(2, 6), (4, 6), (5, 6), (7, 8)]:
    print("---find start ends failed.")

the_start_ends_filtered = filter_start_end_to_longest_pairs([(2, 6), (4, 6), (5, 6), (7, 8)])
if the_start_ends_filtered != [(2, 6), (7, 8)]:
    print("---filter start end to longest pairs failed.")
print("tests done")
