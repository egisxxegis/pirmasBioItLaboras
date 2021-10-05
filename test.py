from main import \
    split_into_reading_frames, \
    find_start_ends, \
    filter_start_end_to_longest_pairs, \
    extract_fragments

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
print("tests done")
