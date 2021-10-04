from main import split_into_reading_frames, find_start_ends

the_split = split_into_reading_frames("AAAGGGTTT")
if len(the_split) != 3 or\
        the_split[0] != "AAAGGGTTT" or\
        the_split[1] != "AAGGGT" or\
        the_split[2] != "AGGGTT":
    print("split into reading frames failed.")

the_start_ends = find_start_ends("Q*MQMM*M**Q")
if the_start_ends[0] != (2, 6) or\
        the_start_ends[1] != (4, 6) or\
        the_start_ends[2] != (5, 6) or\
        the_start_ends[3] != (7, 8):
    print("find start ends failed.")

print("tests done")
