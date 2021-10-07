from Bio import SeqIO
from math import ceil

from sequence_tools import *
from representation_tools import *


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
            print("tai butu baltymų koduojancios sekos")
            for frame_of_fragments in fragments:
                print(frame_of_fragments)
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
