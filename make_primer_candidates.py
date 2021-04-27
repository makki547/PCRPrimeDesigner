import genomics_tools
import argparse
import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser(description = 'Make primer candidates from DNA sequence and start position')

    parser.add_argument('target_sequence_file', type = argparse.FileType('r'))

    parser.add_argument('position', type = int)

    parser.add_argument('output_file', type = str)    
    parser.add_argument('-c', '--complement', action = 'store_true')

    parser.add_argument('-s', '--shortest', type = int, default = 15, dest = 'shortest_length')
    parser.add_argument('-l', '--longest', type = int, default = 30, dest = 'longest_length')        
    
    return parser.parse_args()



if __name__ == '__main__':

    parsed_args = parse_args()

    target_type, target_seq = genomics_tools.read_fasta(parsed_args.target_sequence_file)

    if target_type != 'dna':

        print('Target file is not DNA sequence. Abort.')
        exit()


    position_begin = parsed_args.position
    
    if parsed_args.complement:

        target_seq = genomics_tools.get_complement_sequence(target_seq)
        position_begin = len(target_seq) - position_begin - 1


    primers = []

    for position_end in range(position_begin + parsed_args.shortest_length,
                              position_begin + parsed_args.longest_length + 1):

        
         seq = target_seq[position_begin:position_end]

         primers.append(genomics_tools.PCRPrimer(seq))


    df_primer_candidates = pd.DataFrame(columns = ['length', 'sequence', 'gc_content', 'tm', 'end_with_gc'])
    for primer in primers:

        row = pd.Series([len(primer),
                         str(primer),
                         primer.get_gc_population(),
                         primer.get_tm('nn', {'primer_conc': 2.5E-7}),
                         primer.is_end_with_gc()],
                        index = df_primer_candidates.columns)

        df_primer_candidates = df_primer_candidates.append(row, ignore_index = True)



    df_primer_candidates.to_csv(parsed_args.output_file, index = False)

         


    
    
