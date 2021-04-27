import genomics_tools
import argparse

def parse_args():

    parser = argparse.ArgumentParser(description = 'Get sequece position of unique motif from DNA seq.')

    parser.add_argument('target_sequence_file', type = argparse.FileType('r'))

    parser.add_argument('query', type = str)

    parser.add_argument('-t', '--type', type = str, choices = ['auto', 'dna', 'protein'], default = 'auto', dest = 'sequence_type')
    parser.add_argument('-r', '--right',action = 'store_true', dest = 'get_right_index')
    
    return parser.parse_args()

def pattern_search_dna_vs_dna(target_seq, query_seq):

    try:
        position = genomics_tools.get_sequence_position(target_seq, query_seq)
    except genomics_tools.QueryIsNotFoundError:
        print('Query sequence is not found')
        exit()

    except genomics_tools.QueryIsNotUniqueError:
        print('Query sequence is not unique')
        exit()
    except:
        print('Unknown errror')
        exit()

    return position


def pattern_search_dna_vs_protein(target_seq, query_seq):

    seq_positions = []

    for frame in range(3):

        target_aa = genomics_tools.translate_from_dna_to_aa(target_seq, frame)


        try:
            position = genomics_tools.get_sequence_position(target_aa, query_seq)
            seq_positions.append(position)
            
        except genomics_tools.QueryIsNotFoundError:
            seq_positions.append(None)


        except genomics_tools.QueryIsNotUniqueError:
            seq_positions.append(None)            

        except:
            print('Unknown errror')
            exit()            



    num_unmatched = seq_positions.count(None)
    if num_unmatched != 2:
        print('The query is not unique or not found.')
        exit()

    for frame, position in enumerate(seq_positions):

        if position is not None:
            return frame, position
    
    return None, None

        


if __name__ == '__main__':

    parsed_args = parse_args()

    target_type, target_seq = genomics_tools.read_fasta(parsed_args.target_sequence_file)

    if target_type != 'dna':

        print('Target file is not DNA sequence. Abort.')
        exit()

    
    query_type, query_seq = genomics_tools.convert_char_case_by_type(parsed_args.query, parsed_args.sequence_type)

    if query_type == 'dna':

        
        position = pattern_search_dna_vs_dna(target_seq, query_seq)
        query_length = len(query_seq)

    
    elif query_type == 'protein':

        frame, position = pattern_search_dna_vs_protein(target_seq, query_seq)

        position = position * 3 + frame
        query_length = 3*len(query_seq)

    else:
        exit()

    if parsed_args.get_right_index:

        position = position + query_length - 1

    print('Position = %d' % position)
