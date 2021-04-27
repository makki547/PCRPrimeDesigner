import re
import math

def translate_from_dna_to_aa(dna_seq, frame = 0):

    codons = {'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L',
              'ctt': 'L', 'ctc': 'L', 'cta': 'L', 'ctg': 'L',
              'att': 'I', 'atc': 'I', 'ata': 'I', 'atg': 'M',
              'gtt': 'V', 'gtc': 'V', 'gta': 'V', 'gtg': 'V',
              'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
              'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P',
              'act': 'T', 'acc': 'T', 'aca': 'T', 'acg': 'T',
              'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A',
              'tat': 'Y', 'tac': 'Y', 'taa': '*', 'tag': '*',
              'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q',
              'aat': 'N', 'aac': 'N', 'aaa': 'K', 'aag': 'K',
              'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E',
              'tgt': 'C', 'tgc': 'C', 'tga': '*', 'tgg': 'W',
              'cgt': 'R', 'cgc': 'R', 'cga': 'R', 'cgg': 'R',
              'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R',
              'ggt': 'G', 'ggc': 'G', 'gga': 'A', 'ggg': 'G'}
    
    lower_nuc = dna_seq.lower()
    len_nuc = len(lower_nuc)

    aa_seq = ''
    for codon_index in range(frame, len(lower_nuc), 3):

        if codon_index + 3 > len_nuc: break
        codon = lower_nuc[codon_index:codon_index + 3]

        aa_seq += codons.get(codon, '-')
        
    return aa_seq

def convert_char_case_by_type(sequence, sequence_type):
    
    dna = re.compile(r'[atcgATCG]+')
    if (sequence_type == 'dna') or (sequence_type == 'DNA'):
        sequence_type = 'dna'
        sequence = sequence.lower()

    elif (sequence_type == 'protein') or (sequence_type == 'PROTEIN') or (sequence_type == 'aa') or (sequence_type == 'AA'):

        sequence_type = 'protein'
        sequence = sequence.upper()
        
    else:
        if dna.fullmatch(sequence):
            sequence = sequence.lower()
            sequence_type = 'dna'
        else:
            sequence = sequence.upper()
            sequence_type = 'protein'

    return sequence_type, sequence

def read_fasta(fin, sequence_type = 'auto'):

    comment = re.compile(r'^>')


    
    seq = ''
    for line in fin:

        if comment.match(line):
            continue
        
        seq += line.strip()

    
    return convert_char_case_by_type(seq, sequence_type)

def get_complement_sequence(sequence):

    conv = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(list(map(lambda x: conv[x],list(sequence[::-1]))))    
    

class QueryError(Exception):

    """Query exception for finding unique block."""
    pass

class QueryIsNotFoundError(QueryError):

    """Query is not found in target sequence."""
    pass

class QueryIsNotUniqueError(QueryError):

    """Query is not unique."""
    pass

def get_sequence_position(target, query):

    nmatch = target.count(query)

    if nmatch == 1:

        return target.index(query)
    
    elif nmatch > 1:
        raise QueryIsNotUniqueError()
    else:
        raise QueryIsNotFoundError()        


class PCRPrimer:

    def __init__(self, sequence):

        self.sequence = ''.join(sequence.lower().split())
        pass

    def get_complement(self):

        
        return get_complement_sequence(self.sequence)

    def get_gc_content(self):

        return self.sequence.count('c') + self.sequence.count('g')


    def get_gc_population(self):

        return float(self.get_gc_content()) / float(len(self.sequence))
    
    def is_end_with_gc(self):


        return self.sequence[-1] == 'c' or self.sequence[-1] == 'g'


    def get_tm(self, method = 'wallace', additional_params = None):


        methods = {'wallace': self.get_tm_by_wallace,
                   'nn': self.get_tm_by_nearest_neighbors,
                   'nearest_neighbors': self.get_tm_by_nearest_neighbors}


        return methods.get(method, self.get_tm_by_wallace)(additional_params)


    def get_tm_by_wallace(self, additional_params = None):

        ngc = self.get_gc_content()
        nat = len(self.sequence) - ngc
        
        return float(2*nat + 4*ngc + 35 - 2*len(self.sequence))


    def get_tm_by_nearest_neighbors(self, additional_params):

        def get_di_na():

            for i in range(len(self.sequence) - 1):

                yield self.sequence[i:i+2]

        param = {
            'aa':[-9.1, -24],
            'tt':[-9.1, -24],
            'at':[-8.6, -23.9],
            'ta':[-6, -16.9],
            'ca':[-5.8, -12.9],
            'tg':[-5.8, -12.9],
            'gt':[-6.5, -17.3],
            'ac':[-6.5, -17.3],
            'ct':[-7.8, -20.8],
            'ag':[-7.8, -20.8],
            'ga':[-5.6, -13.5],
            'tc':[-5.6, -13.5],
            'cg':[-11.9, -27.8],
            'gc':[-11.1, -26.7],
            'gg':[-11, -26.6],
            'cc':[-11, -26.6]
        }                

        primer_conc = 2.5E-7
        na_conc = 5.0E-2
        if type(additional_params) is dict:

            primer_conc = float(additional_params.get('primer_conc', primer_conc))
            na_conc = float(additional_params.get('na_conc', na_conc))



        
        ds = 0.0
        dh = 0.0

        for di in get_di_na():

            p = param[di]
            dh += p[0]
            ds += p[1]

              
        return 1000.0 * dh / (-10.8 + ds + 1.987 * math.log(primer_conc/4.0)) - 273.15 + 16.6 * math.log10(na_conc)

    def __add__(self, other):

        if isinstance(other, str):

            return PCRPrimer(self.sequence + other)
        elif isinstance(other, PCRPrimer):
            return PCRPrimer(self.sequence + other.sequence)

        else:
            raise NotImplementedError()

    def __iadd__(self, other):

        if isinstance(other, str):

            self.sequence += ''.join(other.lower().split())

        elif isinstance(other, PCRPrimer):
            self.sequence += other.sequence

        else:
            raise NotImplementedError()        

        return self


    def __len__(self):

        return len(self.sequence)

    def __str__(self):

        return self.sequence
