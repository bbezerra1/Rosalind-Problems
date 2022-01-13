# Count the number of each base in the sequence
def countdna():
    y = open('rosalind_dna (1).txt', 'r').read()
    return y.count('A'), y.count('G'), y.count('C'), y.count('T'), y.count('U')


# Translate a DNA sequence into a RNA seq
def transrna(arq):
    y = open(arq, 'r').read()
    return y.replace('T', 'U')


# Create file
def create(content, file_name):
    f = open(file_name, 'w')
    f.write(content)


# Returns the complement strand
def compstrand():
    y = open('rosalind_revc.txt', 'r').read()
    x = y.maketrans('ACTG', 'TGAC')
    return y[::-1].translate(x)


# Read FASTA files
def fastaread(arq):
    data = {}
    seq = ''
    y = open(arq, 'r')

    line = y.readline().rstrip()
    if not line.startswith('>'):
        raise Exception('Records in FASTA files should start with ">" character')
    while line:
        try:
            if line[0] == '>':
                seq = line[1:].rstrip()
                data[seq] = ''
                line = y.readline().rstrip()
            else:
                data[seq] += line
                line = y.readline().rstrip()
        except EOFError:
            y.close()
            break

    return data


# Counts the GC content for the given sequence
def gccont(v):
    n = v[1]
    return 100 * (n.count('G') + n.count('C')) / len(n)


# Iterates through the FASTA dataset to count and find the max GC content
def hgccont(arq):
    f = fastaread(arq)

    # Just to print the GC content for all sequences
    gccontent = []
    print('SEQUENCE ID     GC CONTENT %')
    for item in f.items():
        gccontent.append(gccont(item))
        print(item[0], gccontent[-1])

    # Find and return the maximum value
    e = max(f.items(), key=gccont)
    print(f'{e[0]} is the sequence with the highest GC content = {gccont(e)}')
    return e[0]


# Returns the number of diferences between the two given seqs
def hamming(arq):
    y = open(arq, 'r')
    f_seq = y.readline().rstrip()
    s_seq = y.readline().rstrip()
    cont = 0
    for i in range(len(f_seq)):
        if f_seq[i] != s_seq[i]:
            cont += 1

    return cont


# Translate the mRNA into the aminoacids seq
def translate(arq):
    y = open(arq, 'r').readline()
    print(y)
    cod = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
           'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
           'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
           'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
           'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
           'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
           'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
           'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
           'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
           'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
           'UAA': '', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
           'UAG': '', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
           'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
           'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
           'UGA': '', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
           'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    frame = 0
    encod = []
    while not frame + 1 >= len(y):
        print(frame)
        for codon in cod.keys():
            print(codon)
            if codon == y[frame:frame + 3]:
                encod.append(cod[codon])
                frame += 3
                print(f'{codon}')

    return ''.join(encod)


# Find all the positions of the given motif into the sequence
def findmotif(arq):
    y = open(arq, 'r')
    seq = y.readline().rstrip()
    motif = y.readline().rstrip()
    loc = []
    frame = 0
    while seq:
        if motif in seq[frame:len(motif) + frame]:
            loc.append(frame + 1)
        if frame - 1 == len(seq) - len(motif):
            break
        frame += 1
        print(frame)

    print(*loc, sep=' ')


# rabbits and fibonacci
def rabbits(n, k):
    na = 1
    nb = 1
    q = 0  # quantidade de individuos

    for i in range(n - 2):
        q = na + nb * k
        nb = na
        na = q

    return q


# rabbits and fibonacci but with mortality, not finished yet
def fib(n):
    print(n)
    if n == 0:
        return 0
    elif n == 1:
        return 1
    else:
        return fib(n - 1) + fib(n - 2)


# Find the most likely common ancestor to the given sequences
def consensus_str(arq):
    data = fastaread(arq)
    a, c, g, t = [], [], [], []
    ls = a, c, g, t
    cs_score = []
    cs_seq = []
    ac = 0

    for sequence in data.values():
        for i, res in enumerate(sequence.lower()):
            if not len(a) == len(sequence):  # alonga as listas de acordo com o tamanho das seqs
                for x in ls:
                    x.append(0)
            if 'c' in res:
                c[i] += 1
            elif 't' in res:
                t[i] += 1
            elif 'a' in res:
                a[i] += 1
            elif 'g' in res:
                g[i] += 1

    for i, y in enumerate(a):
        for x in ls:
            if x[i] > ac:
                ac = x[i]
        cs_score.append(ac)
        ac = 0
        for x in ls:  # essa versao da qualquer uma das opcoes com maior score, msm se haver + d 1
            if x[i] == cs_score[i]:
                if x == a:
                    cs_seq.append('A')
                elif x == c:
                    cs_seq.append('C')
                elif x == t:
                    cs_seq.append('T')
                elif x == g:
                    cs_seq.append('G')
                break

    print(*cs_seq, sep='')
    print('A:', *a,
          '\nC:', *c,
          '\nG:', *g,
          '\nT:', *t,
          sep=' ')


def find_shared_motif(arq):
    data = fastaread(arq)
    seqs = [list(sequence) for sequence in data.values()]
    motifs = [[]]
    i = 0
    sequence_0, sequence_1 = seqs[0], seqs[1]  # just to simplify

    for x, y in [(x, y) for x in zip(sequence_0[::], sequence_0[1::]) for y in zip(sequence_1[::], sequence_1[1::])]:
        print(f'Pairs {"".join(x)} and {"".join(y)} being analyzed...')
        if x == y:
            print(f'Pairs {"".join(x)} and {"".join(y)} match!')
            motifs[i].append(x[0]), motifs[i].append(x[1])
            k = sequence_0.index(x[0]) + 2  # NAO ESTA DEVOLVENDO O NUMERO CERTO
            u = sequence_1.index(y[0]) + 2
            print(k, u)

            # Determines if the rest of the sequence is compatible
            print(f'Starting to elongate the motif {x}...')
            for j, m in enumerate(sequence_1[u::]):
                try:
                    # Checks if the nucleotide is equal for both of the sequences
                    print(f'Analyzing the pair {sequence_0[k + j]}, {m}')
                    if m == sequence_0[k + j]:
                        motifs[i].append(m)
                        print(f'The pair {sequence_0[k + j]}, {m} is equal!')

                    # Stop in the first nonequal residue
                    else:
                        print(f'The pair {sequence_0[k + j]}, {m} is not equal.')
                        break

                except IndexError:
                    print('IndexError, end of the string')
                    
             # Takes the motifs and tests it to the others sequences
            for sequence, motif in[(sequence, motif) for sequence in seqs[2::] for motif in motifs]:
                if ''.join(motif) in sequence:
                    print(f'The motif {motif} was finded in the sequence!')
                else:
                    for i, j in enumerate(motif):
                        motif.drop([j])
                        if ''.join(motif) in sequence:
                            break
                        
                
            
        else:
            i += 1
            motifs.append([])

        return motifs
