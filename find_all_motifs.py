import regex
import sys
import fileinput

def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return "".join([seq_dict[base] for base in reversed(seq)])

TOSF='F[DEAVP][FMILV][DEVRL][MLIEFYAR]'
TOSR='[MLIEFYAR][DEVRL][FMILV]DEAVP]F'

for row in fileinput.input():
    REF_seq = str(regex.split('@', row)[1])
    peptide_length = len(REF_seq)
    patternF = regex.compile(r"((" + TOSF + "){s<=0,i<=0,d<=0})", overlapped=True, concurrent=True, flags=regex.IGNORECASE)
    resultF = patternF.findall(string=REF_seq, overlapped=True, concurrent=True)
    patternR = regex.compile(r"((" + TOSR + "){s<=0,i<=0,d<=0})", overlapped=True, concurrent=True, flags=regex.IGNORECASE)
    resultR = patternR.findall(string=REF_seq, overlapped=True, concurrent=True)
    description = str(regex.split('@', row)[0])
    genome = description.split(' ')[2].split(':')[1]
    peptide_id = description.split(' ')[0].replace('>',"")
    chromosome = description.split(' ')[2].replace('chromosome:',"")
    ensembl_gene_id = description.split(' ')[3].replace('gene:',"")
    ensembl_transcript_id = description.split(' ')[4].replace('transcript:',"")
    gene_biotype = description.split(' ')[5].replace('gene_biotype:',"")
    transcript_biotype = description.split(' ')[6].replace('transcript_biotype:',"")
    # try:
    #     gene_symbol = description.split(' ')[7].replace('gene_symbol:',"")
    # except:
    #     gene_symbol = 'not available'
    # # try:
    #     gene_description = description.split('gene_symbol:')[1].replace("description:", "").split('[Source:')[0]
    # except:
    #     gene_description = 'not available'
    try:
        gene_symbol = description.split('gene_symbol:')[1].replace("gene_symbol:", "").split(' ')[0]
    except:
        gene_symbol = 'not available'

    for i in resultF:
        lst1=[]
        for item in regex.finditer(i[1], REF_seq):
            start = item.start()
            end = item.end()
            TOS = i[1]
            lst1.append(genome)
            lst1.append('F')
            lst1.append(peptide_id)
            lst1.append(chromosome)
            lst1.append(ensembl_gene_id)
            lst1.append(ensembl_transcript_id)
            lst1.append(gene_biotype)
            lst1.append(transcript_biotype)
            lst1.append(gene_symbol)
            #lst1.append(gene_description)
            lst1.append(TOS)
            lst1.append(start)
            lst1.append(end)
            lst1.append(peptide_length)
            lst1.append(round(float(start)/float(peptide_length),2))
            lst1.append('#')
            print(lst1)
    for i in resultR:
            lst1 = []
            for item in regex.finditer(i[1], REF_seq):
                start = item.start()
                end = item.end()
                TOS = i[1]
                lst1.append(genome)
                lst1.append('R')
                lst1.append(peptide_id)
                lst1.append(chromosome)
                lst1.append(ensembl_gene_id)
                lst1.append(ensembl_transcript_id)
                lst1.append(gene_biotype)
                lst1.append(transcript_biotype)
                lst1.append(gene_symbol)
                #lst1.append(gene_description)
                lst1.append(TOS)
                lst1.append(start)
                lst1.append(end)
                lst1.append(peptide_length)
                lst1.append(round(float(start)/float(peptide_length),2))
                lst1.append('#')
                print(lst1)
