from collections import defaultdict
import sys
# arguments should be the length of the mer and hla types
mer_len = int(sys.argv[1]) #for example 8
peptide_len = int(sys.argv[2]) #for example 15

# Obtain a dictionary of score


score_dict = defaultdict(list) #key is the peptide and values are scores
with open(sys.argv[3], 'r') as f:
    for line in f:
        if not line.startswith('allele'):
            tempLine = line.rstrip('\n').split('\t')
            if "netmhcpan" in tempLine[6]:
                score_dict[tempLine[5]].append(float(tempLine[14]))
            elif len(tempLine) < 9:
                score_dict[tempLine[5]].append(float(tempLine[6]))
            else:
                score_dict[tempLine[5]].append(float(tempLine[8]))


transcript = ""
trans_ = ""
counter = 0
mutant_Map = defaultdict(list)  # Key is the transcript name and values are all the possible x-mers
with open(sys.argv[4]) as f:  # open A7-A26G.15.txt
    for line in f:
        seq = ""
        if ">" in line:
            transcript = line.strip().replace(">", "")
            trans_ = transcript.replace("MT.", "").replace("WT.", "")
            counter += 1
        else:
            seq = line.strip()
            counter += 1
        if counter % 2 == 0:
            if "MT." in transcript:
                peptideScore = defaultdict()
                all_mer_len_peptide = [seq[i:i+mer_len] for i in range(mer_len)]

                # print all_mer_len_peptide
                # for i in all_mer_len_peptide:
                #     mutant_Map[trans_].append(i)
                for mer_len_peptide in all_mer_len_peptide:
                    # print score_dict[mer_len_peptide]
                    if len(score_dict[mer_len_peptide]) != 0:
                        # print min(score_dict[mer_len_peptide])
                        peptideScore[mer_len_peptide] = min(score_dict[mer_len_peptide])
                # print peptideScore
                out = min(peptideScore.items(), key=lambda l: l[1])
                out_to_print = [out[0], str(out[1])]
                print '\t'.join(x for x in out_to_print)





# print mutant_Map
# peptide = 'TFRHSVVVPHEPPEVGSDC'
# length = 10
# seqs = [peptide[i:i+length] for i in range(length)]
# # print seqs
# dict = defaultdict(list)
# with open("c://Users/tuyen/Documents/github_repo/Neoepitope_Prediction/A7-A26G/HLA-A_03-01/output_IEDB.19.txt", 'r') as f:
#     for line in f:
#         if not line.startswith('allele'):
#             items = line.rstrip('\n').split('\t')
#             dict[items[5]].append(float(items[8]))
# print dict
# peptideScore = defaultdict()
# for seq in seqs:
#     peptideScore[seq] = min(dict[seq])
# print peptideScore


# print min(peptideScore.items(), key=lambda l: l[1])

# all_scores_flat = [item for sublist in all_scores for item in sublist]
# print all_scores_flat
# print min(all_scores_flat)
