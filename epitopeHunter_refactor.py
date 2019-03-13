from copy import deepcopy
from os import listdir
from os.path import isfile, join, dirname
import os
import sys
from collections import defaultdict
from functions import *
import math
import os
import argparse
from mhc_i.src.predict_binding import Prediction


def parse_args():
    parser = argparse.ArgumentParser(description='Predict neoepitopes')
    parser.add_argument('--hla', required=True,
                        help='REQUIRED. Input the path to the hla file.')
    # parser.add_argument('--peptide', required=True,
    #                     help = 'REQUIRED. Input the name of the peptide.')
    parser.add_argument('--patientID', required=True,
                        help='REQUIRED. Input the patient ID.')

    args = parser.parse_args()
    return args

def process_hla_mer(hla, patientID, mer_len, peptide_len):
    '''
    This function processes for 1 hla and 1 mer
    hla: hla_a_03_01_01_01
    patientID: A7-A26G
    mer_len: '8'
    peptide_len: '15
    :return:
    '''


    # Make set
    IEDBList = makeSet('IEDB.txt')
    item = hla.upper().split("_")
    hla_IEDBlabel = item[0] + "-" + item[1] + "*" + item[2] + ":" + item[3]

    IEDBStr = "False"

    # new_value.append(hla_IEDBlabel)

    if hla_IEDBlabel in IEDBList:
        IEDBStr = "True"

    if IEDBStr == "False":
        hla_IEDBlabel = getClosestHLA(hla_IEDBlabel, IEDBList)

    prediction = Prediction()

    prediction.commandline_input_w_file('IEDB_recommended', hla_IEDBlabel, str(mer_len),
                                        'peptides/' + patientID + '.' + str(peptide_len) + '.txt',
                                        patientID + '/' + hla + '/output_IEDB.' + str(peptide_len) + '.txt')

    score_dict = defaultdict(list)
    with open(patientID + '/' + hla + '/output_IEDB.' + str(peptide_len) + '.txt', 'r') as f:  # for example: HLA-A_24-02\output_IEDB.15.txt
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

    with open(patientID + '/' + hla + '/' + str(peptide_len) + '_min_score.txt', "w") as outfile:
        with open('peptides/' + patientID + '.' + str(peptide_len) + '.txt') as f:  # open A7-A26G.15.txt
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
                        all_mer_len_peptide = [seq[i:i + mer_len] for i in range(mer_len)]
                        for mer_len_peptide in all_mer_len_peptide:
                            if len(score_dict[mer_len_peptide]) != 0:
                                peptideScore[mer_len_peptide] = min(score_dict[mer_len_peptide])
                        out = min(peptideScore.items(), key=lambda l: l[1])
                        out_to_print = [trans_, out[0], str(out[1])]
                        print >>outfile, '\t'.join(x for x in out_to_print)


def main():
    args = parse_args()
    hlas = []
    with open(args.hla, 'r') as f:
        for line in f:
            hlas.append(line.rstrip('\n').split('\t')[1:])
    hlas = [item for sublist in hlas for item in sublist] #make a flat list
    print hlas

    # Make output folder:
    if not os.path.exists(args.patientID):
        os.makedirs(args.patientID)

    for hla in hlas:
        if not os.path.exists(args.patientID + '/' + hla):
            os.makedirs(args.patientID + '/' + hla)
        process_hla_mer(hla, args.patientID, 9, 17)

if __name__ == '__main__':
    main()



























    # for idx in range(len(mer_len)):
    #
    #     prediction.commandline_input_w_file('IEDB_recommended', hla, mer_len[idx], 'peptides/' + args.patientID + '.' + peptide_len[idx] + '.txt',
    #                                         args.patientID + '/' + hla_dirname + '/output_IEDB.' + peptide_len[idx] + '.txt')
    #
    #     score_dict = defaultdict(list)
    #     with open(args.patientID + '/' + hla_dirname + '/output_IEDB.' + peptide_len[idx] + '.txt', 'r') as f:  # for example: HLA-A_24-02\output_IEDB.15.txt
    #         for line in f:
    #             if not line.startswith('allele'):
    #                 tempLine = line.rstrip('\n').split('\t')
    #                 if "netmhcpan" in tempLine[6]:
    #                     score_dict[tempLine[5]].append(float(tempLine[14]))
    #                 elif len(tempLine) < 9:
    #                     score_dict[tempLine[5]].append(float(tempLine[6]))
    #                 else:
    #                     score_dict[tempLine[5]].append(float(tempLine[8]))
    #
    #     transcript = ""
    #     trans_ = ""
    #     counter = 0
    #     # mutant_Map = defaultdict(list)  # Key is the transcript name and values are all the possible x-mers
    #     # with open(sys.argv[4]) as f:  # open A7-A26G.15.txt
    #
    #     outfile = open(outputFilePath + '/' + hla_dirname + '/' + peptide_len[idx] + '_min_score.txt', "w")
    #     with open('peptides/' + args.patientID + '.' + peptide_len[idx] + '.txt') as f:  # open A7-A26G.15.txt
    #         for line in f:
    #             seq = ""
    #             if ">" in line:
    #                 transcript = line.strip().replace(">", "")
    #                 trans_ = transcript.replace("MT.", "").replace("WT.", "")
    #                 counter += 1
    #             else:
    #                 seq = line.strip()
    #                 counter += 1
    #             if counter % 2 == 0:
    #                 if "MT." in transcript:
    #                     peptideScore = defaultdict()
    #                     all_mer_len_peptide = [seq[i:i + int(mer_len[idx])] for i in range(int(mer_len[idx]))]
    #
    #                     # print all_mer_len_peptide
    #                     # for i in all_mer_len_peptide:
    #                     #     mutant_Map[trans_].append(i)
    #                     for mer_len_peptide in all_mer_len_peptide:
    #                         # print score_dict[mer_len_peptide]
    #                         if len(score_dict[mer_len_peptide]) != 0:
    #                             # print min(score_dict[mer_len_peptide])
    #                             peptideScore[mer_len_peptide] = min(score_dict[mer_len_peptide])
    #                     # print peptideScore
    #                     out = min(peptideScore.items(), key=lambda l: l[1])
    #                     out_to_print = [trans_, out[0], str(out[1])]
    #                     print >>outfile, '\t'.join(x for x in out_to_print)
    #         outfile.close()


########################################################################################
########################################################################################

# prediction.commandline_input_w_file('IEDB_recommended', hla, '9', 'peptides/' + args.patientID + '.17.txt', args.patientID + '/' + hla_dirname + '/output_IEDB.17.txt')
# prediction.commandline_input_w_file('IEDB_recommended', hla, '10', 'peptides/' + args.patientID + '.19.txt',
#                                args.patientID + '/' + hla_dirname + '/output_IEDB.19.txt')
#
# prediction.commandline_input_w_file('IEDB_recommended', hla, '11', 'peptides/' + args.patientID + '.21.txt',
#                                args.patientID + '/' + hla_dirname + '/output_IEDB.21.txt')


# ############
# # After IEDB
# ############
#
# # arguments should be the length of the mer and hla types
# for hla in hlas:
#     for idx in range(len(mer_len)):
#         # mer_len = int(sys.argv[1])  # for example 8
#         # peptide_len = int(sys.argv[2])  # for example 15
#
#         # Obtain a dictionary of score
#
#         score_dict = defaultdict(list)  # key is the peptide and values are scores
#         score_file =
#         with open(sys.argv[3], 'r') as f: #for example: HLA-A_24-02\output_IEDB.15.txt
#              for line in f:
#                   if not line.startswith('allele'):
#                        tempLine = line.rstrip('\n').split('\t')
#                        if "netmhcpan" in tempLine[6]:
#                             score_dict[tempLine[5]].append(float(tempLine[14]))
#                        elif len(tempLine) < 9:
#                             score_dict[tempLine[5]].append(float(tempLine[6]))
#                        else:
#                             score_dict[tempLine[5]].append(float(tempLine[8]))
#
#         transcript = ""
#         trans_ = ""
#         counter = 0
#         mutant_Map = defaultdict(list)  # Key is the transcript name and values are all the possible x-mers
#         with open(sys.argv[4]) as f:  # open A7-A26G.15.txt
#              for line in f:
#                   seq = ""
#                   if ">" in line:
#                        transcript = line.strip().replace(">", "")
#                        trans_ = transcript.replace("MT.", "").replace("WT.", "")
#                        counter += 1
#                   else:
#                        seq = line.strip()
#                        counter += 1
#                   if counter % 2 == 0:
#                        if "MT." in transcript:
#                             peptideScore = defaultdict()
#                             all_mer_len_peptide = [seq[i:i + mer_len] for i in range(mer_len)]
#
#                             # print all_mer_len_peptide
#                             # for i in all_mer_len_peptide:
#                             #     mutant_Map[trans_].append(i)
#                             for mer_len_peptide in all_mer_len_peptide:
#                                  # print score_dict[mer_len_peptide]
#                                  if len(score_dict[mer_len_peptide]) != 0:
#                                       # print min(score_dict[mer_len_peptide])
#                                       peptideScore[mer_len_peptide] = min(score_dict[mer_len_peptide])
#                             # print peptideScore
#                             out = min(peptideScore.items(), key=lambda l: l[1])
#                             out_to_print = [trans_, out[0], str(out[1])]
#                             print '\t'.join(x for x in out_to_print)




# ## Get the current directory
# filepath = os.path.dirname(os.path.realpath(__file__)) + "/"
# hlapath = filepath+ "hla"
# # allfiles = [f for f in listdir(hlapath) if isfile(join(hlapath, f)) and f.endswith(".txt")]
# # allfiles = set(allfiles)
# length = 11
# tab = "\t"
# #print allfiles
#
# file_hlas_map =dict() #This is a dictionary where the key is the filename of hla and the values are the hlas. Since I'm doing this for each patient, is the dictionary still necessary?
# # # for file in allfiles:
# file=sys.argv[2] #TODO: This is the hla filename. Make this into an argument.
# # # file_hlas_map[file] = list();
# # # with open(hlapath+ "/"+ file) as f:
# # #      for line in f:
# # #           hlas = line.strip().split("\t")
# # #           file_hlas_map[file].extend(list(set(hlas[1:len(hlas)])))
# #
# with open(hlapath+ "/"+ file) as f:
#      hla_list = [line.rstrip('\n').split('\t')[1:] for line in f]
#      file_hlas_map[file] = [item for sublist in hla_list for item in sublist]
# #
# # # outputFolder =  str(datetime.utcnow()).replace(" ","-").replace(":","-").replace(".","-") + "/"
# outputFolder = 'CRC-AN' #TODO: Rename output folder to be the name of the patient
# # outputFilePath = filepath + "output/" + outputFolder
# outputFilePath = filepath + outputFolder
# # print filepath
# # print outputFilePath
# isomers = ['15','17','19','21'] #TODO: move this somewhere else.
# # # print file_hlas_map
#
# #############################################
# # Core Processing
# #############################################
#
# procs = []


# for index, patient in enumerate(patientSet):
#      processPatient = True
#      while processPatient and index < len(patientSet)  :
#           # if psutil.cpu_percent() < 85.0:
#           hlas = patientHlaMAP[patient]
#           print patient + "-- " + str(hlas)
#           for hla in hlas:
#                #if not os.path.exists(outputFilePath+ patient +"/" + hla.replace(":","-") + "/"):
#                os.makedirs(outputFilePath+ patient +"/" + hla.replace(":","-") + "/")
#                for num in isomers:
#                     filename = "TCGA-" + patient +"_Varscan_variants_filter.pass."+ str(num) +".peptide"
#                     filename = sys.argv[1] + str(num) + ".peptide"
#                     functions.writeInputFile(filepath + "/peptides/", filename, patient,str(num))
#                     proc = subprocess.Popen([sys.executable, 'epitope.py', hla, num , filename,  patient, outputFolder, syfpeithiStr, IEDBStr, netmhcpanStr])
#                     procs.append(proc)
#           processPatient = False
#           #      time.sleep(5)
#           # else:
#           #      time.sleep(60)
#      print patient + tab + str(patientHlaMAP[patient])

# Re-write to just take in 1 patient
# patient = 'A7-A26G'
# hlas = ['hla_a_03_01_01_01', 'hla_a_30_02_04', 'hla_b_15_16_01', 'hla_b_57_03_01', 'hla_c_16_01_01', 'hla_c_07_01_09'] #TODO: this is placeholder for testing
# prediction = Prediction()
# prediction.commandline_input_w_file('IEDB_recommended', 'HLA-A*30:02', '9', '/scratch/tphung3/Neoepitope_Prediction/peptides/A7-A26G.17.txt', 'output_IEDB.17.txt')

# for hla in hlas:
#      os.makedirs(outputFilePath+ patient +"/" + hla.replace(":","-") + "/")
#      for num in isomers:
#           peptide = sys.argv[3] + str(num) + ".peptide"
#           writeInputFile(filepath + "/peptides/", peptide, patient,str(num)) #TODO: What is this output file?
#           proc = subprocess.Popen([sys.executable, 'epitope.py', hla, num , peptide,  patient, outputFolder, syfpeithiStr, IEDBStr, netmhcpanStr])
#           procs.append(proc)
# print patient + tab + str(patientHlaMAP[patient])


# for proc in procs:
#      proc.wait()
#
# for patient in patientSet:
# patient = 'A7-A26G'
# patient = 'CRC-AN'
# fwrite = open(sys.argv[1],'w') #TODO: make this into an argument
# fwrite.write(getHeaderText())
# # hlas = ['hla_a_03_01_01_01', 'hla_a_30_02_04', 'hla_b_15_16_01', 'hla_b_57_03_01', 'hla_c_16_01_01', 'hla_c_07_01_09']
# # hlas = ['HLA-B*57-01', 'HLA-B*15-25', 'HLA-A*30-02', 'HLA-C*16-01', 'HLA-C*07-01', 'HLA-A*03-01']
# hlas = patientHlaMAP[patient]
# for hla in hlas:
#      print hla
#      IEDB_transcriptMap_MT = defaultdict(list)
#      IEDB_transcriptMap_WT = defaultdict(list)
#      IEDB_TranscriptMap_SameSeq = defaultdict(list)
#      transcript_SET_G = set()
#      dataTuple = namedtuple("dataTuple", ("mer","data"))
#
#      for num in isomers:
#           print num
#
#           transcript_SET, mutant_Map, wildType_Map = getMutantWildTypeData(filepath, patient, num)
#           peptideLenStr = "."+ str(num) + ".txt"
#           print outputFilePath
#           outputFolder = getFilePath(outputFilePath,hla,patient)
#           print outputFolder
#           outputIEDBFile= outputFolder + "output_IEDB" + peptideLenStr
#           IEDBMap, IEDBTuples =  getMapwithValuesIEDB(outputIEDBFile)
#
#           ###################################################
#           # Transcript - peptide - min score
#           ###################################################
#
#           peptideSet = set()
#           transcript_SET_G = transcript_SET
#           print len(transcript_SET)
#           sameSeqTupleList= list()
#           for t in transcript_SET:
#                mutantLowestScoreTuple = getLowestScore(mutant_Map[t][0], IEDBTuples )
#                IEDB_transcriptMap_MT[t].append( dataTuple( mer = num, data = mutantLowestScoreTuple))
#                IEDB_transcriptMap_WT[t].append(dataTuple( mer= num, data = getLowestScore(wildType_Map[t][0], IEDBTuples )))
#                IEDB_TranscriptMap_SameSeq[t].append(dataTuple( mer= num, data = getSameSeqScore(mutantLowestScoreTuple, wildType_Map[t][0], IEDBTuples )))
#                peptideSet.update(getPeptides(IEDB_transcriptMap_MT[t],num))
#                peptideSet.update(getPeptides(IEDB_transcriptMap_WT[t],num) )
#                peptideSet.update(getPeptides(IEDB_TranscriptMap_SameSeq[t],num) )
#
#      for transcript in transcript_SET_G:
#           fwrite.write(patient+  tab + hla+ tab + transcript )
#           merList  = list()
#           for num in isomers:
#                position, peptide, bindingScore = getData(IEDB_transcriptMap_MT[transcript], num)
#                merList.append((num,float(bindingScore), peptide))
#                if num == "15":
#                     fwrite.write(tab  + "MT")
#                fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
#           fileData = getFinalMer(merList)
#           mer = fileData[0]
#           score = fileData[1]
#           peptide = fileData[2]
#           fwrite.write(tab+ peptide + tab + mer + tab+ str(score))
#
#           merList  = list()
#           for num in isomers:
#                position, peptide, bindingScore = getData(IEDB_transcriptMap_WT[transcript], num)
#                merList.append((num,float(bindingScore),peptide))
#                if num == "15":
#                     fwrite.write(tab  + "WT")
#                fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
#           fileData = getFinalMer(merList)
#           mer = fileData[0]
#           score = fileData[1]
#           peptide = fileData[2]
#           fwrite.write(tab+ peptide + tab + mer + tab+ str(score))
#
#           merList  = list()
#           for num in isomers:
#                position, peptide, bindingScore = getData(IEDB_TranscriptMap_SameSeq[transcript], num)
#                merList.append((num,float(bindingScore),peptide))
#                if num == "15":
#                     fwrite.write(tab  + "Same_Seq")
#                fwrite.write(tab + str(position) + tab + str(peptide) + tab + str(bindingScore))
#           fileData = getFinalMer(merList)
#           mer = fileData[0]
#           score = fileData[1]
#           peptide = fileData[2]
#           fwrite.write(tab+ peptide + tab + mer + tab+ str(score))
#
#           newHla = getNewHLA(newOldHLAMap, hla)
#           fwrite.write(tab + newHla+"\n")
#
#
# fwrite.close()


#############################################
# End of Program
#############################################
