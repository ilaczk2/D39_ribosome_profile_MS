"""
-----
The project is for identifying sORF, small open-reading frame based on rib-seq data
-----
"""
import numpy as np

genecode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

genecode_antisense = {
    'TAT': 'I', 'TAG': 'I', 'TAA': 'I', 'TAC': 'M',
    'TGT': 'T', 'TGG': 'T', 'TGC': 'T', 'TGA': 'T',
    'TTG': 'N', 'TTA': 'N', 'TTT': 'K', 'TTC': 'K',
    'TCG': 'S', 'TCA': 'S', 'TCT': 'R', 'TCC': 'R',
    'GAT': 'L', 'GAG': 'L', 'GAC': 'L', 'GAA': 'L',
    'GGT': 'P', 'GGG': 'P', 'GGC': 'P', 'GGA': 'P',
    'GTG': 'H', 'GTA': 'H', 'GTT': 'Q', 'GTC': 'Q',
    'GCT': 'R', 'GCG': 'R', 'GCC': 'R', 'GCA': 'R',
    'CAT': 'V', 'CAG': 'V', 'CAC': 'V', 'CAA': 'V',
    'CGT': 'A', 'CGG': 'A', 'CGC': 'A', 'CGA': 'A',
    'CTG': 'D', 'CTA': 'D', 'CTT': 'E', 'CTC': 'E',
    'CCT': 'G', 'CCG': 'G', 'CCC': 'G', 'CCA': 'G',
    'AGT': 'S', 'AGG': 'S', 'AGC': 'S', 'AGA': 'S',
    'AAG': 'F', 'AAA': 'F', 'AAT': 'L', 'AAC': 'L',
    'ATG': 'Y', 'ATA': 'Y', 'ATT': '_', 'ATC': '_',
    'ACG': 'C', 'ACA': 'C', 'ACT': '_', 'ACC': 'W'
}

start_codon_list = ['ATG', 'GTG', 'CTG', 'TTG', 'ATT', 'ATC']
stop_codon_list = ['TAA', 'TAG', 'TGA']

start_codon_antisense_list = ['CAT', 'CAC', 'CAG', 'CAA', 'AAT', 'GAT']
stop_codon_antisense_list = ['ATT', 'ATC', 'ACT']


def genome_extract(gb_file: str):
    """
    return the genome seq at the bottom of the gb file as string.
    -----
    :parameter: a gb file path
    """
    import re

    with open(gb_file, 'r') as f:
        file_split = f.read().split('ORIGIN')[1].split('//')[0]
        genome_seq = re.sub('\d+', '', file_split)
        genome_seq = re.sub('\s+', '', genome_seq)

        print('genome has %i nucleic acids' % len(genome_seq))
    return genome_seq


def wig_file_reader(wig_file_path: str):
    """
    read position and corresponding RPM into a 2D list, [(pos1, RPM1),(pos2, RPM2)]
    returns a 2D list of tuples
    -----
    :param wig_file_path: a wig file path
    """
    pos_rpm_list = []
    with open(wig_file_path, 'r') as f:
        # skip first two lines
        for i in range(2):
            next(f)
        for line in f:
            line_split = line.split('\t')
            pos = int(line_split[0]) - 1
            rpm = float(line_split[1])
            pos_rpm_list.append((pos, rpm))
    return pos_rpm_list


def wig_info_into_array(pos_rpm_list: list, zero_array):
    """
    get all position and rpm data into a numpy zero array, length as genome
    returns a numpy array
    :param pos_rpm_list: list of tuples, [(pos1,rpm1),(pos2,rpm2)...]
    :param zero_array:
    """
    for each in pos_rpm_list:
        try:
            zero_array[each[0]] = each[1]  # some positions in wig file are out of genome seq range
        except IndexError:
            continue
    return zero_array


def annotation_reader(annotation_path: str):
    """
    reads start and end position for each annotated
    returns a 2D list of position tuples
    -----
    :param annotation_path: an annotation file path
    """
    start_end_pos_list_sense = []
    start_end_pos_list_anti_sense = []
    with open(annotation_path, 'r') as f:
        for i in range(1):
            next(f)
        for line in f:
            line_split = line.split('\t')
            start_pos = int(line_split[1])
            end_pos = int(line_split[2])
            strand = line_split[3]
            if strand == '+':
                start_end_pos_list_sense.append((start_pos, end_pos))
            elif strand == '-':
                start_end_pos_list_anti_sense.append((start_pos, end_pos))
    return start_end_pos_list_sense, start_end_pos_list_anti_sense


def intergenetic_rpm_array_getter(num_array, annotation_start_end_pos_list: list):
    """
    zero the annotated gene on num_array, the non_zero parts are intergenetic
    returns a new numpy 1D array
    -----
    :param num_array: 1D numpy array
    :param annotation_start_end_pos_list: list of tuples with start and end position of annotated gene in each tuple
    """
    for each in annotation_start_end_pos_list:
        num_array[each[0]:each[1]] = 0

    return num_array


def main_function(intergenetic_rpm_array, genome_seq, annotated_range_array, rpm_threshold=10):
    """
    sORF finding function for sense strand
    returns two dictionaries, one records start pos and the protein seq that sORF translate into,
    one records start pos and nucleotides for sORF
    -----
    :param intergenetic_rpm_array: numpy array after zero all the annotated genee
    :param genome_seq: string, genome sequence
    :param annotated_range_array: a list of annotated gene position range
    """
    # make sequence in genome sequence uppercase
    genome_seq = genome_seq.upper()
    # find the non zeros rpm from intergenentic rpm array, index them
    non_zero_indices = [j for i in np.nonzero(intergenetic_rpm_array) for j in i]
    unannotated_protein_seq_dictionary = {}  # save result
    unannotated_nt_dictionary = {}
    rpm_statistic = []  # get average rpm for identified sORF
    # iterate each non zero rpm and check if it's higher than threshold
    for each_index in non_zero_indices:

        if intergenetic_rpm_array[each_index] > rpm_threshold:
            checking_start, checking_end = each_index - 10, each_index + 10
            checking_region = genome_seq[each_index - 10: each_index + 10]  # check 20nt regions
            for each_start in start_codon_list:

                # check if any start codon within 20nt regions
                if each_start in checking_region:
                    print(each_index)
                    start_pos = checking_region.find(each_start) + checking_start
                    seq = 'M'  # initialize translation, M as methionine
                    sorf_nt = ''
                    i = 3  # skip start codon position and read from second amino acid
                    j = 6
                    stop = False  # stop codon switch
                    while not stop:  # when nt is not in stop codons list
                        # 3 nucleotides to be translated
                        nt = genome_seq[start_pos + i:start_pos + j]
                        sorf_nt += nt
                        # if read a stop codon, stop translating
                        if nt in stop_codon_list:
                            stop = True
                        else:
                            aa = genecode[nt]
                            seq += aa
                            i += 3
                            j += 3

                    # also check if the position of stop codon falls in the range of #annotated gene
                    # if the position of last nucleotide falls into any annotated range,
                    # it’s not an intergenetic translation
                    # 2<seq<50 are selected
                    rpm_mean = np.mean(intergenetic_rpm_array[start_pos:start_pos + j])
                    if 7 < len(seq) < 50 and start_pos + j - 1 not in annotated_range_array and rpm_mean > 1:
                        unannotated_protein_seq_dictionary[
                            str(start_pos + 1) + '_' + str(start_pos + j) + '_sense'] = seq
                        unannotated_nt_dictionary[str(start_pos + 1) + '_' + str(start_pos + j) + '_sense'] = sorf_nt
                        rpm_statistic.append(rpm_mean)
    rpm_statistic = np.array(rpm_statistic)
    print('mean: %s, deviation: %s' % (np.mean(rpm_statistic), np.std(rpm_statistic)))
    # check if the position of stop codon falls in the range of any annotated gene
    return unannotated_protein_seq_dictionary, unannotated_nt_dictionary


def main_function_antisense(intergenetic_rpm_array, genome_seq, annotated_range_array, rpm_threshold=10):
    """
    sORF finding function for anti-sense strand
    returns two dictionaries, one records start pos and the protein seq that sORF translate into,
    one records start pos and nucleotides for sORF
    -----
    :param intergenetic_rpm_array:
    :param genome_seq:
    :param annotated_range_array:
    :param rpm_threshold:

    """
    # make sequence in genome sequence uppercase
    genome_seq = genome_seq.upper()
    # find the non zeros rpm from intergenentic rpm array, index them
    non_zero_indices = [j for i in np.nonzero(intergenetic_rpm_array) for j in i]
    unannotated_protein_seq_dictionary = {}  # save result
    unannotated_nt_dictionary = {}
    rpm_statistic = []  # get average rpm for identified sORF
    # iterate each non zero rpm and check if it's higher than threshold
    for each_index in non_zero_indices:

        if intergenetic_rpm_array[each_index] > rpm_threshold:
            checking_start, checking_end = each_index - 10, each_index + 10
            checking_region = genome_seq[each_index - 10: each_index + 10]  # check 20nt regions, anti-strand
            for each_start in start_codon_antisense_list:

                # check if any start codon within 20nt regions
                if each_start in checking_region:
                    print(each_index)
                    start_pos = checking_region.find(each_start) + checking_start
                    seq = 'M'  # initialize translation, M as methionine
                    sorf_nt = ''
                    i = 0  # skip start codon position and read from second amino acid
                    j = 3
                    stop = False  # stop codon switch
                    while not stop:  # when nt is not in stop codons list
                        # 3 nucleotides to be translated
                        nt = genome_seq[start_pos - j:start_pos - i][::-1]
                        sorf_nt += nt
                        # if read a stop codon, stop translating
                        if nt in stop_codon_antisense_list:
                            stop = True
                        else:
                            aa = genecode_antisense[nt]
                            seq += aa
                            i += 3
                            j += 3

                    # also check if the position of stop codon falls in the range of #annotated gene
                    # if the position of last nucleotide falls into any annotated range,
                    # it’s not an intergenetic translation
                    # 7<seq<50, rpm mean >1 are selected
                    rpm_mean = np.mean(intergenetic_rpm_array[start_pos - j:start_pos + 3])
                    if 7 < len(seq) < 50 and start_pos - j not in annotated_range_array and rpm_mean > 1:
                        unannotated_protein_seq_dictionary[
                            str(start_pos - j + 1) + '_' + str(start_pos + 3) + '_anti_sense'] = seq
                        unannotated_nt_dictionary[
                            str(start_pos - j + 1) + '_' + str(start_pos + 3) + '_anti_sense'] = sorf_nt
                        rpm_statistic.append(rpm_mean)
    rpm_statistic = np.array(rpm_statistic)
    print('mean: %s, deviation: %s' % (np.mean(rpm_statistic), np.std(rpm_statistic)))
    # check if the position of stop codon falls in the range of any annotated gene
    return unannotated_protein_seq_dictionary, unannotated_nt_dictionary


def annotation_gene_rpm_statistics(annotation_pos_list: list, rpm_array):
    """
    calculate the average rpm for each annotated gene
    returns a 1D array with average rpm for each gene
    -----
    :param annotation_pos_list: [(start1,end1),(start2,end2)...]
    :param rpm_array: numpy array after dumping rpm info from wig file
    """

    average_array = []
    for each in annotation_pos_list:
        annotated_rpm_array = rpm_array[each[0] - 1:each[1] - 1]
        average = np.mean(annotated_rpm_array)
        average_array.append(average)
    average_array = np.array(average_array)
    print('average rpm in annotation genes: %s, deviation: %s' % (np.mean(average_array), np.std(average_array)))
    return average_array


def find_overlapping_sORFs(file_list_in, file_out):
    """
    find all overlapping sORFs from all files, sense and anti-sense strand
    -----
    :param file_list_in:[(pos1,neg1), (pos2,neg2), (pos3,neg3)]
    :return:
    """
    from protein_coverage import fasta_reader
    pos_dict, neg_dict = fasta_reader(file_list_in[0][0]), fasta_reader(file_list_in[0][1])
    pos_dict.update(neg_dict)

    combined_rpm_list = []
    with open(file_out, 'w') as f_out:
        for each in file_list_in:
            with open(each[0], 'r') as f_pos:
                pos_list = [line.split('|')[1] for line in f_pos if line.startswith('>')]
            with open(each[1], 'r') as f_neg:
                neg_list = [line.split('|')[1] for line in f_neg if line.startswith('>')]
            combine_list = pos_list + neg_list
            combined_rpm_list.append(combine_list)
        common_rpm = list(set(combined_rpm_list[0]) & set(combined_rpm_list[1]) & set(combined_rpm_list[2]))
        for each in common_rpm:
            f_out.write('>sp|' + each + '|\n')
            seq = pos_dict[each]
            bins = range(0, len(seq) + 60, 60)
            for i in range(len(bins) - 1):
                f_out.write(seq[bins[i]:bins[i + 1]] + '\n')


def output_result(gb_path: str, wig_path: str, annotation_path: str, output_path: str, strand='sense'):
    """
    output all identified sORFs and their genome coordinates to a fasta file
    -----
    :param gb_path:
    :param wig_path:
    :param annotation_path:
    :param strand: default as sense, could use anti-sense
    :return:
    """
    genome_seq = genome_extract(gb_path)  # extract genome seq from gb file
    pos_rpm_list = wig_file_reader(wig_path)  # extract rpm/position infomation from wig file
    zero_array = np.zeros(len(genome_seq))  # initiate a zero-array, length same as genome seq
    rpm_array = wig_info_into_array(pos_rpm_list, zero_array)  # get rpm into zero array

    # read both sense and anti-sense annotation
    anno_pos_list_sense, anno_pos_list_antisense = annotation_reader(annotation_path)

    # get annotated position range
    annotated_range_array = [i for each in anno_pos_list_sense for i in np.arange(each[0], each[1])]
    annotated_range_antisense_array = [i for each in anno_pos_list_antisense for i in np.arange(each[0], each[1])]

    if strand == 'sense':  # for sense strand
        new_rpm_array = intergenetic_rpm_array_getter(rpm_array, anno_pos_list_sense)
        protein_seq_dict = main_function(new_rpm_array, genome_seq, annotated_range_array)[0]
    elif strand in ['antisense', 'anti-sense', 'reverse']:  # for antisense strand
        new_rpm_antisense_array = intergenetic_rpm_array_getter(rpm_array, anno_pos_list_antisense)
        protein_seq_dict = \
        main_function_antisense(new_rpm_antisense_array, genome_seq, annotated_range_antisense_array)[0]
    else:
        raise ValueError('strand has to be sense or anti-sense')

    # output result to a fasta file
    with open(output_path, 'w') as f_open:
        for each in protein_seq_dict:
            seq = protein_seq_dict[each]
            f_open.write(">sp|" + each + '|' + '\n')
            bins = range(0, len(seq) + 60, 60)
            for i in range(len(bins) - 1):
                f_open.write(seq[bins[i]:bins[i + 1]] + '\n')
    return protein_seq_dict

if __name__ == '__main__':
    import numpy as np

    gb_path = '/home/xshao/lab/Irina/sequence (10).gb'
    wig_path = '/home/xshao/lab/Irina/Wig files/dataset_24550.datRET_Spn.negRPM.wig'
    annotation_path = '/home/xshao/lab/Irina/Mochi_view_CP027540_Genes_annotation.txt'

    output_result(gb_path,wig_path,annotation_path,'irina_test.fasta',strand='antisense')

