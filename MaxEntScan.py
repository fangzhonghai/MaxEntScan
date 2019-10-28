# -*- coding:utf-8 -*-
import pyfaidx
# import vcf
import pandas as pd
import optparse
import sys
import os
import numpy as np
import yaml
import collections
pd.options.mode.chained_assignment = None


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    MaxEntScan::score5ss for human 5' splice sites AND MaxEntScan::score3ss for human 3' splice sites
    This program Packaged the 'MaxEntScan'.
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print usage_message
    sys.exit()


def read_hgmd_vcf(hgmd_file_line):
    line_list = hgmd_file_line.strip().split('\t')
    line_list[1] = int(line_list[1])
    return line_list


def base_complement(dna_base):
    complement_rule = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_rule[dna_base]


def read_genome_fa(genome_fa):
    read_fa = pyfaidx.Fasta(genome_fa)
    return read_fa


def locate_mut(vcf_record, exon_gff):
    if len(exon_gff.loc[(exon_gff['Chrom'] == vcf_record[0]) & (exon_gff['Exon_Start'] <= vcf_record[1]) & (vcf_record[1] <= exon_gff['Exon_End']), ].values):
        mut_loc = 'In_Exon'
    else:
        mut_loc = 'In_Intron'
    return mut_loc


def intron_mut_nearby(vcf_record, exon_gff, offset):
    nearby_judge = len(exon_gff.loc[(exon_gff['Chrom'] == vcf_record[0]) & (exon_gff['Exon_Start'] <= vcf_record[1]+offset) & (vcf_record[1]+offset <= exon_gff['Exon_End']), ].values)
    return nearby_judge


def exon_range(vcf_record, exon_gff, offset):
    exon_start = exon_gff.loc[(exon_gff['Chrom'] == vcf_record[0]) & (exon_gff['Exon_Start'] <= vcf_record[1]+offset) & (vcf_record[1]+offset <= exon_gff['Exon_End']), 'Exon_Start'].values[0]
    exon_end = exon_gff.loc[(exon_gff['Chrom'] == vcf_record[0]) & (exon_gff['Exon_Start'] <= vcf_record[1]+offset) & (vcf_record[1]+offset <= exon_gff['Exon_End']), 'Exon_End'].values[0]
    return exon_start, exon_end


def strand_orientation(vcf_record, exon_gff, exon_start, exon_end):
    strand = exon_gff.loc[(exon_gff['Chrom'] == vcf_record[0]) & (exon_gff['Exon_Start'] == exon_start) & (exon_gff['Exon_End'] == exon_end), 'Strand'].values[0]
    return strand


def pos_strand_exon_mut_ss5_seq(vcf_record, read_fa, exon_end):
    ss5_seq = str(read_fa.get_seq(vcf_record[0], exon_end-2, exon_end+6))
    ss5_seq = ss5_seq[0:3].lower() + ss5_seq[3:5] + ss5_seq[5:].lower()
    ss5_mut_seq = ss5_seq[0:vcf_record[1] - exon_end+2] + vcf_record[3] + ss5_seq[vcf_record[1] - exon_end+3:]
    return ss5_seq, ss5_mut_seq


def neg_strand_exon_mut_ss5_seq(vcf_record, read_fa, exon_start):
    ss5_seq = str(read_fa.get_seq(vcf_record[0], exon_start-6, exon_start+2).reverse.complement)
    ss5_seq = ss5_seq[0:3].lower() + ss5_seq[3:5] + ss5_seq[5:].lower()
    ss5_mut_seq = ss5_seq[0:exon_start - vcf_record[1]+2] + base_complement(vcf_record[3]) + ss5_seq[exon_start - vcf_record[1]+3:]
    return ss5_seq, ss5_mut_seq


def pos_strand_exon_mut_ss3_seq(vcf_record, read_fa, exon_start):
    ss3_seq = str(read_fa.get_seq(vcf_record[0], exon_start - 20, exon_start + 2))
    ss3_seq = ss3_seq[0:18].lower() + ss3_seq[18:20] + ss3_seq[20:].lower()
    ss3_mut_seq = ss3_seq[0:20] + ss3_seq[20:vcf_record[1] - exon_start+20] + vcf_record[3] + ss3_seq[vcf_record[1] - exon_start+21:]
    return ss3_seq, ss3_mut_seq


def neg_strand_exon_mut_ss3_seq(vcf_record, read_fa, exon_end):
    ss3_seq = str(read_fa.get_seq(vcf_record[0], exon_end - 2, exon_end + 20).reverse.complement)
    ss3_seq = ss3_seq[0:18].lower() + ss3_seq[18:20] + ss3_seq[20:].lower()
    ss3_mut_seq = ss3_seq[0:20] + ss3_seq[20:exon_end - vcf_record[1]+20] + base_complement(vcf_record[3]) + ss3_seq[exon_end - vcf_record[1]+21:]
    return ss3_seq, ss3_mut_seq


def pos_strand_intron_mut_ss5_seq(vcf_record, read_fa, exon_end):
    ss5_seq = str(read_fa.get_seq(vcf_record[0], exon_end - 2, exon_end + 6))
    ss5_seq = ss5_seq[0:3].lower() + ss5_seq[3:5] + ss5_seq[5:].lower()
    ss5_mut_seq = ss5_seq[0:3] + ss5_seq[3:vcf_record[1] - exon_end+2] + vcf_record[3] + ss5_seq[vcf_record[1] - exon_end+3:]
    return ss5_seq, ss5_mut_seq


def neg_strand_intron_mut_ss5_seq(vcf_record, read_fa, exon_start):
    ss5_seq = str(read_fa.get_seq(vcf_record[0], exon_start - 6, exon_start + 2).reverse.complement)
    ss5_seq = ss5_seq[0:3].lower() + ss5_seq[3:5] + ss5_seq[5:].lower()
    ss5_mut_seq = ss5_seq[0:3] + ss5_seq[3:exon_start - vcf_record[1]+2] + base_complement(vcf_record[3]) + ss5_seq[exon_start - vcf_record[1]+3:]
    return ss5_seq, ss5_mut_seq


def pos_strand_intron_mut_ss3_seq(vcf_record, read_fa, exon_start):
    ss3_seq = str(read_fa.get_seq(vcf_record[0], exon_start - 20, exon_start + 2))
    ss3_seq = ss3_seq[0:18].lower() + ss3_seq[18:20] + ss3_seq[20:].lower()
    ss3_mut_seq = ss3_seq[0:vcf_record[1] - exon_start+20] + vcf_record[3] + ss3_seq[vcf_record[1] - exon_start+21:]
    return ss3_seq, ss3_mut_seq


def neg_strand_intron_mut_ss3_seq(vcf_record, read_fa, exon_end):
    ss3_seq = str(read_fa.get_seq(vcf_record[0], exon_end - 2, exon_end + 20).reverse.complement)
    ss3_seq = ss3_seq[0:18].lower() + ss3_seq[18:20] + ss3_seq[20:].lower()
    ss3_mut_seq = ss3_seq[0:exon_end - vcf_record[1]+20] + base_complement(vcf_record[3]) + ss3_seq[exon_end - vcf_record[1]+21:]
    return ss3_seq, ss3_mut_seq


def store_info(vcf_record, mut_location, strand, ss_seq, ss_mut_seq, exon_start, exon_end, chrom, pos, mut_place, near_exon_start, near_exon_end, exon_strand, ref, alt, complement_ref, complement_alt, pre_mut_seq, post_mut_seq, splice_site_type):
    chrom.append(vcf_record[0])
    pos.append(vcf_record[1])
    mut_place.append(mut_location)
    near_exon_start.append(exon_start)
    near_exon_end.append(exon_end)
    exon_strand.append(strand)
    ref.append(vcf_record[2])
    alt.append(vcf_record[3])
    complement_ref.append(base_complement(vcf_record[2]))
    complement_alt.append(base_complement(vcf_record[3]))
    pre_mut_seq.append(ss_seq)
    post_mut_seq.append(ss_mut_seq)
    if len(ss_mut_seq) == 9:
        splice_site_type.append('SS5')
    else:
        splice_site_type.append('SS3')
    return chrom, pos, mut_place, near_exon_start, near_exon_end, exon_strand, ref, alt, complement_ref, complement_alt, pre_mut_seq, post_mut_seq, splice_site_type


def write_info(path, sample, chrom, pos, mut_place, near_exon_start, near_exon_end, exon_strand, ref, alt, complement_ref, complement_alt, pre_mut_seq, post_mut_seq, splice_site_type):
    mut_dic = collections.OrderedDict({'#CHROM': chrom, 'POS': pos, 'Mut_Place': mut_place, 'Near_Exon_Start': near_exon_start,
                                       'Near_Exon_End': near_exon_end, 'Exon_Strand': exon_strand, 'REF': ref, 'ALT': alt,
                                       'Complement_REF': complement_ref, 'Complement_ALT': complement_alt, 'Pre_Mut_Seq': pre_mut_seq,
                                       'Post_Mut_Seq': post_mut_seq, 'Splice_Site_Type': splice_site_type})
    mut_df = pd.DataFrame(mut_dic, columns=['#CHROM', 'POS', 'Mut_Place', 'Near_Exon_Start', 'Near_Exon_End', 'Exon_Strand',
                                            'REF', 'ALT', 'Complement_REF', 'Complement_ALT', 'Pre_Mut_Seq', 'Post_Mut_Seq', 'Splice_Site_Type'])
    mut_df.to_csv(path + '/' + sample+'.mut_info.txt', index=False, sep='\t')
    # fp = open(sample+'.mut_info.txt', 'w')
    # fp.write('CHROM'+'\t'+'POS'+'\t'+'Mut_Place'+'\t'+'Near_Exon_Start'+'\t'+'Near_Exon_End'+'\t'+'Exon_Strand'+'\t'+'REF'+'\t'+'ALT'+'\t' +
    #          'Complement_REF'+'\t'+'Complement_ALT'+'\t'+'Pre_Mut_Seq'+'\t'+'Post_Mut_Seq'+'\t'+'Splice_Site_Type'+'\n')
    # for i in range(len(chrom)):
    #     fp.write(chrom[i] + '\t')
    #     fp.write(str(pos[i]) + '\t')
    #     fp.write(mut_place[i] + '\t')
    #     fp.write(str(near_exon_start[i]) + '\t')
    #     fp.write(str(near_exon_end[i]) + '\t')
    #     fp.write(exon_strand[i] + '\t')
    #     fp.write(ref[i] + '\t')
    #     fp.write(alt[i] + '\t')
    #     fp.write(complement_ref[i] + '\t')
    #     fp.write(complement_alt[i] + '\t')
    #     fp.write(pre_mut_seq[i] + '\t')
    #     fp.write(post_mut_seq[i] + '\t')
    #     fp.write(splice_site_type[i] + '\n')
    # fp.close()


def read_mut_info(mut_file):
    all_mut_info = pd.read_table(mut_file)
    uni_mut_info = all_mut_info.drop_duplicates(subset=['#CHROM', 'Near_Exon_Start', 'Near_Exon_End', 'Splice_Site_Type'], keep='first')
    un_rep_info = all_mut_info.drop_duplicates(subset=['#CHROM', 'Near_Exon_Start', 'Near_Exon_End', 'Splice_Site_Type'], keep=False)
    rep_mut_info = uni_mut_info.append(un_rep_info).drop_duplicates(subset=['#CHROM', 'Near_Exon_Start', 'Near_Exon_End', 'Splice_Site_Type'], keep=False)
    return all_mut_info, uni_mut_info, rep_mut_info


def combine_mut_info(all_mut_info, uni_mut_info, rep_mut_info, rep_mut_index, combine_col):
    get_rep_val = all_mut_info.loc[(all_mut_info['#CHROM'] == rep_mut_info.iloc[rep_mut_index, 0]) & (all_mut_info['Near_Exon_Start'] == rep_mut_info.iloc[rep_mut_index, 3]) & (all_mut_info['Near_Exon_End'] == rep_mut_info.iloc[rep_mut_index, 4]) & (all_mut_info['Splice_Site_Type'] == rep_mut_info.iloc[rep_mut_index, 12])][combine_col].values
    combine_val = ';'.join([str(ele) for ele in get_rep_val])
    uni_mut_info.loc[(uni_mut_info['#CHROM'] == rep_mut_info.iloc[rep_mut_index, 0]) & (uni_mut_info['Near_Exon_Start'] == rep_mut_info.iloc[rep_mut_index, 3]) & (uni_mut_info['Near_Exon_End'] == rep_mut_info.iloc[rep_mut_index, 4]) & (uni_mut_info['Splice_Site_Type'] == rep_mut_info.iloc[rep_mut_index, 12]), combine_col] = combine_val


def sub_base(string, sub_pos_list, sub_mut_list):
    new_string = list(string)
    for index, item in enumerate(sub_pos_list):
        new_string[item] = sub_mut_list[index]
    return ''.join(new_string)


def combine_post_mut_seq(all_mut_info, uni_mut_info, rep_mut_info, rep_mut_index):
    combine_pos = []
    combine_base = []
    get_rep = all_mut_info.loc[(all_mut_info['#CHROM'] == rep_mut_info.iloc[rep_mut_index, 0]) & (all_mut_info['Near_Exon_Start'] == rep_mut_info.iloc[rep_mut_index, 3]) & (all_mut_info['Near_Exon_End'] == rep_mut_info.iloc[rep_mut_index, 4]) & (all_mut_info['Splice_Site_Type'] == rep_mut_info.iloc[rep_mut_index, 12])]
    pre_seq = get_rep['Pre_Mut_Seq'].values[0]
    post_seq = get_rep['Post_Mut_Seq'].values
    for j in range(get_rep.shape[0]):
        where_mut = np.array(list(pre_seq)) == np.array(list(post_seq[j]))
        combine_pos.append(np.argwhere(where_mut == False)[0][0])
        combine_base.append(post_seq[j][np.argwhere(where_mut == False)[0][0]])
    real_post_mut_seq = sub_base(pre_seq, combine_pos, combine_base)
    uni_mut_info.loc[(uni_mut_info['#CHROM'] == rep_mut_info.iloc[rep_mut_index, 0]) & (uni_mut_info['Near_Exon_Start'] == rep_mut_info.iloc[rep_mut_index, 3]) & (uni_mut_info['Near_Exon_End'] == rep_mut_info.iloc[rep_mut_index, 4]) & (uni_mut_info['Splice_Site_Type'] == rep_mut_info.iloc[rep_mut_index, 12]), 'Post_Mut_Seq'] = real_post_mut_seq


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-v', '--vcf', dest='vcf_file', default='none', type='string')
    parser.add_option('-s', '--sample', dest='sample', default='none', type='string')
    parser.add_option('-c', '--config', dest='config', default='none', type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default='none', type='string')
    (opts, args) = parser.parse_args()
    vcf_file = opts.vcf_file
    sample = opts.sample
    config = opts.config
    pwd = opts.pwd
    config_yaml = yaml.load(open(config), Loader=yaml.FullLoader)
    maxentscan_path = config_yaml['MaxEntScan_path']
    Read_Genome_Fa = read_genome_fa(config_yaml['hg19'])
    All_Exon = pd.read_csv(config_yaml['GRCh37_gff'], sep='\t')
    All_Exon.columns = ['Chrom', 'Exon_Start', 'Exon_End', 'Strand']
    Uni_Exon = All_Exon.drop_duplicates(subset=['Chrom', 'Exon_Start', 'Exon_End'], keep='first')
    # Vcf_Reader = vcf.Reader(filename=vcf_file)
    CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type = \
        [], [], [], [], [], [], [], [], [], [], [], [], []
    for line in open(vcf_file, 'r'):
        if line[0] == '#':
            continue
        record = read_hgmd_vcf(line)
        if len(record[2]) == 1 and record[2] != '.' and len(record[3]) == 1 and record[3] != '.':
            Mut_Location = locate_mut(record, Uni_Exon)
            if Mut_Location == 'In_Exon':
                Exon_Start, Exon_End = exon_range(record, Uni_Exon, 0)
                Strand = strand_orientation(record, Uni_Exon, Exon_Start, Exon_End)
                if Strand == '+':
                    if record[1] + 2 >= Exon_End:
                        SS5_Seq, SS5_Mut_Seq = pos_strand_exon_mut_ss5_seq(record, Read_Genome_Fa, Exon_End)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS5_Seq, SS5_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                    elif record[1] - 2 <= Exon_Start:
                        SS3_Seq, SS3_Mut_Seq = pos_strand_exon_mut_ss3_seq(record, Read_Genome_Fa, Exon_Start)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS3_Seq, SS3_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                    else:
                        continue
                else:
                    if record[1] - 2 <= Exon_Start:
                        SS5_Seq, SS5_Mut_Seq = neg_strand_exon_mut_ss5_seq(record, Read_Genome_Fa, Exon_Start)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS5_Seq, SS5_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                    elif record[1] + 2 >= Exon_End:
                        SS3_Seq, SS3_Mut_Seq = neg_strand_exon_mut_ss3_seq(record, Read_Genome_Fa, Exon_End)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS3_Seq, SS3_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                    else:
                        continue
            else:
                if intron_mut_nearby(record, Uni_Exon, -6):
                    Exon_Start, Exon_End = exon_range(record, Uni_Exon, -6)
                    Strand = strand_orientation(record, Uni_Exon, Exon_Start, Exon_End)
                    if Strand == '+':
                        SS5_Seq, SS5_Mut_Seq = pos_strand_intron_mut_ss5_seq(record, Read_Genome_Fa, Exon_End)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS5_Seq, SS5_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                if intron_mut_nearby(record, Uni_Exon, 20):
                    Exon_Start, Exon_End = exon_range(record, Uni_Exon, 20)
                    Strand = strand_orientation(record, Uni_Exon, Exon_Start, Exon_End)
                    if Strand == '+':
                        SS3_Seq, SS3_Mut_Seq = pos_strand_intron_mut_ss3_seq(record, Read_Genome_Fa, Exon_Start)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS3_Seq, SS3_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                if intron_mut_nearby(record, Uni_Exon, 6):
                    Exon_Start, Exon_End = exon_range(record, Uni_Exon, 6)
                    Strand = strand_orientation(record, Uni_Exon, Exon_Start, Exon_End)
                    if Strand == '-':
                        SS5_Seq, SS5_Mut_Seq = neg_strand_intron_mut_ss5_seq(record, Read_Genome_Fa, Exon_Start)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS5_Seq, SS5_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
                if intron_mut_nearby(record, Uni_Exon, -20):
                    Exon_Start, Exon_End = exon_range(record, Uni_Exon, -20)
                    Strand = strand_orientation(record, Uni_Exon, Exon_Start, Exon_End)
                    if Strand == '-':
                        SS3_Seq, SS3_Mut_Seq = neg_strand_intron_mut_ss3_seq(record, Read_Genome_Fa, Exon_End)
                        CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, \
                        Post_Mut_Seq, splice_site_type = store_info(record, Mut_Location, Strand, SS3_Seq, SS3_Mut_Seq, Exon_Start, Exon_End, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End,
                                                  Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)
    write_info(pwd, sample, CHROM, POS, Mut_Place, Near_Exon_Start, Near_Exon_End, Exon_Strand, REF, ALT, Complement_REF, Complement_ALT, Pre_Mut_Seq, Post_Mut_Seq, Splice_Site_Type)

    # All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info = read_mut_info('mut_info.txt')
    # for i in range(Rep_Mut_Info.shape[0]):
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'POS')
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'Mut_Place')
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'REF')
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'ALT')
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'Complement_REF')
    #     combine_mut_info(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i, 'Complement_ALT')
    #     combine_post_mut_seq(All_Mut_Info, Uni_Mut_Info, Rep_Mut_Info, i)
    # Uni_Mut_Info.to_csv('mut_info.txt', index=False, sep='\t')

    c1 = "awk '$0 ~ /SS5/' "+pwd+'/'+sample+".mut_info.txt > "+pwd+'/'+sample+".SS5_Mut_Info"
    os.system(c1)
    c2 = "awk '$0 ~ /SS3/' "+pwd+'/'+sample+".mut_info.txt > "+pwd+'/'+sample+".SS3_Mut_Info"
    os.system(c2)
    c3 = "awk '{print $11}' "+pwd+'/'+sample+".SS5_Mut_Info > "+pwd+'/'+sample+".Pre_Mut_SS5"
    os.system(c3)
    c4 = "awk '{print $12}' "+pwd+'/'+sample+".SS5_Mut_Info > "+pwd+'/'+sample+".Post_Mut_SS5"
    os.system(c4)
    c5 = "awk '{print $11}' "+pwd+'/'+sample+".SS3_Mut_Info > "+pwd+'/'+sample+".Pre_Mut_SS3"
    os.system(c5)
    c6 = "awk '{print $12}' "+pwd+'/'+sample+".SS3_Mut_Info > "+pwd+'/'+sample+".Post_Mut_SS3"
    os.system(c6)
    c7 = "perl " + maxentscan_path + "/score5.pl "+pwd+'/'+sample+".Pre_Mut_SS5 > "+pwd+'/'+sample+".Pre_Mut_SS5_MaxEntScore"
    os.system(c7)
    c8 = "perl " + maxentscan_path + "/score5.pl "+pwd+'/'+sample+".Post_Mut_SS5 > "+pwd+'/'+sample+".Post_Mut_SS5_MaxEntScore"
    os.system(c8)
    c9 = "perl " + maxentscan_path + "/score3.pl "+pwd+'/'+sample+".Pre_Mut_SS3 > "+pwd+'/'+sample+".Pre_Mut_SS3_MaxEntScore"
    os.system(c9)
    c10 = "perl " + maxentscan_path + "/score3.pl "+pwd+'/'+sample+".Post_Mut_SS3 > "+pwd+'/'+sample+".Post_Mut_SS3_MaxEntScore"
    os.system(c10)
    os.system("awk '{print $2}' "+pwd+'/'+sample+".Pre_Mut_SS5_MaxEntScore > "+pwd+'/'+sample+".New_Pre_Mut_SS5_MaxEntScore")
    os.system("awk '{print $2}' "+pwd+'/'+sample+".Post_Mut_SS5_MaxEntScore > "+pwd+'/'+sample+".New_Post_Mut_SS5_MaxEntScore")
    os.system("paste "+pwd+'/'+sample+".SS5_Mut_Info "+pwd+'/'+sample+".New_Pre_Mut_SS5_MaxEntScore "+pwd+'/'+sample+".New_Post_Mut_SS5_MaxEntScore > "+pwd+'/'+sample+".SS5_Mut_Info_Final")
    os.system("awk '{print $15-$14}' "+pwd+'/'+sample+".SS5_Mut_Info_Final > "+pwd+'/'+sample+".SS5_Delta_Score")
    os.system("paste "+pwd+'/'+sample+".SS5_Mut_Info_Final "+pwd+'/'+sample+".SS5_Delta_Score > "+pwd+'/'+sample+".Real_SS5_Mut_Info_Final")
    os.system("awk '{print $2}' "+pwd+'/'+sample+".Pre_Mut_SS3_MaxEntScore > "+pwd+'/'+sample+".New_Pre_Mut_SS3_MaxEntScore")
    os.system("awk '{print $2}' "+pwd+'/'+sample+".Post_Mut_SS3_MaxEntScore > "+pwd+'/'+sample+".New_Post_Mut_SS3_MaxEntScore")
    os.system("paste "+pwd+'/'+sample+".SS3_Mut_Info "+pwd+'/'+sample+".New_Pre_Mut_SS3_MaxEntScore "+pwd+'/'+sample+".New_Post_Mut_SS3_MaxEntScore > "+pwd+'/'+sample+".SS3_Mut_Info_Final")
    os.system("awk '{print $15-$14}' "+pwd+'/'+sample+".SS3_Mut_Info_Final > "+pwd+'/'+sample+".SS3_Delta_Score")
    os.system("paste "+pwd+'/'+sample+".SS3_Mut_Info_Final "+pwd+'/'+sample+".SS3_Delta_Score > "+pwd+'/'+sample+".Real_SS3_Mut_Info_Final")
    os.system("printf '#CHROM\tPOS\tMut_Place\tNear_Exon_Start\tNear_Exon_End\tExon_Strand\tREF\tALT\tComplement_REF\tComplement_ALT\t"
              "Pre_Mut_Seq\tPost_Mut_Seq\tSplice_Site_Type\tRef_Seq_MaxEntScore\tAlt_Seq_MaxEntScore\tDelta_Score\n' > "+pwd+'/'+sample+".MaxEntRes.txt")
    os.system("cat "+pwd+'/'+sample+".Real_SS5_Mut_Info_Final >> "+pwd+'/'+sample+".MaxEntRes.txt")
    os.system("cat "+pwd+'/'+sample+".Real_SS3_Mut_Info_Final >> "+pwd+'/'+sample+".MaxEntRes.txt")
