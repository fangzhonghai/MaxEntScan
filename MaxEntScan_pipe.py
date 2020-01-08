# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import sys
import os
import yaml
import glob
import time
pd.options.mode.chained_assignment = None


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    MaxEntScan::score5ss for human 5' splice sites AND MaxEntScan::score3ss for human 3' splice sites
    This program Packaged the 'MaxEntScan'.
    
    python MaxEntScan/MaxEntScan_pipe.py -a MaxEntScan/example/NA12878.txt -s sampleID -p /path/to/work \
    -s sampleID --pro project_in_sge --cols "#Chr,Stop,Ref,Call"
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def split_anno_df(anno_df, columns):
    if len(anno_df[~anno_df[columns[0]].str.startswith('chr')]):
        anno_df[columns[0]] = 'chr' + anno_df[columns[0]]
    anno_df.loc[anno_df[columns[0]] == 'chrMT', columns[0]] = 'chrM_NC_012920.1'
    anno_df_snp = anno_df[(anno_df[columns[2]].map(len) == 1) & (anno_df[columns[3]].map(len) == 1) & (anno_df[columns[2]] != '.') & (anno_df[columns[3]] != '.')].copy()
    anno_df_snp_dic = {}
    for chrom in anno_df_snp[columns[0]].unique():
        anno_df_snp_dic[chrom] = anno_df_snp[anno_df_snp[columns[0]] == chrom]
    return anno_df_snp_dic


def anno_dic2bed(anno_dic, work_path, sample_name):
    for key in anno_dic.keys():
        anno_key = anno_dic[key]
        anno_key.to_csv(work_path + "/" + sample_name + "_" + key + ".bed", index=False, sep='\t')


def write_mes_sh(anno_dic, work_path, sample_name, config):
    config_yaml = yaml.load(open(config), Loader=yaml.FullLoader)
    python = config_yaml['python']
    py = config_yaml['MaxEntScan_py']
    for key in anno_dic.keys():
        fp = open(work_path + "/" + sample_name + "_" + key + ".mes.sh", "w")
        shell = '''#!/bin/bash
{python} {py} -s {sample_name}_{key} -c {config} -v {work_path}/{sample_name}_{key}.bed -p {work_path}
if [ $? -ne 0 ]; then
        echo "{sample_name} {key} MaxEntScan failed."
        exit 1
else
        touch {work_path}/{sample_name}_{key}.mes.sh.check
fi
'''.format(**locals())
        fp.write(shell)
        fp.close()


def run_mes_sh(anno_dic, work_path, sample_name, project):
    for key in anno_dic.keys():
        qsub_sh = "qsub -cwd -l vf=1G,p=1 -P " + project + " " + work_path + "/" + sample_name + "_" + key + ".mes.sh"
        staus = os.system(qsub_sh)
        if staus != 0:
            sys.exit(1)


def wait_mes_res(anno_dic, work_path, sample_name):
    while 1:
        mes_checks = glob.glob(work_path + "/" + sample_name + "_*.mes.sh.check")
        if len(anno_dic.keys()) == len(mes_checks):
            break
        time.sleep(60)


def merge_mes_res(anno_dic, work_path, sample_name):
    mes_res_df = pd.DataFrame()
    for key in anno_dic.keys():
        df = pd.read_csv(work_path + "/" + sample_name + "_" + key + ".MaxEntRes.txt", sep='\t')
        mes_res_df = mes_res_df.append(df)
    mes_res_df['Relative_MaxEntScore'] = abs(mes_res_df['Delta_Score'])/mes_res_df['Ref_Seq_MaxEntScore']
    mes_res_df.to_csv(work_path + "/" + sample_name + ".MaxEntRes.txt", index=False, sep='\t')
    return mes_res_df


def merge_mes_res_2anno(mes_res_df, work_path, sample_name, anno_df, columns):
    mes_res = mes_res_df[['#CHROM', 'POS', 'REF', 'ALT', 'Splice_Site_Type', 'Relative_MaxEntScore']].copy()
    mes_res.columns = columns + ['Splice_Site_Type', 'Relative_MaxEntScore']
    if len(anno_df[~anno_df[columns[0]].str.startswith('chr')]):
        mes_res.loc[mes_res[columns[0]] == 'chrM_NC_012920.1', columns[0]] = 'chrMT'
        mes_res[columns[0]] = mes_res[columns[0]].str.replace('^chr', '')
    merged_mes_res = pd.merge(anno_df, mes_res, on=columns, how='left')
    merged_mes_res.drop_duplicates(inplace=True)
    merged_mes_res.to_csv(work_path + "/" + sample_name + ".MaxEntRes.merge.txt", index=False, sep='\t')


def write_rm_sh(work_path, sample_name):
    fp = open(work_path + "/" + sample_name + ".rm.sh", "w")
    shell = '''#!/bin/bash
rm -f {work_path}/{sample_name}_*.*SS*
'''.format(**locals())
    fp.write(shell)
    fp.close()


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-s', '--sample', dest='sample', default='none', type='string')
    parser.add_option('-c', '--config', dest='config', default='none', type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default='none', type='string')
    parser.add_option('-a', '--anno', dest='anno', default='none', type='string')
    parser.add_option('--pro', dest='pro', default='P18Z15000N0143', type='string')
    parser.add_option('--skip_row', dest='skip_row', default=0, type=int)
    parser.add_option('--cols', dest='cols', default='#Chr,Stop,Ref,Call', type='string')
    (opts, args) = parser.parse_args()
    sample = opts.sample
    config = opts.config
    pwd = opts.pwd
    anno = opts.anno
    pro = opts.pro
    skip_row = opts.skip_row
    cols = opts.cols
    cols_list = cols.split(',')
    anno_DF = pd.read_csv(anno, sep='\t', skiprows=range(skip_row), dtype={cols_list[0]: 'str'}, low_memory=False)
    anno_DF_part = anno_DF[cols_list].copy()
    anno_DF_part_snp_dic = split_anno_df(anno_DF_part, cols_list)
    anno_dic2bed(anno_DF_part_snp_dic, pwd, sample)
    write_mes_sh(anno_DF_part_snp_dic, pwd, sample, config)
    run_mes_sh(anno_DF_part_snp_dic, pwd, sample, pro)
    wait_mes_res(anno_DF_part_snp_dic, pwd, sample)
    mes_res_DF = merge_mes_res(anno_DF_part_snp_dic, pwd, sample)
    merge_mes_res_2anno(mes_res_DF, pwd, sample, anno_DF, cols_list)
    write_rm_sh(pwd, sample)
