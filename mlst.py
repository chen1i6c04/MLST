from Bio.Blast.Applications import NcbiblastnCommandline
import os, shutil
import pandas as pd
from Bio.Blast import NCBIXML


gene = ['AROC', 'DNAN', 'HEMD', 'HISD', 'PURE', 'SUCA', 'THRA']
database = 'database'
senterica = pd.read_csv('senterica.csv', index_col=0)

def main(file, outpath):
    file_name = file.split('/')[-1]
    makeoutfolder(file_name, outpath)
    blast(file)
    st = mlst()
    seq_type = result(st)
    df = pd.DataFrame()
    df['Assembly'] = [file_name.split('_')[0] + '_' + file_name.split('_')[1]]
    df['Sequence type'] = [seq_type]
    df.to_csv('{}.csv'.format(os.path.join(out_folder, 'result')), index=False)

def makeoutfolder(file_name, outpath):
    global out_folder
    out_folder = os.path.join(outpath, file_name.split('/')[-1])
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)

def blast(file):
    for i in gene:
        blastx_cline = NcbiblastnCommandline(query=file, db='{}'.format(os.path.join(database, i)), evalue=0.1, outfmt=5,
                                             out=os.path.join(out_folder, "{}.xml".format(i)))
        blastx_cline()


def mlst():
    df = pd.DataFrame()
    st = []
    locus = []
    allele = []
    identity = []
    hps_length = []
    allele_length = []
    gaps = []
    for i in gene:
        with open(os.path.join(out_folder, "{}.xml".format(i))) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for rec in blast_records:
                count = 0
                for alignment in rec.alignments:
                    if count != 0:
                        pass
                    else:
                        count += 1
                        for hsp in alignment.hsps:
                            if hsp.score / alignment.length == 1 and hsp.gaps == 0:
                                st.append(alignment.hit_def.split('-')[-1])
                            else:
                                st.append('unknown')
                            locus.append(alignment.hit_def.split('-')[0])
                            allele.append(alignment.hit_def)
                            identity.append("%.2f" % float(float(hsp.identities) / float(alignment.length) * 100))
                            hps_length.append(round(hsp.score))
                            allele_length.append(alignment.length)
                            gaps.append(hsp.gaps)


    df['Locus'] = locus
    df['Identity'] = identity
    df['HPS length'] = hps_length
    df['Allele length'] = allele_length
    df['Gaps'] = gaps
    df['Allele'] = allele
    df.to_csv('{}.csv'.format(os.path.join(out_folder, 'mlst_report')), index=False)
    return st

def result(st):
    if 'unknown' in st:
        return 'unknown'
    else:
        st = [int(i) for i in st]
        senterica_ST = senterica.values.tolist()
        for i in senterica_ST:
            if st == i:
                return senterica.index[senterica_ST.index(i)]