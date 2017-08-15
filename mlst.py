from Bio.Blast.Applications import NcbiblastnCommandline
import os, shutil
import pandas as pd
from Bio.Blast import NCBIXML


# gene = ['AROC', 'DNAN', 'HEMD', 'HISD', 'PURE', 'SUCA', 'THRA']
# database = 'database'
# senterica = pd.read_csv('senterica.csv', index_col=0)

def main(file, outpath, species):
    file_name = file.split('/')[-1]

    # === search database of bacterial ===
    database_gene(species)

    # === make out folder ===
    makeoutfolder(file_name, outpath)

    # === run blast ===
    blast(file)

    # === run mlst ===
    st = mlst()

    # === print result ===
    seq_type = result(st)

    df = pd.DataFrame()
    df['Assembly'] = [file_name.split('_')[0] + '_' + file_name.split('_')[1]]
    df['Sequence type'] = [seq_type]
    df.to_csv('{}.csv'.format(os.path.join(out_folder, 'result')), index=False)

def database_gene(species):
    global gene
    global database
    global allele_profile
    database = 'database/' + species
    gene = sorted(set([i.split('.')[0] for i in os.listdir(database)]))
    allele_profile = pd.read_csv('{}.csv'.format(species), index_col=0)

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
                                st.append(0)
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
    if 0 in st:
        return 'unknown'
    else:
        st = [int(i) for i in st]
        allele_profiling = allele_profile.values.tolist()
        for i in allele_profiling:
            if st == i:
                return allele_profile.index[allele_profiling.index(i)]