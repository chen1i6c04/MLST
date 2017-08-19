from Bio.Blast.Applications import NcbiblastnCommandline
import os, shutil
import pandas as pd
from Bio.Blast import NCBIXML

def database_dir(species):
    database = 'database/' + species
    return database

def species_gene(database):
    gene = sorted(set([i.split('.')[0] for i in os.listdir(database)]))
    return gene

def makeoutfolder(in_file, outpath):
    out_folder = os.path.join(outpath, in_file.split('/')[-1])
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)
    return out_folder

def blast(file, gene, out_folder, database):
    for i in gene:
        blastx_cline = NcbiblastnCommandline(query=file, db='{}'.format(os.path.join(database, i)), evalue=0.1, outfmt=5,
                                             out=os.path.join(out_folder, "{}.xml".format(i)))
        blastx_cline()

def mlst(gene, out_folder):
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
                            if hsp.align_length / alignment.length == 1 and hsp.gaps == 0:
                                st.append(alignment.hit_def.split('-')[-1])
                            else:
                                st.append(0)
                            locus.append(alignment.hit_def.split('-')[0])
                            allele.append(alignment.hit_def)
                            identity.append(round((hsp.identities / alignment.length) * 100, 2))
                            hps_length.append(hsp.align_length)
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

def result(st, allele_profile):
    if 0 in st:
        return 'unknown'
    else:
        st = [int(i) for i in st]
        allele_profiling = allele_profile.values.tolist()
        for i in allele_profiling:
            if st == i:
                return allele_profile.index[allele_profiling.index(i)]
            else:
                return 'new type'