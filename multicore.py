import os
import mlst
import pandas as pd

def run(commend):
    in_path, out_path, species = commend
    database = mlst.database_dir(species)
    gene = mlst.species_gene(database)
    profile = mlst.allele_profile(species)
    out_folder = mlst.makeoutfolder(in_path, out_path)
    mlst.blast(in_path, gene, out_folder, database)
    st = mlst.mlst(gene, out_folder)
    report = mlst.result(st, profile)
    df = pd.DataFrame({'Assembly': [in_path.split('/')[-1]], 'Sequence': [report]})
    df.to_csv(os.path.join(out_folder, 'result.csv'))
