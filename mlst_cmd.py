from argparse import ArgumentParser
import mlst, multicore, total_report
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

parser = ArgumentParser()
parser.add_argument('-i', required=True, help='Path for sequence folder locations')
parser.add_argument('-o', required=True, help='Result of output path')
parser.add_argument('-t', type=int, default=2, required=False, help='Number of core be used')
parser.add_argument('-s', required=True, help='Species of bacterial')

def main():
    args = parser.parse_args()
    in_path = args.i
    out_path = args.o
    thread = args.t
    species = args.s
    allele_profile = pd.read_csv('{}.csv'.format(species), index_col=0)

    if os.path.isdir(in_path) == True:
        commend = []
        for i in os.listdir(in_path):
            in_dir = os.path.join(in_path, i)
            commend.append((in_dir, out_path, species, allele_profile))
        with ProcessPoolExecutor(thread) as executor:
            executor.map(multicore.run, commend)
        total_report.collect(out_path)

    else:
        database = mlst.database_dir(species)
        gene = mlst.species_gene(database)
        out_folder = mlst.makeoutfolder(in_path, out_path)
        mlst.blast(in_path, gene, out_folder, database)
        st = mlst.mlst(gene, out_folder)
        report = mlst.result(st, allele_profile)
        df = pd.DataFrame({'Assembly': [in_path.split('/')[-1]], 'Sequence': [report]})
        df.to_csv(out_folder + '/result.csv')

if __name__ == '__main__':
    main()