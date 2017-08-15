from argparse import ArgumentParser
import mlst, result_collect
import os
from concurrent.futures import ProcessPoolExecutor

parser = ArgumentParser()
parser.add_argument('-i', required=True, help='Path for sequence folder locations')
parser.add_argument('-o', required=True, help='Result of output path')
parser.add_argument('-t', type=int, default=2, required=False, help='Number of core be used')
parser.add_argument('-s', required=True, help='Species of bacterial')

def main():
    args = parser.parse_args()
    in_folder = args.i
    out_folder = args.o
    thread = args.t
    species = args.s

    path = []
    for i in os.listdir(in_folder):
        path.append((os.path.join(in_folder, i), out_folder, species))

    with ProcessPoolExecutor(thread) as executor:
        executor.map(multiple_work, path)

    result_collect.collect(out_folder)

def multiple_work(i):
    seq_file, out_dir, species = i
    mlst.main(seq_file, out_dir, species)


if __name__ == '__main__':
    main()