import os
import pandas as pd

def collect(out_folder):
    all_result = []
    for i in os.listdir(out_folder):
        report = pd.read_csv(os.path.join(out_folder, i) + '/result.csv', index_col=0)
        all_result.append(report)

    df = pd.DataFrame()
    for i in all_result:
        df = pd.concat([df, i])

    df.to_csv('all_result.csv')
