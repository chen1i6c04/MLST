import pandas as pd
import os

def collect(out_path):
    all_result = []
    for i in os.listdir(out_path):
        report = pd.read_csv(os.path.join(out_path, i) + '/result.csv', index_col=0)
        all_result.append(report)

    df = pd.DataFrame()
    for i in all_result:
        df = pd.concat([df, i])

    df.to_csv('all_result.csv', index=False)