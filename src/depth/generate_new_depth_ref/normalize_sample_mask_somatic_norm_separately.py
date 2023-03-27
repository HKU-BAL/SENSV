import argparse
import pandas as pd
import os

def mask_unknown(df, mask_table_path):
    tmp_mask_table = pd.read_csv(mask_table_path, dtype={'chr': str})
    df['match'] = df[["chr", "start"]].apply(lambda x: hash(tuple(x)), axis=1)
    tmp_mask_table['match'] = tmp_mask_table[["chr", "start"]].apply(lambda x: hash(tuple(x)), axis=1)
    df = df[(~df['match'].isin(tmp_mask_table['match'])) | ( (df['match'].isin(tmp_mask_table['match'])) &  ((df['chr']=='X')|(df['chr']=='Y')) )]
    df.reset_index(drop=True, inplace=True)
    return df


def main(args=None):
    parser = argparse.ArgumentParser(description='Normalize sample depth')
    parser.add_argument('input_depth', help='sample depth csv.')
    parser.add_argument('mask_table_path', help='mask table csv.')  
    args = parser.parse_args()
    sample_depth = pd.read_csv(args.input_depth, header=None, names=["index", "start", "depth", "chr"],
                               dtype={'chr': str})
    norm_mask_sample = mask_unknown(sample_depth, args.mask_table_path)
    mean_depth_somatic = norm_mask_sample[(norm_mask_sample['chr']!='X')&(norm_mask_sample['chr']!='Y')].depth.astype(float).mean()
    mean_depth_sex = norm_mask_sample[(norm_mask_sample['chr']=='X')|(norm_mask_sample['chr']=='Y')].depth.astype(float).mean()
    normalizer = lambda row: float(row['depth']) * 4 / mean_depth_somatic if (row['chr']!='X' and row['chr']!='Y') else float(row['depth']) * 4 / mean_depth_sex
    
    #print(norm_mask_sample.apply(normalizer,axis=1))
    #print(norm_mask_sample['depth'])
    norm=norm_mask_sample.apply(normalizer,axis=1)
    #print(norm_mask_sample['norm'])
    #print(norm_mask_sample['depth'])

    input_name = args.input_depth.rsplit(".",1)[0]
    norm.to_csv(os.path.join(input_name + '_mask_somatic_norm.csv'), header=["x"])


if __name__ == "__main__":
    main()
