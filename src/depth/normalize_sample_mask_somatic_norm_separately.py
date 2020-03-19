import argparse
import pandas as pd


def mask_unknown(df, mask_table_path):
    tmp_mask_table = pd.read_csv(mask_table_path, dtype={'chr': str})
    df['match'] = df[["chr", "start"]].apply(lambda x: hash(tuple(x)), axis=1)
    tmp_mask_table['match'] = tmp_mask_table[["chr", "start"]].apply(lambda x: hash(tuple(x)), axis=1)
    df = df[(~df['match'].isin(tmp_mask_table['match'])) | ((df['match'].isin(tmp_mask_table['match'])) & ((df['chr'] == 'X') | (df['chr'] == 'Y')))]
    df.reset_index(drop=True, inplace=True)

    return df


def main(args=None):
    input_depth, mask_table_path = args.input_depth, args.mask_table_path

    sample_depth = pd.read_csv(
        input_depth,
        header=None,
        names=["index", "start", "depth", "chr"],
        dtype={'chr': str}
    )

    norm_mask_sample = mask_unknown(sample_depth, mask_table_path)
    mean_depth_somatic = norm_mask_sample[(norm_mask_sample['chr'] != 'X') & (norm_mask_sample['chr'] != 'Y')].depth.astype(float).mean()
    mean_depth_sex = norm_mask_sample[(norm_mask_sample['chr'] == 'X') | (norm_mask_sample['chr'] == 'Y')].depth.astype(float).mean()

    def normalizer(row):
        is_sex_chromosome = bool(row['chr'] != 'X' and row['chr'] != 'Y')
        denominator = float(mean_depth_somatic if is_sex_chromosome else mean_depth_sex)

        return float(row['depth']) * 4 / denominator

    norm = norm_mask_sample.apply(normalizer, axis=1)

    input_name = input_depth.rsplit(".", 1)[0]
    output_file_path = f'{input_name}_mask_somatic_norm.csv'
    # print(f'output_file_path: #{output_file_path}#')
    norm.to_csv(output_file_path, header=["x"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Normalize sample depth')
    parser.add_argument('input_depth', help='sample depth csv.', type=str)
    parser.add_argument('mask_table_path', help='mask table csv.', type=str)
    args = parser.parse_args()

    main(args)
