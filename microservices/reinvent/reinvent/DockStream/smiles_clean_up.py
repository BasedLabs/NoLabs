import argparse


def smiles_cleaner(input_smiles: str, output_smiles: str):
    with open(input_smiles, 'r') as f:
        with open(output_smiles, 'w+') as new_f:
            for line in f.readlines():
                fixed_line = line.split(' ')[0]
                new_f.write(f'{fixed_line}\n')


if __name__ == '__main__':
    # take user specified input parameters to run the benchmarking script
    parser = argparse.ArgumentParser(description='Implements entry point to SMILES clean-up')
    parser.add_argument('-input_smiles', type=str, required=True, help='The path to a SMILES file')
    parser.add_argument('-output_smiles', type=str, required=True, help='The output path for the cleaned-up SMILES file')
    args = parser.parse_args()

    smiles_cleaner(args.input_smiles, args.output_smiles)
