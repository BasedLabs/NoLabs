#!/usr/bin/env python
#  coding=utf-8

import os
import pandas as pd
import argparse

import rdkit.Chem as Chem

from dockstream.utils.smiles import to_smiles


if __name__ == "__main__":

    # get the input parameters and parse them
    parser = argparse.ArgumentParser(description="Implements simple translator taking an SDF file and spitting out SMILES.")
    parser.add_argument("-sdf", type=str, default=None, help="A path a SDF file.")
    parser.add_argument("-smi", type=str, default=None, required=False, help="A path an output text file.")
    parser.add_argument("-csv", type=str, default=None, required=False, help="A path an output CSV file.")
    parser.add_argument("-keep_stereo", action="store_true", help="If set, exported SMILES contain stereo-information.")
    parser.add_argument("-tags2columns", type=str, nargs='+', default=None, required=False,
                        help="A list of strings for which tags should be transformed into columns.")
    args = parser.parse_args()

    if args.sdf is None or not os.path.isfile(args.sdf):
        raise Exception("Parameter \"-sdf\" must be a relative or absolute path to valid sdf file.")
    if args.smi is None and args.csv is None:
        raise Exception("At least one of the \"-smi\" or \"-csv\" output paths must be set.")

    molecules = []
    for mol in Chem.SDMolSupplier(args.sdf):
        if mol is None:
            continue
        molecules.append(mol)

    # write out
    # ---------
    if args.smi is not None:
        with open(args.smi, 'w') as smi_file:
            for mol in molecules:
                smi_file.write(to_smiles(mol, isomericSmiles=args.keep_stereo) + "\n")

    if args.csv is not None:
        data_buffer = []
        columns = ["Name", "SMILES"]
        tags2columns = []
        if args.tags2columns is not None:
            tags2columns = args.tags2columns
            columns = columns + tags2columns
        for mol in molecules:
            # add default columns for this row
            row = [mol.GetProp("_Name"),
                   to_smiles(mol, isomericSmiles=args.keep_stereo)]

            # add selected columns for this row (if specified)
            for tag in tags2columns:
                try:
                    row.append(mol.GetProp(tag))
                except KeyError:
                    row.append(None)

            data_buffer.append(row)
        df_writeout = pd.DataFrame(data_buffer, columns=columns)
        df_writeout.to_csv(path_or_buf=args.csv,
                           sep=',',
                           na_rep='',
                           header=True,
                           index=False,
                           mode='w',
                           quoting=None)
