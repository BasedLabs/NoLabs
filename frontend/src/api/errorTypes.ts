export type ErrorResponse = {
  errors: string[],
  error_code: number
}

export enum ErrorCodes {
  fix_pdb_error = 1,
  conformations_log_websocket_undefined = 2,
  amino_acid_solubility_run_error = 3,
  no_amino_acids = 4,
  experiment_id_not_found = 5,
  amino_acid_localisation_run_error = 6,
  protein_design_run_error = 7,
  protein_design_update_metadata_error = 8,
}
