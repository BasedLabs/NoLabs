# Changelog of `DockStream`

## 1.0.0 - 2021-08-02 (RC)
### Added
- Minimal environment definitions.

### Fixed
- Fixed bug in `Schrodinger`/`LigPrep` sublist generation ('0' was interpreted as '1').

### Internal
- Added additional unit test for `Schrodinger`/`LigPrep` CLI application.

## 0.2.2 - 2020-04-30
### Added
- Extended input overwrite capabilities for `docker.py` (CSV specification).
- New parsing capability to load SDF input with `DockStream` nomenclature to ensure enumerations are properly assigned.
- Exposed additional command-line parameters for `Schrodinger`/`LigPrep`.
- Added execution of `Schrodinger`/LigPrep` on AWS ("relative paths" bug on Schrodinger's side - workaround).

### Fixed
- Hotfix to improve `Corina` result parsing.

### Internal
- Show-case of how to execute `Schrodinger`/`LigPrep` on AWS including a token guard.
- Show-case of how to execute `Schrodinger`/`Glide` on AWS including a token guard.

## 0.2.1 - 2021-04-16
### Added
- Additional measures to protect agains failing ligands at the embedding stage.
- Added `OpenEye`/`OMEGA` embedder.

### Fixed
- Hotfix to increase stability of `RDkit` ligand preparator.

### Internal
- Added unit test to cover changes in parameters for `LigPrep` when using `EPIK`.
- Improved logging for embedding stage.
- Reverted file system handling for `Glide`.

## 0.2.0 - 2021-03-31
### Added
- New backend `Hybrid` (`OpenEye`) added (will replace the API version shortly).
- Added stereo-isomer support using `RDkit`.
- Added stereo-isomer support using `Corina`.
- Exposed `Corina`'s "-d" options.
- Added support for `Ligprep`'s filtering capabilities.
- Enhanced stability for parallel execution mode for `Hybrid`, `Gold` and `Glide`.
- Benchmarking script added.
- Analysis script added.

### Fixed
- Fixed instability with protonations when using `RDkit` ligand embedding.
- Fixed logging bug within `Ligprep`/`Schrodinger`.
- Hot-fix for problem with CSV write-out when no pose was accepted.
- Improved logging for `Ligprep`/`Schrodinger`.

### Internal
- Introduced stereo-enumeration factory class.
- Adde `pydantic`-based parsing.
- Optimized log messages.
- Removed most versions from environment specifications.

## 0.1.4 - 2020-11-03
### Added
- Added support for parallelization for `Ligprep`/`Schrodinger` embedding.
- Added `best_per_ligand` write-out mode for conformers.
- Added CSV-option for entry point `sdf2smiles.py`.

### Fixed
- Fixed bug for `OpenEye` backend, when `max_compounds_per_subjob` was set to 0.
- Fixed bug in `best_per_enumeration` write-out mode for CSV results.
- Fixed irregularity with `Gold` conformer ordering (fitness versus score).
- Fixed naming bug in `OpenEye` backend.

### Internal
- Extended `AutoDock Vina` unit tests for tautomer compounds.
- Small update of `Ligprep`/`Schrodinger` unit tests.
- Updated docking backend unit tests to fully cover enumerations.
- Added unit test for `Glide`/`Schrodinger` constraints.

## 0.1.3 - 2020-09-29
### Added
- Added `OpenBabel`/`AutoDock Vina` target preparation (generating `PDBQT` files).
- Added `OpenBabel`/`AutoDock Vina` box extraction (based on XYZ coordinate ranges of template ligand).
- Added `AutoDock Vina` backend and result parser.
- Added possibility to specify binary path for external binary executions.

### Fixed
- Fixed issue with logfile logging for `Glide` when time was exceeded.
- Fixed issue with `OpenBabel` binary call when environment was not loaded.

### Internal
- Implemented hard over-write of `POSE_OUTTYPE` to be the only supported type.
- Redesign of internal dictionary usage.
- Addition of unit tests.
- Added logging the version number (and started tagging versions).
- Changed internal structure to be "package ready".
- Improved the tag adding (and extended it for the ligand preparation step).

## 0.1.2 - 2020-09-09
### Added
- Integration of a "progress bar" to docking jobs.
- Added support for parameter `max_compounds_per_subjob` to reign in (sub-)lists too long (especially with `Glide`).
- Added option to only output the best scores per enumeration.
- Added `ligprep` to available ligand embedding techniques.

### Fixed
- Fixed issues with "internal alignment" and integrated a fail-safe version.
- Made all receptor paths `lists` rather than simple strings to streamline the interface prior to the implementation of `ensemble docking`.

### Internal
- Clean-up of "OpenEye" result parser.
- Refactored result parsers.
- Improved `Gold` feedback for docking.
- Restructuring of internal "Ligand" handling.
- Refactored "docker.py" entry point.
- Added result parser output checks to respective unit tests.
- Refactored some methods of the ligand preparation tools.
- Updated example configuration files.

## 0.1.1 - 2020-08-17
### Added
- Added possibility to change the logfile path.
- Added parameter to change the time limit per compound for `Glide` docking.
- Added support to set a prefix for the output files.
- Refactored and extended tagging system for all backends (adds "smiles" and "original_smiles" now).
- Added option to only output the best poses per enumeration.
- Added transformation support (using SMIRKS and `OpenEye`) for ligand preparation.
- Added `Glide`/`Schrodinger` token guard.
- Added support of "SDF" files as input.
- Added support of arbitrary names to the parsing of CSV files as input.
- Added "-debug" parameter to entry points.

### Fixed
- Critical fix for score aggregation if enumerations were used.
- Made "Corina" embedding much more stable.
- Several minor bug fixes in the logging write-out.

### Internal
- Added "ChangeLog.md".
- Changed entry point structure (refactored and harmonized).
- Fixed issue with "min"/"max" docking scoring directions (Gold/CCDC backend).
