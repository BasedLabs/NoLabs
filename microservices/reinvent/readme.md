Running reinvent4 with docking contains of multiple steps

Start from this file
https://github.com/MolecularAI/REINVENT4/blob/main/configs/toml/scoring_components_example.toml

there is an attachment with .zip archive of files with examples (but they are with Glide which is paid software)
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5#Sec23

You will need some outsource DockStream plugin, like autodock vina

1) generate a pdbqt file using this tutorial for receptor
https://autodock-vina.readthedocs.io/en/latest/docking_basic.html
You can install tools from conda https://anaconda.org/hcc/adfr-suite (python 2.7!!!!)
before that install conda-forge::xorg-libice
2) You can find a /dockstream_rl_direct_uncs.json file in archive from biomedcentral.com. Basically docking contains from two sections
"ligand_preparation" and "docking_runs". You'll need to adjust both of them
3) Ligand preparations - for DockStream use RDKit - https://github.com/MolecularAI/DockStream/blob/master/examples/docking/parallelized_AutodockVina_enumerated.json
4) You have to merge the file from the point 2 with the point 3 and prepared pdbqt file
5) Change the constrants here AutodockVina_enums.py