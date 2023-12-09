protein_design_root=src/ai/custom_models/protein_design
rfdiff_root=src/ai/custom_models/protein_design/RFdiffusion

# Install RFdiffusion
#git clone https://github.com/RosettaCommons/RFdiffusion.git $protein_design_root/RFdiffusion

#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt

# Optional:
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt

# original structure prediction weights
#wget -P $rfdiff_root/models http://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt

conda env create -f $rfdiff_root/env/SE3nv.yml

conda activate SE3nv
pip install --no-cache-dir -r $rfdiff_root/env/SE3Transformer/requirements.txt
python $rfdiff_root/env/SE3Transformer/setup.py install
pip install -e $rfdiff_root # install the rfdiffusion module from the root of the repository

conda run -n SE3nv pip install flask --upgrade