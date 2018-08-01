#############################################
# Trying smri processing
#############################################
cd ${structural}/open_ms/open_ms_data/programs; 

# new may need upwards of 70G
Rnosave process_cs_smri.R -t 1-30 -N PROCESS \
	-l mem_free=10G,h_vmem=11G

Rnosave create_predictor_df.R -N DF_MODEL \
	-l mem_free=20G,h_vmem=21G -hold_jid PROCESS

Rnosave open_ms_fit_model.R -t 3-4 -N MS_MODEL \
	-l mem_free=100G,h_vmem=101G

Rnosave predict_lesion.R -t 3-4 -N PREDICT \
	-l mem_free=50,h_vmem=51G -hold_jid_ad MS_MODEL

Rnosave create_prediction_df.R -t 1-4 -N OUTPREDICT \
	-hold_jid_ad PREDICT

Rnosave calculate_dice.R -N DICE

# Potentially exclude brain stem?
Rnosave plot_image_segmentations.R -N PLOT_RED

###################################
# OASIS Pipeline
###################################
Rnosave run_oasis.R -t 1-30 -N OASIS_PREPROC

Rnosave open_ms_fit_oasis.R -N OASIS_MODEL \
	-hold_jid OASIS_PREPROC

# only needed for #2 - increased mem
Rnosave predict_oasis_lesion.R -t 1-3 -N OASIS_PRED \
	-hold_jid OASIS_MODEL 
	# -l mem_free=30,h_vmem=31G 

Rnosave oasis_create_prediction_df.R -t 1-3 \
	-N OASIS_OUTPREDICT -hold_jid_ad OASIS_PRED


Rnosave oasis_calculate_dice.R -N OASIS_DICE


