cd /Volumes/DATA_LOCAL/Projects/open_ms_data/cross_sectional/processed
id=patient08
while read id; do
	echo ${id} ;
	mkdir -p ${id} ;
	iddir="${open_ms}/cross_sectional/raw/${id}/prenorm" ;
	hdr="$jdenig:${iddir}"
	rsync --progress "${hdr}"/*N4_noneck_reduced_winsor_regtoFLAIR_brain_N4_regtoMNI.nii.gz ${id} ;
	rsync --progress "${hdr}"/*N4_noneck_reduced_winsor_regtoFLAIR_brain_N4_regtoMNI_quantile.nii.gz ${id} ;
	rsync --progress "${hdr}"/*N4_noneck_reduced_winsor_regtoFLAIR_brain_N4_regtoMNI_trimmedz.nii.gz ${id} ;
	rsync --progress "${hdr}"/GOLD*N4_noneck_reduced_winsor_regtoFLAIR_regtoMNI.nii.gz ${id} ;
done < patients.txt