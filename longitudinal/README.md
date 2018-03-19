## Longitudinal MR image database of Multiple Sclerosis patients with white matter lesion change segmentation



This archive contains longitudinal Magnetic Resonance (MR) images 
of Multiple Sclerosis (MS) patients with corresponding ground 
truth segmentations of white matter lesion changes. Each patient 
has been imaged twice on two separate occasions. The images are 
grouped into folders for each of the patients. Each patient's 
folder consists of:
- Co-registered and N4 corrected T1-weighted, T2-weighted and FLAIR images for both MR studies 
- Brain mask 
- White matter lesion change mask
- RAW images in their original space 
- Intra-study transform parameters and transform parameters from RAW to common space

To preserve patients' privacy **all images have been defaced**. 
For more details please see the *References*.

*Published on: http://lit.fe.uni-lj.si/tools*



### Folder structure ###

- `patientX/`
    - `patientX_brainmask.nii.gz` (Brain mask)
    - `patientX_gt.nii.gz` (Lesion change mask)
    - `patientX_studyY_modalityreg.nii.gz` (Spatially co-registered and n4 corrected T1, T2, FLAIR images)
    - `patientX_studyY_FLAIR_to_common_space.txt` (Affine transform parameters to register original FLAIR images into common reference space)
    - `raw/`
        - `patientX_studyY_modality.nii.gz` (RAW images in their original space)
        - `patientX_studyY_modality_intrastudy_to_FLAIR.txt` (Intra-study affine transform parameters for registering T1 and T2 images into FLAIR space)




### License ###

The database is released under the Creative-Commons Attribution (CC-BY) license. Please cite the references below in any published work that uses this database.



### References ###

- Lesjak, Žiga, Franjo Pernuš, Boštjan Likar, and Žiga Špiclin. "Validation of White-Matter Lesion Change Detection Methods on a Novel Publicly Available MRI Image Database." Neuroinformatics (2016): 1-18.
