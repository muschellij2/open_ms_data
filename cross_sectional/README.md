## 3D MR image database of Multiple Sclerosis patients with white matter lesion segmentations



This archive contains Magnetic Resonance (MR) images of Multiple 
Sclerosis (MS) patients with corresponding consensus based ground
truth segmentations of white matter lesions. The images are grouped 
into folders for each of the patients. Each patient's folder consists of:

- Co-registered and bias corrected T1-weighted (T1W), contrast enhanced T1-weighted (T1WKS), T2-weighted (T2W) and FLAIR images 
- Brain mask 
- Consensus white matter lesion segmentations
- RAW images in their original space 
- Intra-study transform parameters

To preserve patients' privacy **all images have been defaced**. 
For more details please see the *References*.

*Published on: http://lit.fe.uni-lj.si/tools*



### Folder structure ###

- `patientXX/`
    - `patientXX_brainmask.nii.gz` (Brain mask)
    - `patientXX_consensus_gt.nii.gz` (Consensus white matter lesion segmentations)
    - `patientXX_modality.nii.gz` (Spatially co-registered and bias corrected T1, T1KS, T2, FLAIR images)
    - `patientXX_modality_to_FLAIR.txt` (Affine transform parameters to register images into FLAIR space)
    - `raw/`
        - `patientXX_modality.nii.gz` (RAW images in their original space)



### License ###

The database is released under the Creative-Commons Attribution (CC-BY) license. Please cite the references below in any published work that uses this database.



### References ###

- Lesjak, Žiga, Franjo Pernuš, Boštjan Likar, and  Žiga Špiclin. "A Novel Public MR Image Dataset of Multiple Sclerosis Patients With Lesion Segmentations Based on Multi-rater Consensus", Neuroinformatics (2017), DOI: 10.1007/s12021-017-9348-7
