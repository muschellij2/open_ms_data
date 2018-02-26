
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multiple Sclerosis Lesion Data

The data in this repository is from the Laboratory of Imaging
Technologies (<http://lit.fe.uni-lj.si/tools.php?lang=eng>).

## Cross-sectional data

The cross-sectional data is from the “3D MR image database of Multiple
Sclerosis patients with white matter lesion segmentations” section of
the website.

The database is released under the Creative-Commons Attribution (CC-BY)
license. Please cite the references below in any published work that
uses this
database.

![](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by.svg)

### References

The data is described in: Lesjak, Žiga, et al. “A novel public MR image
dataset of multiple sclerosis patients with lesion segmentations based
on multi-rater consensus.” Neuroinformatics (2017): 1-13.

### Demographics

| id        | age | sex | ms\_type | edss | criteria      |
| :-------- | --: | :-- | :------- | ---: | :------------ |
| patient01 |  31 | F   | RR       |  1.5 | McDonald 2010 |
| patient02 |  33 | M   | CIS      |  0.0 | McDonald 2010 |
| patient03 |  37 | F   | NA       |   NA | NA            |
| patient04 |  25 | M   | SP       |  6.5 | McDonald 2005 |
| patient05 |  33 | F   | RR       |  3.5 | McDonald 2005 |
| patient06 |  37 | F   | SP       |  4.0 | McDonald 2005 |
| patient07 |  53 | F   | RR       |  0.5 | McDonald 2010 |
| patient08 |  41 | M   | RR       |  5.0 | McDonald 2005 |
| patient09 |  40 | F   | RR       |  2.0 | McDonald 2010 |
| patient10 |  64 | F   | RR       |  2.0 | McDonald 2005 |
| patient11 |  29 | M   | RR       |  2.0 | McDonald 2010 |
| patient12 |  39 | F   | RR       |   NA | McDonald 2005 |
| patient13 |  26 | M   | RR       |  2.0 | McDonald 2010 |
| patient14 |  42 | M   | RR       |  4.0 | McDonald 2005 |
| patient15 |  57 | F   | PR       |  6.5 | McDonald 2005 |
| patient16 |  42 | F   | RR       |  4.0 | McDonald 2005 |
| patient17 |  27 | F   | RR       |  0.0 | McDonald 2010 |
| patient18 |  60 | F   | RR       |  2.0 | McDonald 2005 |
| patient19 |  47 | F   | RR       |  1.0 | McDonald 2005 |
| patient20 |  37 | F   | RR       |   NA | McDonald 2010 |
| patient21 |  33 | F   | RR       |  4.5 | McDonald 2005 |
| patient22 |  30 | F   | RR       |  2.0 | McDonald 2005 |
| patient23 |  39 | F   | RR       |  5.0 | McDonald 2005 |
| patient24 |  43 | M   | RR       |  0.0 | McDonald 2005 |
| patient25 |  35 | F   | RR       |  3.5 | McDonald 2005 |
| patient26 |  40 | F   | RR       |  2.0 | McDonald 2005 |
| patient27 |  39 | F   | RR       |  2.5 | McDonald 2005 |
| patient28 |  39 | F   | RR       |  0.0 | McDonald 2005 |
| patient29 |  26 | F   | CIS      |  3.0 | McDonald 2010 |
| patient30 |  54 | F   | RR       |  1.5 | McDonald 2005 |

### Raw Data

|           |                                                     |                                                                     |                                                |                                                      |                                                |
| :-------- | :-------------------------------------------------- | :------------------------------------------------------------------ | :--------------------------------------------- | :--------------------------------------------------- | :--------------------------------------------- |
| patient01 | [FLAIR](cross_sectional/raw/patient01/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient01/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient01/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient01/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient01/T2W.nii.gz) |
| patient02 | [FLAIR](cross_sectional/raw/patient02/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient02/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient02/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient02/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient02/T2W.nii.gz) |
| patient03 | [FLAIR](cross_sectional/raw/patient03/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient03/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient03/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient03/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient03/T2W.nii.gz) |
| patient04 | [FLAIR](cross_sectional/raw/patient04/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient04/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient04/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient04/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient04/T2W.nii.gz) |
| patient05 | [FLAIR](cross_sectional/raw/patient05/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient05/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient05/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient05/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient05/T2W.nii.gz) |
| patient06 | [FLAIR](cross_sectional/raw/patient06/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient06/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient06/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient06/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient06/T2W.nii.gz) |
| patient07 | [FLAIR](cross_sectional/raw/patient07/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient07/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient07/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient07/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient07/T2W.nii.gz) |
| patient08 | [FLAIR](cross_sectional/raw/patient08/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient08/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient08/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient08/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient08/T2W.nii.gz) |
| patient09 | [FLAIR](cross_sectional/raw/patient09/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient09/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient09/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient09/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient09/T2W.nii.gz) |
| patient10 | [FLAIR](cross_sectional/raw/patient10/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient10/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient10/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient10/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient10/T2W.nii.gz) |
| patient11 | [FLAIR](cross_sectional/raw/patient11/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient11/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient11/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient11/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient11/T2W.nii.gz) |
| patient12 | [FLAIR](cross_sectional/raw/patient12/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient12/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient12/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient12/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient12/T2W.nii.gz) |
| patient13 | [FLAIR](cross_sectional/raw/patient13/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient13/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient13/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient13/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient13/T2W.nii.gz) |
| patient14 | [FLAIR](cross_sectional/raw/patient14/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient14/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient14/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient14/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient14/T2W.nii.gz) |
| patient15 | [FLAIR](cross_sectional/raw/patient15/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient15/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient15/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient15/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient15/T2W.nii.gz) |
| patient16 | [FLAIR](cross_sectional/raw/patient16/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient16/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient16/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient16/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient16/T2W.nii.gz) |
| patient17 | [FLAIR](cross_sectional/raw/patient17/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient17/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient17/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient17/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient17/T2W.nii.gz) |
| patient18 | [FLAIR](cross_sectional/raw/patient18/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient18/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient18/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient18/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient18/T2W.nii.gz) |
| patient19 | [FLAIR](cross_sectional/raw/patient19/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient19/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient19/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient19/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient19/T2W.nii.gz) |
| patient20 | [FLAIR](cross_sectional/raw/patient20/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient20/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient20/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient20/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient20/T2W.nii.gz) |
| patient21 | [FLAIR](cross_sectional/raw/patient21/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient21/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient21/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient21/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient21/T2W.nii.gz) |
| patient22 | [FLAIR](cross_sectional/raw/patient22/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient22/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient22/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient22/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient22/T2W.nii.gz) |
| patient23 | [FLAIR](cross_sectional/raw/patient23/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient23/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient23/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient23/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient23/T2W.nii.gz) |
| patient24 | [FLAIR](cross_sectional/raw/patient24/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient24/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient24/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient24/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient24/T2W.nii.gz) |
| patient25 | [FLAIR](cross_sectional/raw/patient25/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient25/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient25/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient25/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient25/T2W.nii.gz) |
| patient26 | [FLAIR](cross_sectional/raw/patient26/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient26/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient26/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient26/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient26/T2W.nii.gz) |
| patient27 | [FLAIR](cross_sectional/raw/patient27/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient27/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient27/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient27/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient27/T2W.nii.gz) |
| patient28 | [FLAIR](cross_sectional/raw/patient28/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient28/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient28/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient28/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient28/T2W.nii.gz) |
| patient29 | [FLAIR](cross_sectional/raw/patient29/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient29/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient29/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient29/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient29/T2W.nii.gz) |
| patient30 | [FLAIR](cross_sectional/raw/patient30/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/raw/patient30/consensus_gt.nii.gz) | [T1](cross_sectional/raw/patient30/T1W.nii.gz) | [T1Post](cross_sectional/raw/patient30/T1WKS.nii.gz) | [T2](cross_sectional/raw/patient30/T2W.nii.gz) |

### Slightly Processed/Coregistered Data

These images were co-registered to the FLAIR and bias corrected and
include a brain mask (in the FLAIR
space.)

|           |                                                                         |                                                               |                                                                               |                                                          |                                                                |                                                          |
| :-------- | :---------------------------------------------------------------------- | :------------------------------------------------------------ | :---------------------------------------------------------------------------- | :------------------------------------------------------- | :------------------------------------------------------------- | :------------------------------------------------------- |
| patient01 | [Brain\_Mask](cross_sectional/coregistered//patient01/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient01/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient01/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient01/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient01/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient01/T2W.nii.gz) |
| patient02 | [Brain\_Mask](cross_sectional/coregistered//patient02/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient02/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient02/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient02/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient02/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient02/T2W.nii.gz) |
| patient03 | [Brain\_Mask](cross_sectional/coregistered//patient03/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient03/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient03/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient03/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient03/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient03/T2W.nii.gz) |
| patient04 | [Brain\_Mask](cross_sectional/coregistered//patient04/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient04/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient04/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient04/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient04/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient04/T2W.nii.gz) |
| patient05 | [Brain\_Mask](cross_sectional/coregistered//patient05/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient05/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient05/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient05/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient05/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient05/T2W.nii.gz) |
| patient06 | [Brain\_Mask](cross_sectional/coregistered//patient06/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient06/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient06/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient06/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient06/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient06/T2W.nii.gz) |
| patient07 | [Brain\_Mask](cross_sectional/coregistered//patient07/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient07/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient07/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient07/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient07/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient07/T2W.nii.gz) |
| patient08 | [Brain\_Mask](cross_sectional/coregistered//patient08/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient08/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient08/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient08/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient08/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient08/T2W.nii.gz) |
| patient09 | [Brain\_Mask](cross_sectional/coregistered//patient09/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient09/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient09/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient09/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient09/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient09/T2W.nii.gz) |
| patient10 | [Brain\_Mask](cross_sectional/coregistered//patient10/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient10/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient10/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient10/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient10/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient10/T2W.nii.gz) |
| patient11 | [Brain\_Mask](cross_sectional/coregistered//patient11/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient11/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient11/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient11/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient11/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient11/T2W.nii.gz) |
| patient12 | [Brain\_Mask](cross_sectional/coregistered//patient12/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient12/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient12/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient12/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient12/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient12/T2W.nii.gz) |
| patient13 | [Brain\_Mask](cross_sectional/coregistered//patient13/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient13/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient13/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient13/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient13/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient13/T2W.nii.gz) |
| patient14 | [Brain\_Mask](cross_sectional/coregistered//patient14/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient14/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient14/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient14/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient14/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient14/T2W.nii.gz) |
| patient15 | [Brain\_Mask](cross_sectional/coregistered//patient15/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient15/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient15/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient15/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient15/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient15/T2W.nii.gz) |
| patient16 | [Brain\_Mask](cross_sectional/coregistered//patient16/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient16/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient16/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient16/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient16/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient16/T2W.nii.gz) |
| patient17 | [Brain\_Mask](cross_sectional/coregistered//patient17/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient17/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient17/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient17/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient17/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient17/T2W.nii.gz) |
| patient18 | [Brain\_Mask](cross_sectional/coregistered//patient18/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient18/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient18/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient18/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient18/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient18/T2W.nii.gz) |
| patient19 | [Brain\_Mask](cross_sectional/coregistered//patient19/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient19/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient19/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient19/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient19/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient19/T2W.nii.gz) |
| patient20 | [Brain\_Mask](cross_sectional/coregistered//patient20/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient20/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient20/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient20/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient20/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient20/T2W.nii.gz) |
| patient21 | [Brain\_Mask](cross_sectional/coregistered//patient21/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient21/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient21/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient21/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient21/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient21/T2W.nii.gz) |
| patient22 | [Brain\_Mask](cross_sectional/coregistered//patient22/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient22/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient22/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient22/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient22/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient22/T2W.nii.gz) |
| patient23 | [Brain\_Mask](cross_sectional/coregistered//patient23/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient23/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient23/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient23/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient23/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient23/T2W.nii.gz) |
| patient24 | [Brain\_Mask](cross_sectional/coregistered//patient24/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient24/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient24/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient24/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient24/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient24/T2W.nii.gz) |
| patient25 | [Brain\_Mask](cross_sectional/coregistered//patient25/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient25/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient25/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient25/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient25/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient25/T2W.nii.gz) |
| patient26 | [Brain\_Mask](cross_sectional/coregistered//patient26/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient26/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient26/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient26/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient26/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient26/T2W.nii.gz) |
| patient27 | [Brain\_Mask](cross_sectional/coregistered//patient27/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient27/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient27/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient27/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient27/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient27/T2W.nii.gz) |
| patient28 | [Brain\_Mask](cross_sectional/coregistered//patient28/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient28/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient28/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient28/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient28/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient28/T2W.nii.gz) |
| patient29 | [Brain\_Mask](cross_sectional/coregistered//patient29/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient29/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient29/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient29/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient29/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient29/T2W.nii.gz) |
| patient30 | [Brain\_Mask](cross_sectional/coregistered//patient30/brainmask.nii.gz) | [FLAIR](cross_sectional/coregistered//patient30/FLAIR.nii.gz) | [Gold\_Standard](cross_sectional/coregistered//patient30/consensus_gt.nii.gz) | [T1](cross_sectional/coregistered//patient30/T1W.nii.gz) | [T1Post](cross_sectional/coregistered//patient30/T1WKS.nii.gz) | [T2](cross_sectional/coregistered//patient30/T2W.nii.gz) |

## Longitudinal data

The longitudinal data is from the “Longitudinal MR image database of
Multiple Sclerosis patients with white matter lesion change
segmentation” section of the website.

This archive contains longitudinal Magnetic Resonance (MR) images of
patients with Multiple Sclerosis (MS) with corresponding ground truth
segmentations of white matter lesion changes. Each patient has been
imaged twice on two separate occasions. The images are grouped intro
folders for each of the patients.

The database is released under the Creative-Commons Attribution (CC-BY)
license. Please cite the references below in any published work that
uses this
database.

![](https://mirrors.creativecommons.org/presskit/buttons/88x31/svg/by.svg)

### Description

This archive contains longitudinal Magnetic Resonance (MR) images of
Multiple Sclerosis (MS) patients with corresponding ground truth
segmentations of white matter lesion changes. Each patient has been
imaged twice on two separate occasions. The images are grouped into
folders for each of the patients. Each patient’s folder consists of: -
Co-registered and N4 corrected T1-weighted, T2-weighted and FLAIR images
for both MR studies - Brain mask - White matter lesion change mask - RAW
images in their original space - Intra-study transform parameters and
transform parameters from RAW to common space

### Demographics

| id        | sex | age | ms\_type | days\_between\_studies |
| :-------- | :-- | --: | :------- | ---------------------: |
| patient01 | F   |  20 | NA       |                    203 |
| patient02 | F   |  29 | RR       |                    283 |
| patient03 | M   |  19 | RR       |                    671 |
| patient04 | F   |  33 | NA       |                    155 |
| patient05 | F   |  35 | RR       |                    380 |
| patient06 | M   |  47 | RR       |                    473 |
| patient07 | F   |  40 | SP       |                    317 |
| patient08 | F   |  41 | RR       |                    305 |
| patient09 | M   |  20 | NA       |                    116 |
| patient10 | M   |  50 | NA       |                    236 |
| patient11 | F   |  31 | RR       |                    405 |
| patient12 | M   |  19 | NA       |                     81 |
| patient13 | F   |  49 | RR       |                    329 |
| patient14 | F   |  29 | RR       |                    723 |
| patient15 | F   |  41 | RR       |                    168 |
| patient16 | F   |  41 | RR       |                    164 |
| patient17 | F   |  25 | RR       |                    321 |
| patient18 | F   |  44 | RR       |                    524 |
| patient19 | F   |  19 | RR       |                    228 |
| patient20 | F   |  50 | RR       |                    454 |

### Raw Data

|           |   |                                                         |                                                          |                                                     |                                                     |
| :-------- | :- | :------------------------------------------------------ | :------------------------------------------------------- | :-------------------------------------------------- | :-------------------------------------------------- |
| patient01 | 1 | [Gold\_Standard](longitiudinal/raw/patient01/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient01/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient01/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient01/study1_T2W.nii.gz) |
| patient01 | 2 | [Gold\_Standard](longitiudinal/raw/patient01/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient01/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient01/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient01/study2_T2W.nii.gz) |
| patient02 | 1 | [Gold\_Standard](longitiudinal/raw/patient02/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient02/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient02/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient02/study1_T2W.nii.gz) |
| patient02 | 2 | [Gold\_Standard](longitiudinal/raw/patient02/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient02/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient02/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient02/study2_T2W.nii.gz) |
| patient03 | 1 | [Gold\_Standard](longitiudinal/raw/patient03/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient03/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient03/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient03/study1_T2W.nii.gz) |
| patient03 | 2 | [Gold\_Standard](longitiudinal/raw/patient03/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient03/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient03/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient03/study2_T2W.nii.gz) |
| patient04 | 1 | [Gold\_Standard](longitiudinal/raw/patient04/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient04/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient04/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient04/study1_T2W.nii.gz) |
| patient04 | 2 | [Gold\_Standard](longitiudinal/raw/patient04/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient04/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient04/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient04/study2_T2W.nii.gz) |
| patient05 | 1 | [Gold\_Standard](longitiudinal/raw/patient05/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient05/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient05/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient05/study1_T2W.nii.gz) |
| patient05 | 2 | [Gold\_Standard](longitiudinal/raw/patient05/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient05/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient05/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient05/study2_T2W.nii.gz) |
| patient06 | 1 | [Gold\_Standard](longitiudinal/raw/patient06/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient06/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient06/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient06/study1_T2W.nii.gz) |
| patient06 | 2 | [Gold\_Standard](longitiudinal/raw/patient06/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient06/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient06/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient06/study2_T2W.nii.gz) |
| patient07 | 1 | [Gold\_Standard](longitiudinal/raw/patient07/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient07/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient07/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient07/study1_T2W.nii.gz) |
| patient07 | 2 | [Gold\_Standard](longitiudinal/raw/patient07/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient07/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient07/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient07/study2_T2W.nii.gz) |
| patient08 | 1 | [Gold\_Standard](longitiudinal/raw/patient08/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient08/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient08/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient08/study1_T2W.nii.gz) |
| patient08 | 2 | [Gold\_Standard](longitiudinal/raw/patient08/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient08/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient08/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient08/study2_T2W.nii.gz) |
| patient09 | 1 | [Gold\_Standard](longitiudinal/raw/patient09/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient09/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient09/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient09/study1_T2W.nii.gz) |
| patient09 | 2 | [Gold\_Standard](longitiudinal/raw/patient09/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient09/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient09/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient09/study2_T2W.nii.gz) |
| patient10 | 1 | [Gold\_Standard](longitiudinal/raw/patient10/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient10/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient10/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient10/study1_T2W.nii.gz) |
| patient10 | 2 | [Gold\_Standard](longitiudinal/raw/patient10/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient10/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient10/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient10/study2_T2W.nii.gz) |
| patient11 | 1 | [Gold\_Standard](longitiudinal/raw/patient11/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient11/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient11/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient11/study1_T2W.nii.gz) |
| patient11 | 2 | [Gold\_Standard](longitiudinal/raw/patient11/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient11/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient11/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient11/study2_T2W.nii.gz) |
| patient12 | 1 | [Gold\_Standard](longitiudinal/raw/patient12/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient12/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient12/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient12/study1_T2W.nii.gz) |
| patient12 | 2 | [Gold\_Standard](longitiudinal/raw/patient12/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient12/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient12/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient12/study2_T2W.nii.gz) |
| patient13 | 1 | [Gold\_Standard](longitiudinal/raw/patient13/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient13/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient13/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient13/study1_T2W.nii.gz) |
| patient13 | 2 | [Gold\_Standard](longitiudinal/raw/patient13/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient13/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient13/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient13/study2_T2W.nii.gz) |
| patient14 | 1 | [Gold\_Standard](longitiudinal/raw/patient14/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient14/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient14/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient14/study1_T2W.nii.gz) |
| patient14 | 2 | [Gold\_Standard](longitiudinal/raw/patient14/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient14/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient14/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient14/study2_T2W.nii.gz) |
| patient15 | 1 | [Gold\_Standard](longitiudinal/raw/patient15/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient15/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient15/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient15/study1_T2W.nii.gz) |
| patient15 | 2 | [Gold\_Standard](longitiudinal/raw/patient15/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient15/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient15/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient15/study2_T2W.nii.gz) |
| patient16 | 1 | [Gold\_Standard](longitiudinal/raw/patient16/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient16/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient16/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient16/study1_T2W.nii.gz) |
| patient16 | 2 | [Gold\_Standard](longitiudinal/raw/patient16/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient16/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient16/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient16/study2_T2W.nii.gz) |
| patient17 | 1 | [Gold\_Standard](longitiudinal/raw/patient17/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient17/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient17/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient17/study1_T2W.nii.gz) |
| patient17 | 2 | [Gold\_Standard](longitiudinal/raw/patient17/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient17/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient17/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient17/study2_T2W.nii.gz) |
| patient18 | 1 | [Gold\_Standard](longitiudinal/raw/patient18/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient18/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient18/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient18/study1_T2W.nii.gz) |
| patient18 | 2 | [Gold\_Standard](longitiudinal/raw/patient18/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient18/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient18/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient18/study2_T2W.nii.gz) |
| patient19 | 1 | [Gold\_Standard](longitiudinal/raw/patient19/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient19/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient19/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient19/study1_T2W.nii.gz) |
| patient19 | 2 | [Gold\_Standard](longitiudinal/raw/patient19/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient19/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient19/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient19/study2_T2W.nii.gz) |
| patient20 | 1 | [Gold\_Standard](longitiudinal/raw/patient20/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient20/study1_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient20/study1_T1W.nii.gz) | [T2](longitiudinal/raw/patient20/study1_T2W.nii.gz) |
| patient20 | 2 | [Gold\_Standard](longitiudinal/raw/patient20/gt.nii.gz) | [FLAIR](longitiudinal/raw/patient20/study2_FLAIR.nii.gz) | [T1](longitiudinal/raw/patient20/study2_T1W.nii.gz) | [T2](longitiudinal/raw/patient20/study2_T2W.nii.gz) |

### Slightly Processed/Coregistered Data

|           |   |                                                                   |                                                                       |                                                                    |                                                               |                                                               |
| :-------- | :- | :---------------------------------------------------------------- | :-------------------------------------------------------------------- | :----------------------------------------------------------------- | :------------------------------------------------------------ | :------------------------------------------------------------ |
| patient01 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient01/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient01/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient01/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient01/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient01/study1_T2W.nii.gz) |
| patient01 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient01/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient01/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient01/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient01/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient01/study2_T2W.nii.gz) |
| patient02 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient02/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient02/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient02/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient02/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient02/study1_T2W.nii.gz) |
| patient02 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient02/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient02/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient02/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient02/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient02/study2_T2W.nii.gz) |
| patient03 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient03/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient03/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient03/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient03/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient03/study1_T2W.nii.gz) |
| patient03 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient03/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient03/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient03/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient03/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient03/study2_T2W.nii.gz) |
| patient04 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient04/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient04/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient04/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient04/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient04/study1_T2W.nii.gz) |
| patient04 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient04/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient04/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient04/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient04/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient04/study2_T2W.nii.gz) |
| patient05 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient05/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient05/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient05/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient05/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient05/study1_T2W.nii.gz) |
| patient05 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient05/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient05/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient05/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient05/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient05/study2_T2W.nii.gz) |
| patient06 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient06/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient06/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient06/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient06/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient06/study1_T2W.nii.gz) |
| patient06 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient06/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient06/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient06/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient06/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient06/study2_T2W.nii.gz) |
| patient07 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient07/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient07/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient07/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient07/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient07/study1_T2W.nii.gz) |
| patient07 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient07/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient07/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient07/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient07/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient07/study2_T2W.nii.gz) |
| patient08 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient08/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient08/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient08/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient08/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient08/study1_T2W.nii.gz) |
| patient08 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient08/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient08/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient08/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient08/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient08/study2_T2W.nii.gz) |
| patient09 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient09/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient09/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient09/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient09/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient09/study1_T2W.nii.gz) |
| patient09 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient09/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient09/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient09/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient09/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient09/study2_T2W.nii.gz) |
| patient10 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient10/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient10/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient10/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient10/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient10/study1_T2W.nii.gz) |
| patient10 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient10/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient10/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient10/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient10/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient10/study2_T2W.nii.gz) |
| patient11 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient11/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient11/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient11/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient11/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient11/study1_T2W.nii.gz) |
| patient11 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient11/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient11/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient11/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient11/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient11/study2_T2W.nii.gz) |
| patient12 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient12/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient12/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient12/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient12/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient12/study1_T2W.nii.gz) |
| patient12 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient12/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient12/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient12/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient12/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient12/study2_T2W.nii.gz) |
| patient13 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient13/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient13/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient13/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient13/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient13/study1_T2W.nii.gz) |
| patient13 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient13/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient13/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient13/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient13/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient13/study2_T2W.nii.gz) |
| patient14 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient14/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient14/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient14/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient14/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient14/study1_T2W.nii.gz) |
| patient14 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient14/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient14/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient14/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient14/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient14/study2_T2W.nii.gz) |
| patient15 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient15/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient15/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient15/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient15/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient15/study1_T2W.nii.gz) |
| patient15 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient15/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient15/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient15/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient15/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient15/study2_T2W.nii.gz) |
| patient16 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient16/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient16/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient16/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient16/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient16/study1_T2W.nii.gz) |
| patient16 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient16/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient16/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient16/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient16/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient16/study2_T2W.nii.gz) |
| patient17 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient17/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient17/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient17/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient17/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient17/study1_T2W.nii.gz) |
| patient17 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient17/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient17/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient17/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient17/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient17/study2_T2W.nii.gz) |
| patient18 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient18/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient18/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient18/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient18/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient18/study1_T2W.nii.gz) |
| patient18 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient18/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient18/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient18/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient18/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient18/study2_T2W.nii.gz) |
| patient19 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient19/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient19/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient19/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient19/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient19/study1_T2W.nii.gz) |
| patient19 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient19/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient19/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient19/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient19/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient19/study2_T2W.nii.gz) |
| patient20 | 1 | [Gold\_Standard](longitiudinal/coregistered//patient20/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient20/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient20/study1_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient20/study1_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient20/study1_T2W.nii.gz) |
| patient20 | 2 | [Gold\_Standard](longitiudinal/coregistered//patient20/gt.nii.gz) | [Brain\_Mask](longitiudinal/coregistered//patient20/brainmask.nii.gz) | [FLAIR](longitiudinal/coregistered//patient20/study2_FLAIR.nii.gz) | [T1](longitiudinal/coregistered//patient20/study2_T1W.nii.gz) | [T2](longitiudinal/coregistered//patient20/study2_T2W.nii.gz) |

### References

The data is described in: Lesjak, Žiga, Franjo Pernuš, Boštjan Likar,
and Žiga Špiclin. “Validation of White Matter Lesion Change Detection
Methods on a Novel Publicly Available MRI Image Database.”
Neuroinformatics (2016): 1-18.
