# Regularization_parametre_in_reconstruction_of_cparallel-MR-image
Analyzing various regularization parameter  for k-space based parallel MR image reconstruction

#Undersampled data
How to Run
```
load data*
micsplkspacesubsample(k,nACS,s)
mask=mask(:,:,1)
kspace=kspace (varies) /max(abs(Kspace(:)))
main_Ad_spirit(k,[3,3],20,mask)
main_Ad_grappa(k,[3,3],20,mask)
```
Documentation here https://github.com/saradindusengupta/Regularization_parametre_in_reconstruction_of_cparallel-MR-image/blob/master/Document.pdf

s=2 and 3 , nACS = 32,24 and 18 variably.

*dataset is available at large amount in 3D raw data in high volume , not possible to push. The name of the dataset varies from
FLAIR_17
T2Weighted7
2DacqSWI3slice
GEphantom
T2spine_test2_kspace
Spine3
FLAIR_data*
