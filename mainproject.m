k=K/max(abs(K(:)));
[subkspace,ACS,acqACS,acq,mask]  = micsplkspacesubsample(k,26,3);
mask=mask(:,:,1);
[I_recon,Istk,Kstk,E]=main_Ad_spirit(k, [5,5],8,26,3);figure;plot(E);
save result@spine3;
clear all;