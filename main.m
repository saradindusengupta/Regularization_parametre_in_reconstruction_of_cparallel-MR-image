function mainproject(abbbcd,sc)
k=abbbcd/max(abs(abbbcd(:)));
[subkspace,ACS,acqACS,acq,mask]  = micsplkspacesubsample(k,sc,3);
mask=mask(:,:,1);
[I_recon,Istk,Kstk,E]=main_Ad_spirit(k, [5,5], 20, sc,3);figure;plot(E);