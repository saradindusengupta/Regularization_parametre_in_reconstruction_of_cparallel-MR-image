function[I_recon,Istk,Kstk,E]=main_Ad_spirit(DATA, kSize, nIter, mask)

% INPUTS
% kSize       : SPIRiT kernel size
% nIter       : number of iteration; phantom requires twice as much as the brain.
% mask        : mask can be uniform or random
% lambda      : Tykhonov regularization in the calibration
% T           : Wavelet soft-thresholding regularization in the reconstruction

% OUTPUTS
% I_recon    : Reconstructed Image (SoS combined)
% Istk       : Reconstructed image (Coil-wise)
% Kstk       : Reconstructed k-space (Coil-wise)
% Recon_err  : Reconstruction Error

bet = 1; 
im  = ifft2c(DATA);

%-----------------------UnderSampling--------------------------
[DATA, CalibSize, scale_fctr]=get_SampledData(DATA,mask);

% im_dc = ifft2c(DATAcomp);
im = im/scale_fctr;

%-----------------------Perform Calibration--------------------------

[GOP]=tsvdgcvCalib(DATA,CalibSize,kSize);

%--------------------------Reconstruction------------------------

[Kstk,E] = AdSPIRiT_recon(DATA, GOP, nIter, DATA, im,bet);
%  [Kstk, E] = AdWSPIRiT_recon2(DATA, GOP, nIter, DATA, T, im, bet);
Istk = ifft2c(Kstk);

I_recon=sos(Istk);




%------------------- calling function ----------
function [res]=NRMSE(I1,I2)
L2_error=I1-I2;
res=norm(L2_error(:),2)./norm(I1(:),2);