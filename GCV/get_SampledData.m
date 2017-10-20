function[DATA, CalibSize, scale_fctr]=get_SampledData(DATA,mask)
% 
C=size(DATA,3);
%-----------------------UnderSampling--------------------------
% multiply with sampling matrix to get undersampled data

DATA = DATA.*repmat(mask,[1,1,C]); 

% to estimate the largest region in the center of k-space 
% that is fully sampled in order to use for calibraion

[CalibSize, dcomp] = getCalibSize(mask);

% scale the data such that the zero-filled density compensated      
% k-space norm is 1. This is useful in order to use similar         
% regularization penalty values for different problems.             

DATAcomp = DATA.*repmat(dcomp,[1,1,C]);
scale_fctr = norm(DATAcomp(:))/sqrt(C)/20;
DATA = DATA/scale_fctr;
DATAcomp = DATAcomp/scale_fctr;