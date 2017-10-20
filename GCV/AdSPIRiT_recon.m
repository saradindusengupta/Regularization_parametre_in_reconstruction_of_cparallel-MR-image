function [K, E] = AdSPIRiT_recon(Kz, GOP, nIter, x0, im,bet)
%
%
% res = pocsSPIRIT(y, GOP, nIter, x0, T, show)
%
% Implementation of the Cartesian, POCS l1-SPIRiT reconstruction
%
% INPUTS:
%		data    -   Undersampled k-space data. Make sure that empty entries are zero
%                    or else they will not be filled.
%		GOP     -	the SPIRiT kernel operator obtained by calibration
%		nIter   -	Maximum number of iterations
%		x0      -   initial guess
%		T       -   wavlet threshodling parameter
% 		
% OUTPUTS:
%		res     - Full k-space matrix
%
% (c) Michael Lustig 2007
%
acq_ind= find(abs(Kz)>0);
% find the closest diadic size for the images
[sx,sy,nc] = size(Kz);
W          = Wavelet('Daubechies',6,4);
mask       = (Kz==0);
K          = x0;
bet=1;
for n=1:nIter
    Err = zeros(size(Kz));
    %------------ Apply (G-I)*x + x for reconstructing initially----
    Gk   = GOP*K;
	Kcon = (K + GOP*K ).*(mask);
    K    = (bet* Kcon) + Kz; 
    Err (acq_ind) = Kz(acq_ind)-Gk(acq_ind);
    cerr(n) = norm(Err(:));
	K    = K.*mask + Kz; % fix the data (POCS)
    Istk = ifft2c(K);
    I_recon=sos(Istk);
    E(n) = NRMSE(sos(im),I_recon);
end
% figure; plot(er); hold on; plot(ep,'r');
% figure; plot(cerr);


function x = softThresh(y,t)
% apply joint sparsity soft-thresholding 
absy = sqrt(sum(abs(y).^2,3));
unity = y./(repmat(absy,[1,1,size(y,3)])+eps);

res = absy-t;
res = (res + abs(res))/2;
x = unity.*repmat(res,[1,1,size(y,3)]);

function[X]=getDiadic(K)
% zpad to the closest diadic 
[sx,sy,nc] = size(K);
ssx = 2^ceil(log2(sx)); 
ssy = 2^ceil(log2(sy));
ss = max(ssx, ssy);
X= zpad(ifft2c(K),ss,ss,nc); 

function[T1,l1res,l1pert,delt] = Update_Thresh(wres,wpert,T)
l1res  = norm(sum(wres,3),1);
l1pert = norm(sum(wpert,3),1);
delt   = abs(l1res-l1pert);
Cres   = get_cov(wres);
Cpert  = get_cov(wpert);
res    = sum(sqrt(Cres(:)));
pert   = 1/(T^2)*sum(sqrt(Cpert(:)));
T1     = sqrt(res/(delt+pert));

%------------------- calling function ----------
function [res]=NRMSE(I1,I2)
L2_error=I1-I2;
res=norm(L2_error(:),2)./norm(I1(:),2);

function[eX]=expectatn(X)
[h,bin]=hist(X(:),1001);
h=h/sum(h);
eX=sum(bin.*h);
