function [kernel,rawkernel] = calibrateL(AtA, kSize, nCoil, coil, sampling)

if nargin < 6
	sampling = ones([kSize,nCoil]);
end

dummyK = zeros(kSize(1),kSize(2),nCoil);
dummyK((end+1)/2,(end+1)/2,coil) = 1;
idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);
Aty = AtA(:,idxY);
Aty = Aty(idxA);
AtA = AtA(idxA,:); 
AtA =  AtA(:,idxA);

kernel = sampling*0;
[U S V]=svd(AtA);
lambda = l_curve(U,diag(S),Aty,'tsvd');% gcv 
rawkernel= tsvd(U,diag(S),V,Aty,lambda);
kernel(idxA) = rawkernel; 
