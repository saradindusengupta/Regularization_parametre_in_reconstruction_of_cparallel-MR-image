function[GOP]=l1_Calib(DATA,CalibSize,kSize,lambda)

% ------------Perform Calibration------------
[Nfe, Npe,C] = size(DATA); 

% to get the calibration data 
kCalib = crop(DATA,[CalibSize,C]);
kernel = zeros([kSize,C,C]);

% Formation of Calibration Matrix
[AtA] = dat2AtA(kCalib, kSize);

%Calculation of coil coefficients

for n=1:C
	kernel(:,:,:,n) = calibrategcv(AtA,kSize,C,n);
    [kernel,rawkernel] = calibrategcv(AtA, kSize, nCoil, coil, sampling)
end

% Formation of SPIRiT operator
GOP = SPIRiT(kernel, 'fft',[Nfe,Npe]);