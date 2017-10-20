function[GOP]=tsvdLCalib(DATA,CalibSize,kSize)

% ------------Perform Calibration------------
[Nfe, Npe,C] = size(DATA); 

% to get the calibration data 
kCalib = crop(DATA,[CalibSize,C]);
kernel = zeros([kSize,C,C]);

% Formation of Calibration Matrix
[AtA] = dat2AtA(kCalib, kSize);

%Calculation of coil coefficients

for n=1:C
    kernel(:,:,:,n) = calibrateL(AtA, kSize,C,n);

end

% Formation of SPIRiT operator
GOP = SPIRiT(kernel, 'fft',[Nfe,Npe]);