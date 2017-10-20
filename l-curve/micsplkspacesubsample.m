function [subkspace,ACS,acqACS,acq,mask]  = micsplkspacesubsample(kspace,nACS,s)
%kspace is the channel data for a slice
%nACS is the number of ACS line
%s is the subsampling factor

[Npe,Nfe,C]=size(kspace); 

% The index for ACS lines in the PE direction
ACS =transpose( [floor(Npe/2 - nACS/2 + 1):floor(Npe / 2 + nACS/2)]); 

% The index for subsampled PE lines in the row direction
acq = 1:s:Npe;
if ismember(ACS(1) - 1, acq)
    ACS = union(ACS(1) - 1, ACS);
end    
if ismember(ACS(end) + 1, acq)
    ACS = union(ACS, ACS(end) + 1);
end    

% The index of neighbouring points within ACS lines.
acqACS = intersect(acq, ACS);

% The index of subsampled PE and ACS lines
acq = union(acq, ACS);

% Subsampling k-space data along row-direction. ACS lines are included in the subsampled data. The size of PE lines is incremented by one in order to reconstruct data in the final block
subkspace = zeros(Npe, Nfe, C);
mask=zeros(size(subkspace));
for index = 1: C
    subkspace(acq, :, index) = kspace(acq, :, index);
    mask(acq, :, index)=1;
end
return;
