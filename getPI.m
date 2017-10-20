function[PI, CHI,PE_index]=getPI(subkspace,s,b,indexing,acqACS)

%subkspace is undersampled k-space
%s is the sampling factor 
%b is the number  of acquired neighbours. b should be even
%indexing can be either linear or circular 1:linear  2 :circular
%acqACS is the k-space row indices following the pattern of acquired lines within ACS


[Npe,Nfe,C]=size(subkspace);%in subkspace , an additional unacquired line is included in the end
Npe=Npe-1;

%find coil indices for each coil
% COIL_index=findcoilindex(C,indexing); %size of COIL_index will be  C X C. Each row

% Acquired neighbour indices for each k-space row 
PE_index=twoacqneighbourindex(s,Npe);%size of PE_index will be Npe X 2

lim=( (b/2)-1 ) /2;
subkspacepad=padarray(subkspace,[lim,lim]);
subkspacepad=subkspacepad(lim+1:end-lim,:,:);
[P,Q,C]=size(subkspacepad);
PI=[];
for i=lim+1:Q-lim,
     for j=acqACS(1:end-1),
            temp=subkspacepad(PE_index(j,:),i-lim:i+lim,:);
            [x,y,z]=size(temp);
            PI = [PI; reshape(temp, 1, x*y*z)];
     end
end

CHI=PI'*PI;