function[PE_index]=twoacqneighbourindex(s,Npe)
%s is the sampling factor 
%Npe is the number of phase encoding lines
acq=1:s:Npe;
PE_index=zeros(Npe,2);
PE_index(acq(1:end-1),1)=acq(1:end-1);
PE_index(acq(1:end-1),2)=acq(2:end);
for r=1:s-1,
    PE_index(r+(1:s:acq(end)-1),:)=PE_index(acq(1:end-1),:);
end
l=length((1+acq(end-1)):Npe);
v=PE_index(acq(end-1),:);
v=repmat(v,l,1);
PE_index((1+acq(end-1)):end,:)=v;