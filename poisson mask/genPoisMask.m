function [mask] = genPoisMask(siz,samper,dist)
sper     = samper;
attempt  = 0;
maxatmpt = 10;
mskstk   = [];
interf   = [];
vec      =@(x) x(:);
for k=1:1,
    while(1)
        mask    = zeros(siz);
        rows    = round(siz(1)/2)-round(siz(1)*dist/2):round(siz(1)/2)+round(siz(1)*dist/2);
        clms    = round(siz(2)/2)-round(siz(2)*dist/2):round(siz(2)/2)+round(siz(2)*dist/2);
        spacing = 1.2;
        nPts    = prod(siz)*samper-numel(rows)*numel(clms);
        pts     = poissonDisc(siz,spacing,round(nPts));
        for i=1:size(pts,1), mask(round(pts(i,1)),round(pts(i,2)))=1; end
        mask(rows,clms) = 1;
        fprintf('sampling percentage = %f\n',numel(find(mask))/prod(siz));
        if abs(numel(find(mask))/prod(siz)-sper)>.01, samper=2*sper-numel(find(mask))/prod(siz);
        else break; end
        attempt = attempt+1;
        if attempt==maxatmpt, disp('maximum attempts tried..! quitting'); return; end
    end
    mskstk(:,:,k) = mask;
    interf(k)     = max(vec(abs(ifft2(mask))));
    if k==15,
        ind  = find(interf==min(interf(:)));
        mask = mskstk(:,:,ind(1)); 
    end
end
% if abs(numel(find(mask))/prod(siz)-samper)>.01, disp('nsamples is more/less by .01%'); end
% fprintf('sampling percentage = %f\n',numel(find(mask))/prod(siz));