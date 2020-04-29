function [minD maxD] = maxminD(pos)
% finds minimum interstation distance minD and maximum array aperture maxD
maxD=-Inf;
minD=+Inf;
N=size(pos,1);
for ii=1:N
  for jj=ii+1:N
    tmp=norm(pos(ii,:)-pos(jj,:));
    if(tmp>maxD)
      maxD=tmp;
    end
    if(tmp<minD)
      minD=tmp;
    end
  end
end