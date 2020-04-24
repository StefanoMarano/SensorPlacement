function plotArray(pos)

N=size(pos,1);
pos=pos(:,1:2);

plot(pos(:,1),pos(:,2),'ro','markersize',12,'MarkerFaceColor',[1 0 0])

%axis('square')
title(sprintf('Array layout N=%d',N))
xlabel('x [m]')
ylabel('y [m]')
for nn=1:N
%     t_h = text(pos(nn,1),pos(nn,2),sprintf('   %d',nn))
%     set(t_h,'fontsize',14)
end
axis('equal')
axis('square')
grid on


[minD maxD] = maxminD(pos);
ArrayCenter=0.5*[min(pos(:,1))+max(pos(:,1)) min(pos(:,2))+max(pos(:,2))];

xlim(ArrayCenter(1) + maxD*0.6*[-1 +1])
ylim(ArrayCenter(2) + maxD*0.6*[-1 +1])