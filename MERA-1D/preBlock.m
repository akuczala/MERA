function newh = preBlock(h,d)

u = eye(d^2,2*d);
u = reshape(u,d,d,2*d); udag = u;
%h on the left
legs = {[3 4 -3],[1 2 3 4],[1 2 -1],[-2 -4]};
newh = 0.5*ncon({u,h,udag,eye(2*d)},legs,[1 2 3 4]);
%h on the right
legs = {[-1 -3],[3 4 -4],[1 2 3 4],[1 2 -2]};
newh = newh + 0.5*ncon({eye(2*d),u,h,udag},legs,[1 2 3 4]);
%h in the middle
legs = {[1 5 -3],[6 4 -4],[2 3 5 6],[1 2 -1],[3 4 -2]};
newh = newh + ncon({u,u,h,udag,udag},legs,[2 1 5 6 4 3]);
end