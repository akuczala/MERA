function newh = ascend2C(w,u,wdag,udag,h) %preserves hermicity of H
%Links inspection 8/11/2014
tensors = {w,w,u,h,udag,wdag,wdag};
%center h
legLinks = {[1 2 3 -3],[4 5 6 -4],[7 8 3 4],[9 10 7 8],[9 10 11 12],[1 2 11 -1],[12 5 6 -2]};
sequence = [3 4 7 8 9 10 11 1 2 12 5 6];
newh = ncon(tensors,legLinks,sequence);
%newh = newh/3; %this is wrong
end