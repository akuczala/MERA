function newh = updateH(w,u,wdag,udag,h) %preserves hermicity of H
%Links inspection 8/11/2014
%9/9/14 link sequence optimized with netcon
tensors = {w,w,u,h,udag,wdag,wdag};
%center h
legLinks = {[1 2 3 -3],[4 5 6 -4],[7 8 3 4],[9 10 7 8],[9 10 11 12],[1 2 11 -1],[12 5 6 -2]};
%sequence = [3 4 7 8 9 10 11 1 2 12 5 6];
sequence = [1 2 5 6 9 10 7 8 4 12 3 11];
newh = ncon(tensors,legLinks,sequence);
%left h
legLinks = {[1 2 3 -3],[8 11 12 -4],[4 9 3 8],[6 5 2 4],[5 9 7 10],[1 6 7 -1],[10  11 12 -2]};
%sequence = [5 4 9 6 7 1 2 3 11 12 8 10];
sequence = [11 12 7 5 6 4 9 1 2 3 8 10];
newh = newh + ncon(tensors,legLinks,sequence);
%right h
legLinks = {[3 2 1 -4],[12 11 8 -3],[9 4 8 3],[5 6 4 2],[9 5 10 7],[7 6 1 -2],[12  11 10 -1]};
%sequence = [5 4 9 6 7 1 2 3 11 12 8 10];
sequence = [11 12 7 5 6 4 9 1 2 3 8 10];
newh = newh + ncon(tensors,legLinks,sequence);
%newh = newh/3; %this is wrong/shenglong
end