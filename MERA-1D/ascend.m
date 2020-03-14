function A = ascend(w,u,wdag,udag)
d = length(u);
%Links inspection 8/11/2014
tensors = {w,w,u,udag,wdag,wdag};
%center h
legLinks = {[1 2 3 -7],[4 5 6 -8],[-3 -4 3 4],[-1 -2 11 12],[1 2 11 -5],[12 5 6 -6]};
sequence = [3 4 11 1 2 12 5 6];
A = ncon(tensors,legLinks,sequence)/3;
%left h
legLinks = {[1 -3 3 -7],[8 11 12 -8],[-4 9 3 8],[-2 9 7 10],[1 -1 7 -5],[10  11 12 -6]};
sequence = [9 7 1 3 11 12 8 10];
A = A + ncon(tensors,legLinks,sequence)/3;
%right h
legLinks = {[3 -4 1 -8],[12 11 8 -7],[9 -3 8 3],[9 -1 10 7],[7 -2 1 -6],[12  11 10 -5]};
sequence = [9 7 1 3 11 12 8 10];
A = A + ncon(tensors,legLinks,sequence)/3;
out = A;
Amat = reshape(A,d^4,d^4);
evals = eig(Amat);
disp(sort(evals));
end