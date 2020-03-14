function newrho = descend2C(w,u,wdag,udag,rho) %preserves hermiticity of rho

%fine grained rho
%tensors = {u,w,w,rho,wdag,wdag,udag};
%fixed tensor list 8/11/14
tensors = {udag,wdag,wdag,rho,w,w,u};
%middle diagram
%legLinks = {[-1 -2 3 4],[1 2 3 7],[4 5 6 8],[7 8 9 10],[1 2 11 9],[12 5 6 10],[-3 -4 11 12]};
%updated 8/12/14, swapped top/bottom of input/output
legLinks = {[-3 -4 3 4],[1 2 3 7],[4 5 6 8],[9 10 7 8],[1 2 11 9],[12 5 6 10],[-1 -2 11 12]};
sequence = [7 3 1 2 9 11 12 10 5 6 8 4];
newrho = ncon(tensors,legLinks,sequence);

end
