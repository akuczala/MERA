function newrho = updateRho(w,u,wdag,udag,rho) %preserves hermiticity of rho

%fine grained rho
%tensors = {u,w,w,rho,wdag,wdag,udag};
%fixed tensor list 8/11/14
tensors = {udag,wdag,wdag,rho,w,w,u};
%9/9/14 replaced sequences with sequences given by netcon
%middle diagram
%legLinks = {[-1 -2 3 4],[1 2 3 7],[4 5 6 8],[7 8 9 10],[1 2 11 9],[12 5 6 10],[-3 -4 11 12]};
%updated 8/12/14, swapped top/bottom of input/output
legLinks = {[-3 -4 3 4],[1 2 3 7],[4 5 6 8],[9 10 7 8],[1 2 11 9],[12 5 6 10],[-1 -2 11 12]};
%sequence = [7 3 1 2 9 11 12 10 5 6 8 4];
sequence = [1 2 5 6 8 10 7 9 11 12 3 4];
newrho = ncon(tensors,legLinks,sequence);
%left diagram
%updated as above 8/12/14
legLinks = {[-4 6 2 3],[1 -3 2 7],[3 4 5 8],[9 10 7 8],[1 -1 11 9],[12 4 5 10],[-2 6 11 12]};
%sequence = [7 2 1 9 11 6 12 10 4 5 3 8];
sequence = [4 5 8 10 12 9 11 1 7 2 3 6];
newrho = newrho + ncon(tensors,legLinks,sequence);
%right diagram
legLinks = {[6 -3 3 2],[5 4 3 7],[2 -4 1 8],[9 10 7 8],[5 4 12 9],[11 -2 1 10],[6 -1 12 11]};
%sequence = [8 2 1 10 11 6 12 9 4 5 3 7];
sequence = [4 5 7 9 12 10 11 1 8 2 3 6];
newrho = newrho + ncon(tensors,legLinks,sequence);
newrho = newrho/3;
%contract rho with D, see if eigenoperator (debug)
%tensors = {D,rho}; legLinks = {[-1 -2 -3 -4 1 2 3 4],[1 2 3 4]},sequence = [1 2 3 4];
%disp(ncon(tensors,legLinks,sequence)-rho)

end
