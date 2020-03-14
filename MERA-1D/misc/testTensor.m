d = 2;
u = eye(d*d,d*d); u = reshape(u,d,d,d,d);
udag = conj(u);
randmat = rand(d*d*d,d)+ rand(d*d*d,d)*i;
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d,:));
w = reshape(w,d,d,d,d);
wdag = conj(w);
rho = rand(d,d,d,d);
%w = rand();
%rho = rand(_);

%Descending operator
tensors = {u,w,w,wdag,wdag,udag};
%middle diagram
legLinks = {[-1 -2 3 4],[1 2 3 -5],[4 5 6 -6],[1 2 7 -7],[8 5 6 -8],[-3 -4 7 8]};
sequence = [3 1 2 7 8 5 6 4];
D = ncon(tensors,legLinks,sequence);
%left diagram
legLinks = {[-2 6 2 3],[1 -1 2 -5],[3 4 5 -6],[1 -3 7 -7],[8 4 5 -8],[-4 6 7 8]};
sequence = [2 1 7 6 8 4 5 3];
D = D + ncon(tensors,legLinks,sequence);
%right diagram
legLinks = {[6 -1 3 2],[5 4 3 -5],[2 -2 1 -6],[5 4 8 -7],[7 -4 1 -8],[6 -3 8 7]};
sequence = [2 1 7 6 8 4 5 3];
D = D + ncon(tensors,legLinks,sequence);
D = D/3;
%right diagram
Dmat = reshape(D,d^4,d^4); %make matrix out of D tensor
[rho,evals] = eigs(Dmat,1); %find dominant eigenvector of D
rho = reshape(rho,d,d,d,d); %reshape eigvector into rho

%fine grained rho
tensors = {u,w,w,rho,wdag,wdag,udag};
%middle diagram
legLinks = {[-1 -2 3 4],[1 2 3 7],[4 5 6 8],[7 8 9 10],[1 2 11 9],[12 5 6 10],[-3 -4 11 12]};
sequence = [7 3 1 2 9 11 12 10 5 6 8 4];
newrho = ncon(tensors,legLinks,sequence);
%left diagram
legLinks = {[-2 6 2 3],[1 -1 2 7],[3 4 5 8],[7 8 9 10],[1 -3 11 9],[12 4 5 10],[-4 6 11 12]};
sequence = [7 2 1 9 11 6 12 10 4 5 3 8];
newrho = newrho + ncon(tensors,legLinks,sequence);
%right diagram
legLinks = {[6 -1 3 2],[5 4 3 7],[2 -2 1 8],[7 8 9 10],[5 4 12 9],[11 -4 1 10],[6 -3 12 11]};
sequence = [8 2 1 10 11 6 12 9 4 5 3 7];
newrho = newrho + ncon(tensors,legLinks,sequence);
newrho = newrho/3;
%contract rho with D, see if eigenoperator
tensors = {D,rho}; legLinks = {[-1 -2 -3 -4 1 2 3 4],[1 2 3 4]},sequence = [1 2 3 4];
disp(ncon(tensors,legLinks,sequence)-rho)




