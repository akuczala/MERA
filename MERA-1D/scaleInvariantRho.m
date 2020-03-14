%Computes scale invariant rho by finding eigenoperator of descending
%operator
function newrho = scaleInvariantRho(w,u,wdag,udag,rho) %generates hermitian rho
d = length(w);
%8/13/14 MUST transpose for rho to be eigenoperator
%rng(0); %fix seed
rng('shuffle') %shuffle seed
opts.v0 = reshape(rho,d^4,1); %initial Lanczos vector
opts.isreal = 0;
opts.issym = 0;
rng('shuffle'); %shuffle seed
%DFunc = @(x) reshape(transpose(reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^2,d^2)),d^4,1);
DFunc = @(x) reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^4,1);

[newrho,evals] = eigs(DFunc,d^4,1,'lm',opts); %find dominant eigenvector of D
%[V,evals] = eig(Dmat);
%newrho = V(:,1);
%max(max(Dmat*rho - rho))
%disp('Rho eigenvalue');
%disp(diag(evals))
newrho = reshape(newrho,d^2,d^2);

newrho = newrho/trace(newrho); %normalize, divide out phase

%DEBUG: round rho (keeping, since it reduced randomness)
newrho = roundTensor(newrho,14);

newrho = reshape(newrho,d,d,d,d); %reshape eigvector into rho
toc;
%DEBUGGING
%ev = diag(evals);
% for j=1:4
%     disp('num');
%     disp(j);
%     disp('eigv');
%     disp(ev(j));
%     eop = reshape(V(:,j),d,d,d,d);
%     %eop = newrho;
%     %Deop = ncon({D,eop},{[1 2 3 4 -1 -2 -3 -4],[1 2 3 4]},[1 2 3 4]);
%     Deop = updateRho(w,u,wdag,udag,eop);
%     %eig eqn
%     disp('eig eqn');
%     max(max(max(max(Deop - ev(j)*eop))))
%     %trace
%     disp('trace in/out')
%     ncon({eop},{[1 2 1 2]},[1 2])
%     %trace after
%     ncon({Deop},{[1 2 1 2]},[1 2])
%     %hermiticity
%     disp('herm in/out')
%     disp(checkHermiticity(eop));
%     %and after
%     disp(checkHermiticity(Deop));
% end
end
