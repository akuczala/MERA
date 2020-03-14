%Computes scale invariant rho by finding eigenoperator of descending
%operator
function newrho = SIRho(w,v,u,wconj,vconj,uconj,rho) %generates hermitian rho
d = length(w);
%8/13/14 MUST transpose for rho to be eigenoperator
%rng(0); %fix seed
rng('shuffle') %shuffle seed
opts.v0 = reshape(rho,d^8,1); %initial Lanczos vector
opts.isreal = 0;
opts.issym = 0;
rng('shuffle'); %shuffle seed
%DFunc = @(x) reshape(transpose(reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^2,d^2)),d^4,1);
DFunc = @(x) reshape(descend(w,v,u,wconj,vconj,uconj,reshape(x,d,d,d,d,d,d,d,d)),d^8,1);

[rho,evals] = eigs(DFunc,d^8,1,'lm',opts); %find dominant eigenvector of D
%[V,evals] = eig(Dmat);
%rho = V(:,1);
%max(max(Dmat*rho - rho))
%disp('Rho eigenvalue');
%disp(diag(evals))
rho = reshape(rho,d^4,d^4);

rho = rho/trace(rho); %normalize, divide out phase

%DEBUG: round rho (keeping, since it reduced randomness)
rho = roundTensor(rho,14);

newrho = reshape(rho,d,d,d,d,d,d,d,d); %reshape eigvector into rho
toc;
%disp(sort(diag(evals)));
%DEBUGGING
% ev = diag(evals);
% for j=1:1
%     disp('num');
%     disp(j);
%     disp('eigv');
%     disp(ev(j));
%     eop = reshape(V(:,j),d,d,d,d);
%     Deop = ncon({D,eop},{[1 2 3 4 -1 -2 -3 -4],[1 2 3 4]},[1 2 3 4]);
%     %eig eqn
%     max(max(max(max(Deop - ev(j)*eop))))
%     %trace
%     ncon({eop},{[1 2 1 2]},[1 2])
%     %trace after
%     ncon({Deop},{[1 2 1 2]},[1 2])
%     %hermiticity
%     disp(checkHermiticity(eop));
%     %and after
%     disp(checkHermiticity(Deop));
% end
end
