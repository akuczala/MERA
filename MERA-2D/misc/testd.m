d = 2; d1 = d; d2 = 2;
rng('shuffle');
randmat = rand(d1^5,d2)+rand(d1^5,d2)*i;
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d2,:));
w = reshape(w,d1,d1,d1,d1,d1,d2);
wconj = conj(w);

randmat = rand(d1^4,d1^2)+rand(d1^4,d1^2)*i;
[U,S,V] = svd(randmat);
v = ctranspose(U(1:d1^2,:));
v = reshape(v,d1,d1,d1,d1,d1,d1);
vconj = conj(v);

randmat = rand(d1^4,d1^4)+rand(d1^4,d1^4)*i;
[U,S,V] = svd(randmat);
u = ctranspose(U(1:d1^4,:));
u = reshape(u,d1,d1,d1,d1,d1,d1,d1,d1);
uconj = conj(u);

randmat = rand(d2^4,d2^4)+rand(d2^4,d2^4)*i;
rho = randmat;
%rho = randmat*randmat'; rho = rho/trace(rho);
rho = reshape(rho,d2,d2,d2,d2,d2,d2,d2,d2);

randmat = rand(d1^4,d1^4)+rand(d1^4,d1^4)*i;
h = randmat;
%h = randmat*randmat';
h = reshape(h,d1,d1,d1,d1,d1,d1,d1,d1);

tic;
newrho = descend(w,v,u,wconj,vconj,uconj,rho);
toc;
disp(trace(reshape(rho,d2^4,d2^4)));
disp(trace(reshape(newrho,d1^4,d1^4)));
disp(trace(reshape(newrho,d1^4,d1^4))-trace(reshape(rho,d2^4,d2^4)));
newrhomat = reshape(newrho,d1^4,d1^4);
checkHerm2(rho);
checkHerm2(newrho);
diff = newrhomat - transpose(conj(newrhomat));
diff = max(diff(:))
diff = newrhomat - transpose(conj(newrhomat));
diff = max(diff(:))

%eig(newrhomat)
%check eigenvalues of D and A
newrho = SIRho(w,v,u,wconj,vconj,uconj,rho);

rng('shuffle') %shuffle seed
opts.v0 = reshape(h,d^8,1); %initial Lanczos vector
opts.isreal = 0;
opts.issym = 0;
rng('shuffle'); %shuffle seed
AFunc = @(x) reshape(ascend(w,v,u,wconj,vconj,uconj,reshape(x,d,d,d,d,d,d,d,d)),d^8,1);

[V,evals] = eigs(AFunc,d^8,5,'lm',opts);
disp(sort(diag(evals)));
