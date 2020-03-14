rng('shuffle');
d=2;
[w,u,wdag,udag] = initTensors(d,d,false);
rng('shuffle') %shuffle seed
r = rand(d^2,d^2)+rand(d^2,d^2)*i;
rho = r*r';
rho =r;
r = rand(d^2,d^2)+rand(d^2,d^2)*i;
h = r*r';
h =r;
opts.v0 = reshape(rho,d^4,1); %initial Lanczos vector
opts.isreal = 0;
opts.issym = 0;
rng('shuffle'); %shuffle seed
%DFunc = @(x) reshape(transpose(reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^2,d^2)),d^4,1);
DFunc = @(x) reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^4,1);
[V,evals] = eigs(DFunc,d^4,10,'lm',opts); %find dominant eigenvector of D
disp(sort(diag(evals)));

opts.v0 = reshape(h,d^4,1); %initial Lanczos vector
opts.isreal = 0;
opts.issym = 0;
rng('shuffle'); %shuffle seed
%DFunc = @(x) reshape(transpose(reshape(updateRho(w,u,wdag,udag,reshape(x,d,d,d,d)),d^2,d^2)),d^4,1);
AFunc = @(x) reshape(updateH(w,u,wdag,udag,reshape(x,d,d,d,d)),d^4,1);
[V,evalsa] = eigs(AFunc,d^4,10,'lm',opts); %find dominant eigenvector of D
disp(sort(diag(evalsa)));
