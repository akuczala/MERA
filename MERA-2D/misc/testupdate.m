d = 2; d1 = d; d2 = 3;
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

[neww,newv,newu,newwconj,newvconj,newuconj] = updateMERA(w,v,u,wconj,vconj,uconj,h,rho);

checkuUnitarity(u)
checkuUnitarity(newu)

checkvUnitarity(v)
checkvUnitarity(newv)

checkwUnitarity(w)
checkwUnitarity(neww)