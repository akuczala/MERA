d = 2; d1 = d; d2 = 3;

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

randmat = rand(d1^4,d1^4)+rand(d1^4,d1^4)*i;
h = randmat;
h = randmat*randmat';
h = reshape(h,d1,d1,d1,d1,d1,d1,d1,d1);

tic;
newh = ascend(w,v,u,wconj,vconj,uconj,h);
toc;
disp(trace(reshape(h,d1^4,d1^4)));
disp(trace(reshape(newh,d2^4,d2^4)));
disp(trace(reshape(newh,d2^4,d2^4))-trace(reshape(h,d1^4,d1^4)));
newhmat = reshape(newh,d2^4,d2^4);
checkHerm2(h)
checkHerm2(newh)

%eig(newrhomat)
