function [w,v,u,wconj,vconj,uconj] = initTensors(d1,d2)

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

end