d = 4; d1 = d; d2 = d;

randmat = rand(d1^5,d2)+rand(d1^5,d2)*i;
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d2,:));
w = reshape(w,d1,d1,d1,d1,d1,d2);
wconj = conj(w);

randmat = rand(d1^4,d2^2)+rand(d1^4,d2^2)*i;
[U,S,V] = svd(randmat);
v = ctranspose(U(1:d2^2,:));
v = reshape(v,d1,d1,d1,d1,d2,d2);
vconj = conj(v);

randmat = rand(d1^4,d2^4)+rand(d1^4,d2^4)*i;
[U,S,V] = svd(randmat);
u = ctranspose(U(1:d2^4,:)); % d^3 x chi1
u = reshape(u,d1,d1,d1,d1,d2,d2,d2,d2);
uconj = conj(u);

randmat = rand(d1^4,d2^4)+rand(d1^4,d2^4)*i;
rho = randmat;
%rho = randmat*randmat'; rho = rho/trace(rho);
rho = reshape(rho,d1,d1,d1,d1,d2,d2,d2,d2);

tensors = {wconj,wconj,wconj,wconj,vconj,vconj,vconj,vconj,uconj,u,v,v,v,v,w,w,w,w,rho};
legLinks = {[1 2 3 4 5 61 ],[6 7 8 9 10 62 ],[11 12 13 14 15 63 ],[16 17 18 19 20 64 ],[21 22 23 24 4 13 ],[25 26 27 28 9 18 ],[29 30 31 32 2 6 ],[33 34 35 36 12 16 ],[-1 -2 -3 -4 31 24 36 27 ],[-5 -6 -7 -8 45 46 47 48 ],[21 22 23 46 49 50 ],[25 26 48 28 51 52 ],[29 30 45 32 53 54 ],[33 34 35 47 55 56 ],[1 53 3 49 5 57 ],[54 7 8 51 10 58 ],[11 55 50 14 15 59 ],[56 17 52 19 20 60 ],[57 58 59 60 61 62 63 64 ]};
seq = [21 22 23 1 3 5 29 30 32 2 53 7 8 10 25 26 28 9 51 11 14 15 17 19 20 33 34 35 16 56 12 55 59 60 63 64 18 52 58 62 6 54 57 61 4 13 49 50 45 46 47 48 24 27 31 36];
tic;
newrho = ncon(tensors,legLinks,seq);
toc;
disp(trace(reshape(rho,d1^4,d2^4)));
disp(trace(reshape(newrho,d1^4,d2^4)));
newrhomat = reshape(newrho,d^4,d^4);
diff = newrhomat - transpose(conj(newrhomat));
diff = max(diff(:))

legLinks = {[1 2 3 4 5 61 ],[6 7 8 9 10 62 ],[11 12 13 14 15 63 ],[16 17 18 19 20 64 ],[-1 -3 23 24 4 13 ],[25 26 27 28 9 18 ],[29 30 31 32 2 6 ],[33 34 35 36 12 16 ],[-2 38 -4 40 31 24 36 27 ],[-6 38 -8 40 45 46 47 48 ],[-5 -7 23 46 49 50 ],[25 26 48 28 51 52 ],[29 30 45 32 53 54 ],[33 34 35 47 55 56 ],[1 53 3 49 5 57 ],[54 7 8 51 10 58 ],[11 55 50 14 15 59 ],[56 17 52 19 20 60 ],[57 58 59 60 61 62 63 64 ]};
seq =  [1 3 5 29 30 32 2 53 7 8 10 25 26 28 9 51 11 14 15 17 19 20 33 34 35 16 56 12 55 59 60 63 64 18 52 58 62 6 54 57 61 45 47 48 27 31 36 38 40 46 49 50 4 13 23 24];
tic;
newrho = ncon(tensors,legLinks,seq);
toc;
disp(trace(reshape(rho,d1^4,d2^4)));
disp(trace(reshape(newrho,d1^4,d2^4)));
newrhomat = reshape(newrho,d^4,d^4);
diff = newrhomat - transpose(conj(newrhomat));
diff = max(diff(:))