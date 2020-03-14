d = 4;
n = 50;
r = rand(d*d,d*d)+rand(d*d,d*d)*i;
rho = r*ctranspose(r);
rho = rho/trace(rho);
rho = reshape(rho,d,d,d,d);
r = rand(d*d,d*d)+rand(d*d,d*d)*i;
h = r*ctranspose(r);
h = h/trace(h);
h = -h;
disp('h eigs');
eig(h);
h = reshape(h,d,d,d,d);
[w,u,wdag,udag] = initTensors(d,d,false);
energy = 1:n;
entropy = 1:n;
for j = 1:n
[w,u,wdag,udag] = updateMERAonce(w,u,wdag,udag,rho,h);
rho0 = updateRho(w,u,wdag,udag,rho);
energy(j) = vev2(h,rho0);
entropy(j) = entanglementEntropy(rho0,d*d);
disp(reshape(rho0,d*d,d*d));
end
plot(real(entropy))