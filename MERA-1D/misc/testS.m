r = rand(16,16) + i*rand(16,16);
rho = r*ctranspose(r);
disp(trace(rho));
rho = reshape(rho,4,4,4,4);
rhod = updateRho(w1,u1,conj(w1),conj(u1),rho);
disp(trace(reshape(rhod,16,16)));