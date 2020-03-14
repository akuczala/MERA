function eigstatetest(w,u,wdag,udag,rho,rho3)
d = length(rho);
Z = [1 0; 0 -1]; X = [0 1; 1 0]; id = [1 0; 0 1];
XL = kron(id,X); ZL = kron(Z,Z);
op = ncon({ZL,XL,ZL},{[-1 -4],[-2 -5],[-3 -6]},[]);
out = {};
tensors = {op,w,w,u,udag,wdag,wdag,rho};
leglinks = {[5 7 8 4 10 9 ],[3 4 11 1 ],[12 14 15 2 ],[10 9 11 12 ],[7 8 6 13 ],[3 5 6 -1 ],[13 14 15 -2 ],[1 2 -3 -4 ]};
seq = [14 15 7 8 5 6 9 10 3 4 11 12 13 1 2];
out{1} = ncon(tensors,leglinks,seq);



tensors = {op,w,w,u,udag,wdag,wdag,rho};
leglinks = {[10 12 13 9 11 14 ],[3 15 5 1 ],[6 14 8 2 ],[9 11 5 6 ],[10 12 4 7 ],[3 15 4 -1 ],[7 13 8 -2 ],[1 2 -3 -4 ]};
seq = [3 15 10 12 7 13 9 11 6 8 14 4 5 1 2];
out{2} = ncon(tensors,leglinks,seq);

tensors = {w,w,w,u,u,udag,udag,wdag,wdag,wdag,op,rho3};
leglinks = {[1 2 3 22 ],[6 10 12 21 ],[16 19 20 23 ],[4 7 3 6 ],[13 17 12 16 ],[4 8 5 9 ],[14 17 15 18 ],[1 2 5 -1 ],[9 11 15 -2 ],[18 19 20 -3 ],[8 11 14 7 10 13 ],[22 21 23 -4 -5 -6 ]};
seq = [1 2 19 20 6 15 11 14 8 9 13 17 4 7 10 12 16 18 3 5 21 22 23];
out{3} = ncon(tensors,leglinks,seq);

disp('avg/std after operating');
disp(mean(mean(mean(mean(rho./out{1})))));
disp(mean(mean(mean(mean(rho./out{2})))));
disp(mean(reshape(rho3./out{3},d^6,1)));
disp(std(reshape(rho./out{1},d^4,1)));
disp(std(reshape(rho./out{2},d^4,1)));
disp(std(reshape(rho3./out{3},d^6,1)));
disp('expectation');
disp(ncon({out{1}},{[1 2 1 2]},[1:2]));
disp(min(abs(reshape(rho./out{1},d^4,1))));
end

