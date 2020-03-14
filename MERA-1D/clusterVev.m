function out = clusterVev(w,u,wdag,udag,rho,rho3s)
Z = [1 0; 0 -1]; X = [0 1; 1 0]; id = [1 0; 0 1];
XL = kron(id,X); ZL = kron(Z,Z);
clusterOp = ncon({ZL,XL,ZL},{[-1 -4],[-2 -5],[-3 -6]},[]);
tensors = {clusterOp,w,w,u,udag,wdag,wdag,rho};
legLinks = {[5 7 8 4 10 9 ],[3 4 11 1 ],[12 14 15 2 ],[10 9 11 12 ],[7 8 6 13 ],[3 5 6 16 ],[13 14 15 17 ],[1 2 16 17 ]};
seq = [14 15 2 17 7 8 5 6 13 16 9 10 12 1 3 4 11];
out = 1:3;

a = ncon(tensors,legLinks,seq);
out(1) = a(1,1);
legLinks = {[15 13 12 16 14 11 ],[3 17 5 1 ],[6 11 8 2 ],[16 14 5 6 ],[15 13 4 7 ],[3 17 4 9 ],[7 12 8 10 ],[1 2 9 10 ]};


seq = [13 15 7 12 14 16 6 8 11 2 10 4 9 1 3 5 17];
b =  ncon(tensors,legLinks,seq);
out(2) = b(1,1);
tensors = {w,w,w,u,u,udag,udag,wdag,wdag,wdag,clusterOp,rho3s};
legLinks = {[1 2 3 22 ],[6 10 12 21 ],[16 19 20 23 ],[4 7 3 6 ],[13 17 12 16 ],[4 8 5 9 ],[14 17 15 18 ],[1 2 5 26 ],[9 11 15 25 ],[18 19 20 24 ],[8 11 14 7 10 13 ],[22 21 23 26 25 24 ]};
seq = [15 19 20 23 24 18 25 16 17 12 21 10 11 13 14 8 9 4 6 7 5 26 1 2 3 22];
out(3) = ncon(tensors,legLinks,seq);
end