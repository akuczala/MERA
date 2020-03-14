function out = normtest(w,u,wdag,udag,rho,rho3)

out = [1:2];

tensors = {w,w,u,udag,wdag,wdag,rho};
leglinks = {[3 14 5 1 ],[6 8 9 2 ],[12 13 5 6 ],[12 13 4 7 ],[3 14 4 10 ],[7 8 9 11 ],[1 2 10 11 ]};
seq = [8 9 12 13 6 7 2 11 4 10 1 3 5 14];
nrm= ncon(tensors,leglinks,seq);
out(1) = nrm(1,1);

tensors = {w,w,w,u,u,udag,udag,wdag,wdag,wdag,rho3};
leglinks = {[1 2 3 16 ],[6 23 8 15 ],[10 13 14 17 ],[4 21 3 6 ],[22 11 8 10 ],[4 21 5 7 ],[22 11 9 12 ],[1 2 5 20 ],[7 23 9 19 ],[12 13 14 18 ],[16 15 17 20 19 18 ]};
seq = [4 21 13 14 11 22 10 12 17 18 9 19 8 15 23 6 7 5 20 1 2 3 16];
nrm = ncon(tensors,leglinks,seq);
out(2) = nrm(1,1);
disp(transpose(out));

end

