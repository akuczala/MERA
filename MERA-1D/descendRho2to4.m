function out = descendRho2to4(w,wdag,u,udag,rho2)
    tensors = {udag,wdag,wdag,rho2,w,w,u};
    links = {[1 2 -9 -10],[-7 -8 1 3],[2 -11 -12 4],[5 6 3 4], ...
        [-1 -2 7 5],[8 -5 -6 6],[ -3 -4 7 8]};
    out = ncon(tensors,links,[1 2 3 4 5 6 7 8]);
end