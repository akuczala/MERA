%fixed order of indices 8/19/14
function out = vev2(op,rho)
    out = ncon({rho,op},{[1 2 3 4 5 6 7 8],[5 6 7 8 1 2 3 4]},[1:8]);
end