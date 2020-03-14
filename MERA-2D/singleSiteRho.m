function rho1 = singleSiteRho(rho)
    rho1 = ncon({rho},{[-1 1 2 3 -2 1 2 3]},[1:3]);
    rho1 = rho1 + ncon({rho},{[1 -1 2 3 1 -2 2 3]},[1:3]);
    rho1 = rho1 + ncon({rho},{[1 2 -1 3 1 2 -2 3]},[1:3]);
    rho1 = rho1 + ncon({rho},{[1 2 3 -1 1 2 3 -2]},[1:3]);
    rho1 = rho1/4;
end