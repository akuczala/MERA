function rho1 = singleSiteRho(rho)
    rho1 = 0.5*ncon({rho},{[-1 1 -2 1]},[1]);
    rho1 = rho1 + 0.5*ncon({rho},{[1 -1 1 -2]},[1]);
end