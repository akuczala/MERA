%analytical ground state of transverse field ising model
function out = tfimgs(g)
    out = 2/pi*(1+g)*ellipticE(4*g/(1+g)^2);
end