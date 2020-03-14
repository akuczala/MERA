function out = ellipticE(m)
    [K E] = ellipke(m);
    out = E;
end