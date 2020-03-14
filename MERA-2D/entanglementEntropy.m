function ee = entanglementEntropy(rho,d)

[U,S,V] = svd(reshape(rho,d,d));
svalues = diag(S);
ee = -sum(svalues.*log2(svalues));%uses log base 2
end