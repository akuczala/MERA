function ee = vnEntropy(rho,d)

rhomat = reshape(rho,d,d);
evals = eig(rhomat);
ee = -sum(evals.*log2(evals));%uses log base 2
end