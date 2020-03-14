totsteps = 300;
g = [ferrog parag];
energies = 1:length(g);
entropies = 1:length(g);
mags = 1:length(g);
for j = 1:length(g)
    conv{j} = MERA(g(j),4,totsteps,1);
    energies(j) = conv{j}{1}(totsteps);
    entropies(j) = conv{j}{2}(totsteps);
    mags(j) = conv{j}{4}(totsteps);
end
actual = arrayfun(@(x) tfimgs(x),g)
semilogy(g,energies+actual);