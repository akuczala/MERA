totSteps = 200;
maxLayers = 6;
energies = zeros(maxLayers,totSteps);
entropies = zeros(maxLayers,totSteps);
nLayers = [3 3 3 3 3 3];
for i = 1:maxLayers 
    out = MERA(3,nLayers(i),totSteps);
    energies(i,:) = real(out{1});
    entropies(i,:) = real(out{2});
end