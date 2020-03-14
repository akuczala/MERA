out = {};
index = 1;
maxsteps = 3000;
maxlayers = 6;
maxbond = 4;
energies = zeros(maxlayers*(maxbond-1),maxsteps);
energygrid = zeros(maxlayers,maxbond-1);
for layers = 1:maxlayers
    for bond_dim = 2:2
        chi = ones(1,layers)*bond_dim
        out{index} = MERA(1,chi,maxsteps,1,false);
        energies(index,:) = [real(out{index}{1}),zeros(1,maxsteps-length(out{index}{1}))];
        energygrid(layers,bond_dim-1) = real(out{index}{1}(length(out{index}{1})))+tfimgs(1);
        index = index + 1;
    end
end