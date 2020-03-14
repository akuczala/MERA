g = [1 1.1 1.15 1.2 1.25 1.3 1.35 1.4];
n = length(g);
mags = zeros(n,300);
parfor j = 1:n
    run{j} = MERA(g(j),2,300,1);
    mags(j,:) = real(run{j}{4});
end
plot(transpose(mags));