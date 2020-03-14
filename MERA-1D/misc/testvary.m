
d = 16;
d2 = 32;
rng(1);
r = (rand(d2,d2) + i*rand(d2,d2));
h1 = r*ctranspose(r);
h1 = -h1/trace(h1);
r = (rand(d,d) + i*rand(d,d));
h2 = r*ctranspose(r);
h2 = -h2/trace(h2);
h2 = eye(d);
r = (rand(d,d) + i*rand(d,d));
h3 = r*ctranspose(r);
h3 = -h3/trace(h3);
%h3 = eye(d);
r = (rand(d,d) + i*rand(d,d));
h4 = r*ctranspose(r);
h4 = -h4/trace(h4);
h4 = eye(d);
n = 200;
nruns = 5;
mins = zeros(nruns,n);
tr = zeros(nruns,n);
for run = 1:nruns
    maxchange = 1:n;
    rng('shuffle');
    r = rand(d,d);
    [A,S,V] = svd(r);
    r = rand(d2,d2);
    [B,S,V] = svd(r);
    B = B(:,1:d);
    for j = 1:n
        %disp(j);
        f = h2*ctranspose(A)*ctranspose(B)*h1*B;
        %oldA = A;
        newA = vary(A,f);
        
        f = A*h2*ctranspose(A)*ctranspose(B)*h1;
        newB = vary(B,f);
        %tr(run,j) = trace(A*ctranspose(A));
        %maxchange(j) = max(max(abs(A-oldA)));
        A = newA;
        B = newB;
        mins(run,j) = trace(ctranspose(A)*ctranspose(B)*h1*B*A*h2);
    end
    disp(run);
end
%semilogy(transpose(mins)+50.816);
%semilogy(real(trace(B*C)-transpose(mins)));
plot(real(transpose(mins)))
disp(mins(:,n));