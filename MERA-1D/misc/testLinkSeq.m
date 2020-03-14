d = 4;
rng(1);
[w,u,wdag,udag] = initTensors(d,d,false);
r = rand(d*d,d*d); rho = r*ctranspose(r);
rho = reshape(rho,d,d,d,d);
tensors = {udag,wdag,wdag,rho,w,w,u};
legLinks = {[-3 -4 3 4],[1 2 3 7],[4 5 6 8],[9 10 7 8],[1 2 11 9],[12 5 6 10],[-1 -2 11 12]};
sequence = [7 3 1 2 9 11 12 10 5 6 8 4];
%sequence = netcon(legLinks,0,2,1,1);
sequence = [1 2 5 6 8 10 7 9 11 12 3 4];
%sequence = randperm(12);

tic
ncon(tensors,legLinks,sequence);
time = toc

s = repmat('%d ',1,length(sequence));
s(end)=[]; %Remove trailing comma

disp(sprintf(['Answer: [' s ']'], sequence))