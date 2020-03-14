%for j=1:length(links)
seq = {};
for j=1:length(links)
    legLinks = links{j};
    seq{j} = netcon(legLinks,0,2,1,1);
    disp(j);
end
s = repmat('%d ',1,length(seq{j}));
s(end)=[]; %Remove trailing comma
disp(sprintf(['[' s '], ...\n'], seq{:}))
%sprintf(s, seq{:})

