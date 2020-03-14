%legLinks = {[-4 6 2 3],[1 -3 2 7],[3 4 5 8],[9 10 7 8],[1 -1 11 9],[12 4 5 10],[-2 6 11 12]};

seq = netcon(legLinks,0,2,1,1);
s = repmat('%d ',1,length(seq));
s(end)=[]; %Remove trailing comma

disp(sprintf(['[' s ']'], seq))