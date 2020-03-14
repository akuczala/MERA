function out = MERA(g,chi,totSteps,stepSize,SIMERA,start)
%SEED FIXED
seed = 3;
rng(seed);
%rng('shuffle');
ansprec = 1e-9;
threesite = false;
%disable input checking in ncon
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;
%SIMERA toggles scale invariance (infinite MERA)
initId = false; %sets all initial tensors to identity
givenInitial = false;
if length(start) > 0
    givenInitial = true;
end
  
prec = 14; %precision in rounding, if used
initialBlocking = false;
rounding = false; %rounding has no significant effect at g=3
Eoffs = 4;
%totSteps = 5;
%layers = 3;
Elist = 1:totSteps; %list of energies at each iteration
Slist = 1:totSteps; %entropies
clist = 1:totSteps; %central charges
Mlist = 1:totSteps; %magnetizations 
stablist = zeros(totSteps,3);
%parameters
if initialBlocking
    d = 4;
else
    d = 2; %site dimension
end
d = 2; %site dimension
%chi = [4 4 4 4 4]; %transitional dimension
layers = length(chi);
chiFix = chi(length(chi)); %scale invariant dimension / top dimension
chiTop = 1; %number of eigs to keep at top
%initialize spin operators
Z = [1 0; 0 -1]; X = [0 1; 1 0]; id = [1 0; 0 1];
Zpot = [-1 0 0; 0 2 0; 0 0 -1]; Xpot = [0 1 0; 0 0 1; 1 0 0];
idpot = eye(3);
Y = [0 -i; i 0]; XL = kron(id,X); ZL = kron(Z,Z);
Magnet = ncon({Z,id},{[-1 -3],[-2 -4]},[]); %TFIM
MagnetX = ncon({X,id},{[-1 -3],[-2 -4]},[]); %TFIM
%Magnet = kron(kron(Z,Z),kron(id,id)); %1D gadget
%Magnet = kron(Zpot,idpot); %Potts
%Magnet = reshape(Magnet,d,d,d,d);

hsite = ncon({Z,Z},{[-1 -3],[-2 -4]},[]);
%hsite = hsite + ncon({Y,Y},{[-1 -3],[-2 -4]},[])/4;
%hsite = hsite + ncon({X,X},{[-1 -3],[-2 -4]},[])/4;

%transverse field ising model
hsite = hsite - g*MagnetX;
%hsite = Magnet;
%hsite = -2*kron(Z,Z); %ferromagnet
%hsite = -kron(X,id); %paramagnet
%hsite = kron(X,X)+g*kron(Y,Y); %XY
%hsite = -(1+g)*kron(Z,Z) - (1-g)*kron(Y,Y);
%hsite = -ncon({kron(X,X),kron(id,id)},{[-1 -3],[-2 -4]},[]);%gadget cluster H
%hsite = kron(Zpot,idpot) - kron(Xpot,Xpot') - kron(Xpot',Xpot); %Potts
%hbond = ncon({kron(Y,X),kron(Z,Z)},{[-1 -3],[-2 -4]},[]);

%hbond = hbond - ncon({kron(Z,Z),kron(Y,id)},{[-1 -3],[-2 -4]},[]);
%hsite = hsite + g*hbond;
hsite = reshape(hsite,d*d,d*d);
hsite = hsite - Eoffs*eye(d*d); %shift energy
hsite = reshape(hsite,d,d,d,d);
%block h0 into effective 4 dim sites
if initialBlocking
    h0 = preBlock(hsite,2);
    Magnet = preBlock(Magnet,2);
else
    h0 = hsite;
end
hAvg = eye(chiFix^2,chiFix^2);
hAvg = reshape(hAvg,chiFix,chiFix,chiFix,chiFix);
r = rand(chiFix^2,chiFix^2) + rand(chiFix^2,chiFix^2)*i;
rhoFix = r*r';
rhoFix = reshape(rhoFix,chiFix,chiFix,chiFix,chiFix);
r = rand(chiFix^3,chiFix^3) + rand(chiFix^3,chiFix^3)*i;
if threesite
    rhoFix3s = r*r';
    rhoFix3s = reshape(rhoFix3s,chiFix^3,chiFix^3);
end
laststab = 0;
%inital MERA matrices
%UV layer
%rng(seed);
[w0,u0,w0dag,u0dag] = initTensors(d,chi(1),initId);
if givenInitial
    w0 = start{1};u0 = start{2};
    w0dag = conj(w0); u0dag = conj(u0);
end
%transitional layers
for tau = 1:layers
   
    if tau < layers
        nextChi = chi(tau+1);
    else
        nextChi = chiFix;
    end
    %RESET SEED
    %rng(seed)
    [w{tau},u{tau},wdag{tau},udag{tau}] = initTensors(chi(tau),nextChi,initId);
    if givenInitial
        w{tau} = start{tau*2+1};u{tau} = start{tau*2 + 2};
        wdag{tau} = conj(w{tau}); udag{tau} = conj(u{tau});
    end
end


%fixed point layer
if SIMERA
    [wFix,uFix,wdagFix,udagFix] = initTensors(chiFix,chiFix,initId);
    if givenInitial
        wFix = start{2*layers+3};uFix = start{2*layers+4};
        wdagFix = conj(wFix); udagFix = conj(uFix);
    end
else
    top = eye(chiFix*chiFix,chiTop);
    if givenInitial
        top = start{2*layers+3};
    end
end
%initial tensors
%rho2 = eye(chi1^2); rho2 = reshape(rho2,chi1,chi1,chi1,chi1);
%rho1 = eye(chi1^2); rho1 = reshape(rho1,chi1,chi1,chi1,chi1);
%h1 = eye(chi1^2); h1 = reshape(h1,chi1,chi1,chi1,chi1);
%rho0 = eye(d^2); rho0 = reshape(rho0,d,d,d,d);
lastE = 0;
lastS = 0;
%hupd = h;
%rho = eye(chi^2); rho = reshape(rho,chi,chi,chi,chi);
%h2 = eye(chi1^2); h2 = reshape(h2,chi1,chi1,chi1,chi1);
%hTop = eye(chi^2); hTop = reshape(h2,chi,chi,chi,chi);
%last = {hupd,u,w,rho,h2,u1,w1,rho1,h1,u0,w0,rho0};
%new = last;
%diffvec = 1:length(last);
%rho = cell(layers);
%LOOP
T = 0;
for it= 1:totSteps
    tic
    disp(it);
    
    %Scale invariant/ top layer density matrix
    if SIMERA
        disp('SI rho');
        rhoFix = scaleInvariantRho(wFix,uFix,wdagFix,udagFix,rhoFix);
        if threesite
            rhoFix3s = scaleInvariantRho3s(wFix,uFix,wdagFix,udagFix,rhoFix3s);
        end
    else
        disp('layer');
        
        rhoFix = top*ctranspose(top);
        rhoFix = rhoFix/chiTop;
        rhoFix = rhoFix/trace(rhoFix);
        rhoFix = reshape(rhoFix,chiFix,chiFix,chiFix,chiFix);
    end
    disp('update rho');
    %descend MERA updating rho
    for tau = fliplr(1:layers)
        if tau == layers
            rho{tau} = updateRho(w{tau},u{tau},wdag{tau},udag{tau},rhoFix);
            if threesite
                rho3s{tau} = descend3c(rhoFix3s,w{tau},u{tau},wdag{tau},udag{tau});
            end
        else
            rho{tau} = updateRho(w{tau},u{tau},wdag{tau},udag{tau},rho{tau+1});
            if threesite
                rho3s{tau} = descend3c(rho3s{tau+1},w{tau},u{tau},wdag{tau},udag{tau});
            end
        end
    end
    rho0 = updateRho(w0,u0,w0dag,u0dag,rho{1});
    if threesite
        rho03s = descend3c(rho3s{1},w0,u0,w0dag,u0dag);
    end
    toc
    disp('update MERA');
    %ascend MERA updating u,w and h
    tic
    [w0,u0,w0dag,u0dag] = updateMERA(w0,u0,w0dag,u0dag,rho{1},h0,stepSize,T); %update layer 0 MERA
    
    for tau = 1:layers
        if(tau == 1)
            h{tau} = updateH(w0,u0,w0dag,u0dag,h0); %update layer 1 hamiltonian
        else
            h{tau} = updateH(w{tau-1},u{tau-1},wdag{tau-1},udag{tau-1},h{tau-1});
        end
        if(tau == layers)
            [w{tau},u{tau},wdag{tau},udag{tau}] = updateMERA(w{tau},u{tau},wdag{tau},udag{tau},rhoFix,h{tau},stepSize,T);
        else
            [w{tau},u{tau},wdag{tau},udag{tau}] = updateMERA(w{tau},u{tau},wdag{tau},udag{tau},rho{tau+1},h{tau},stepSize,T);
        end
    end
    toc
    %get average h for SI layers / update top h
    hTop = updateH(w{layers},u{layers},wdag{layers},udag{layers},h{layers});
    tic
    if SIMERA
        %approximate layer average
        %hAvg = hTop +updateH(wFix,uFix,wdagFix,udagFix,hAvg)/3;
        %hTop = hAvg;
        %sum a few levels
        %hTop2 = updateH(wFix,uFix,wdagFix,udagFix,hTop);
        %hTop3 = updateH(wFix,uFix,wdagFix,udagFix,hTop2);
        %hTop4 = updateH(wFix,uFix,wdagFix,udagFix,hTop3);
        %hTop = hTop + hTop2/3 + hTop3/9 + hTop4/27;
        %disp('update scale invariant layer');
        %hupd = updateH(wFix,uFix,wdagFix,udagFix,hFix);
        %hFix = hTop + hupd/3; %update scale invariant h
        %h = h2; %DEBUG
        %hTop2 = updateH(wFix,uFix,wdagFix,udagFix,hTop);
        %hTop3 = updateH(wFix,uFix,wdagFix,udagFix,hTop2);
        %hTop4 = updateH(wFix,uFix,wdagFix,udagFix,hTop3);
        %hAvg = hTop + updateH(wFix,uFix,wdagFix,udagFix,hAvg);
        %%%%%%%%%%%%% hupd/3 term does nothing for real valued u,w
        %hTop = hTop/3 + hTop2/9 + hTop3/27 + hTop4/81; %next-order terms in sum
        %next order terms decrease E - does not consistently improve accuracy
        %of matlab.matmatlab.matenergyreshape(
    end
    
    %update scale invariant tensors / top tensor
    if SIMERA
        [wFix,uFix,wdagFix,udagFix] = updateMERA(wFix,uFix,wdagFix,udagFix,rhoFix,hTop,stepSize,T);
        %disp(wFix)
    else
        %update top tensor, fill with lowest energy eigenvalues of H_(T-1)
        %rng(0);
        %[V,D] = eigs(reshape(h{layers},chiFix^2,chiFix^2),chiTop);
        [V,D] = eig(reshape(hTop,chiFix^2,chiFix^2));
        evals = real(diag(D)); %take real part so sorting works
        [esort,sorti] = sort(evals);
        %minIndex = find(evals == min(evals));
        minIndex = sorti(1);
        minIndex2 = sorti(2);
        %disp(evals);
        %disp(minIndex);
        %disp(minIndex2);
        %disp(evals);
        %disp(min(evals));
        %disp(evals(minIndex)-min(evals));
        if chiTop == 2
            V = horzcat(V(:,minIndex),V(:,minIndex2));
        end
        if chiTop == 1
            V = V(:,minIndex);
        end
        if chiTop == 3
            V = horzcat(V(:,minIndex),V(:,minIndex2), ...
                V(:,sorti(3)));
        end
        if chiTop == 4
            V = horzcat(V(:,minIndex),V(:,minIndex2), ...
                V(:,sorti(3)),V(:,sorti(4)));
        end
        %disp(evals);
        %V = V(:,minIndex2);
        %rng('shuffle');
        %[V,D] = eig(reshape(h{layers},chiFix^2,chiFix^2));
        %V = transpose(V(1:chiFix,:));
        top = V;
        %top = roundTensor(top,14);
    end
    toc
    
    %round tensors to precision
    %seems not to help get better consistency of gsE between runs
    if rounding
        u0 = roundTensor(u0,prec);u0dag = conj(u0);
        w0 = roundTensor(w0,prec);w0dag = conj(w0);
        rho0 = roundTensor(rho0,prec);
        for tau = 1:layers
            u{tau} = roundTensor(u{tau},prec); udag{tau} = conj(u{tau});
            w{tau} = roundTensor(w{tau},prec); wdag{tau} = conj(w{tau});
            rho{tau} = roundTensor(rho{tau},prec);
        end
    end
    
    %energy expectation value
    E = vev2(h0,rho0) + Eoffs;
    %magnetisation
    M = vev2(Magnet,rho0);
    disp('E');
    disp(E);
    %show energies at each layer
    %Evec = 1:layers;
    %for tau = 1:layers
    %    Evec(tau) = log(real(vev2(h{tau},rho{tau})/E))/log(3);
    %end
    %disp(Evec);
    
    %disp('Efix')
    %disp(log(real(vev2(hAvg,rhoFix)/E))/log(3));
    disp('de');
    dE = E-lastE;
    disp(dE);
    disp('accuracy');
    disp(log(real(E+tfimgs(1)))/log(10));
    %disp(log(real(E+tfimgs(1))));
    disp('M');
    disp(M);
    
    %debug
    %disp(eig(reshape(rho0,d*d,d*d)));
    %disp(eig(reshape(rho{1},4*4,4*4)));
    %disp(eig(reshape(rho{2},4*4,4*4)));
    %disp(eig(reshape(rho{3},4*4,4*4)));
    %disp(eig(reshape(rhoFix,chiFix^2,chiFix^2)))
    %plot1 = stem(eig(reshape(rhoFix,chiFix^2,chiFix^2)));
    if(SIMERA)
    S2 = entanglementEntropy(rhoFix,chiFix*chiFix);
    S1 = entanglementEntropy(singleSiteRho(rhoFix),chiFix)
    else
    S2 = entanglementEntropy(rho0,d*d);
    S1 = entanglementEntropy(singleSiteRho(rho0),d)
    end
    %S0 = entanglementEntropy(rho0,d*d);
    cc = 3*(S2-S1);
    disp('c');
    disp(cc);
    %plotcmatrix(...
    %    [reshape(rho0,d^2,d^2),...
    %    reshape(rho{1},chi(1)^2,chi(1)^2),...
    %    reshape(rho{2},chi(2)^2,chi(2)^2),...
    %    reshape(rhoFix,chiFix^2,chiFix^2)]')
    %image(255*real([sort(eig(reshape(rho0,d^2,d^2))),...
    %    sort(eig(reshape(rho{1},chi(1)^2,chi(1)^2))),...
    %    sort(eig(reshape(rho{2},chi(2)^2,chi(2)^2))),...
    %    sort(eig(reshape(rhoFix,chiFix^2,chiFix^2)))]'))
    
    %disp(reshape(rho{1},chi(1)^2,chi(1)^2))
    
    %image(real(reshape(rhoFix,chiFix^2,chiFix^2)))
    S = S2;
    dS = S-lastS;
    %disp('dS');
    %disp(dS);
    disp('T');
    disp(T);
    %cc = 0;
    %show largest changes in tensors
    %new = {hupd,u,w,rho,h2,u1,w1,rho1,h1,u0,w0,rho0};
    %for dispi = 1:length(last)
    %    diff = abs(new{dispi}-last{dispi});
    %    diffvec(dispi) = max(diff(:));
    %end
    %disp(diffvec);
    %last = new;
    lastE = E;
    lastS = S;
    Elist(it) = E;
    Slist(it) = S;
    clist(it) = cc;
    Mlist(it) = M;
    semilogy(real(Elist+tfimgs(1)));
    %attempt #2 to measure cluster operator
    %eigstatetest(w0,u0,w0dag,u0dag,rho{1},rho3s{1});
    
    %stablist(it,:) = clusterVev(w0,u0,w0dag,u0dag,rho{1},rho3s{1});
    %dstab = (stablist(it,1)-laststab);
    %laststab = stablist(it,1);
    %disp(transpose(stablist(it,:)));
    %normtest(w0,u0,w0dag,u0dag,rho{1},rho3s{1});
    if(abs(dE) < ansprec) %stop when E reaches precision
        Elist = Elist(1:it);
        Slist = Slist(1:it);
        clist = clist(1:it);
        Mlist = Mlist(1:it);
        
    end
    T = T/1.5;
end
if SIMERA
    out = {Elist,Slist,clist,Mlist,stablist,{rho0,rho,rhoFix}, ...
        {w0,w,wFix},{u0,u,uFix},{h0,h,hTop}};
    %disp(singleSiteScaling(wFix,wdagFix,5));
else
    out = {Elist,Slist,clist,Mlist,{rho0,rho,rhoFix}, ...
        {w0,w},{u0,u},{h0,h,hTop},top};
end

end

