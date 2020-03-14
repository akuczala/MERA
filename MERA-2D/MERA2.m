function out = MERA2(g,layers,totSteps,stepSize,SIMERA)
%SEED FIXED
seed = 3;
rng(seed);
%rng('shuffle');
ansprec = 1e-6;
%disable input checking in ncon
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;
%SIMERA toggles scale invariance (infinite MERA)
initId = false; %sets all initial tensors to identity
Eoffs = 4;
Elist = 1:totSteps;
Slist = 1:totSteps;
clist = 1:totSteps;
Mlist = 1:totSteps;
stablist = zeros(totSteps,9);
d = 2 %site dimension
chi = [2]; %transitional dimension
chiFix = 2; %scale invariant dimension / top dimension
chiTop = 1; %number of eigs to keep at top

%initialize
Z = [1 0; 0 -1]; X = [0 1; 1 0]; id = [1 0; 0 1];
Magnet = kron(kron(Z,id),kron(id,id)); %TFIM
%Magnet = reshape(Magnet,d,d,d,d);

hsite = -ncon({Z,Z,id,id},{[-1 -5],[-2 -6],[-3 -7],[-4,-8]});
hsite = hsite - ncon({Z,id,Z,id},{[-1 -5],[-2 -6],[-3 -7],[-4,-8]});
hsite = hsite + g*ncon({X,id,id,id},{[-1 -5],[-2 -6],[-3 -7],[-4,-8]});
%hsite = -kron(kron(X,id),kron(id,id)); %para
%hsite = -kron(kron(Z,Z),kron(id,id))...
%    - kron(kron(Z,id),kron(Z,id)); %ferro
%hsite = ncon({X,Z,Z,X},{[-1 -5],[-2 -6],[-3 -7],[-4 -8]},[]); %Wen TC
%hsite = hsite + ...
%    g*ncon({Z,id,id,id},{[-1 -5],[-2 -6],[-3 -7],[-4 -8]},[]); %para term
hsite = reshape(hsite,d^4,d^4);

%hsite = -kron(kron(Z,X),kron(Z,X)); 
hsite = hsite - Eoffs*eye(d^4); %shift energy
h0 = reshape(hsite,d,d,d,d,d,d,d,d);

r = rand(chiFix^4,chiFix^4) + rand(chiFix^4,chiFix^4)*i;
rhoFix = r*r';
rhoFix = reshape(rhoFix,chiFix,chiFix,chiFix,chiFix, ...
    chiFix,chiFix,chiFix,chiFix);

%inital MERA matrices
%UV layer
%rng(seed);
[w0,v0,u0,w0conj,v0conj,u0conj] = initTensors(d,chi(1));
%transitional layers
for tau = 1:layers
   
    if tau < layers
        nextChi = chi(tau+1);
    else
        nextChi = chiFix;
    end
    [w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau}] = initTensors(chi(tau),nextChi);

end

%fixed point layer
if SIMERA
    [wFix,vFix,uFix,wconjFix,vconjFix,uconjFix] = initTensors(chiFix,chiFix);
else
    top = eye(chiFix^4,chiTop);
end

lastE = 0;
lastS = 0;

%LOOP
for it= 1:totSteps
    tic
    disp(it);
    
    %Scale invariant/ top layer density matrix
    if SIMERA
        disp('SI rho');
        rhoFix = SIRho(wFix,vFix,uFix,wconjFix,vconjFix,uconjFix,rhoFix);
    else
        disp('layer');
        rhoFix = top*ctranspose(top);
        rhoFix = rhoFix/chiTop;
        rhoFix = rhoFix/trace(rhoFix);
        rhoFix = reshape(rhoFix,chiFix,chiFix,chiFix,chiFix, ...
            chiFix,chiFix,chiFix,chiFix);
    end
    disp('update rho');
    %descend MERA updating rho
    for tau = fliplr(1:layers)
        if tau == layers
            rho{tau} = descend(w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau},rhoFix);
        else
            rho{tau} = descend(w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau},rho{tau+1});
        end
    end
    rho0 = descend(w0,v0,u0,w0conj,v0conj,u0conj,rho{1});
    toc
    disp('update MERA');
    %ascend MERA updating u,w and h
    tic
    [w0,v0,u0,w0conj,v0conj,u0conj] = updateMERAmult(w0,v0,u0,w0conj,v0conj,u0conj,h0,rho{1},stepSize); %update layer 0 MERA
    
    for tau = 1:layers
        if(tau == 1)
            h{tau} = ascend(w0,v0,u0,w0conj,v0conj,u0conj,h0); %update layer 1 hamiltonian
        else
            h{tau} = ascend(w{tau-1},v{tau-1},u{tau-1},wconj{tau-1},vconj{tau-1},uconj{tau-1},h{tau-1});
        end
        if(tau == layers)
            [w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau}] ...
                = updateMERAmult(w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau},h{tau},rhoFix,stepSize);
        else
            [w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau}]...
                = updateMERAmult(w{tau},v{tau},u{tau},wconj{tau},vconj{tau},uconj{tau},h{tau},rho{tau+1},stepSize);
        end
    end
    
    %get average h for SI layers / update top h
    hTop = ascend(w{layers},v{layers},u{layers},wconj{layers},vconj{layers},uconj{layers},h{layers});
    toc
    tic
    disp('top/SI');
    %update scale invariant tensors / top tensor
    if SIMERA
        [wFix,vFix,uFix,wconjFix,vconjFix,uconjFix] = updateMERAmult(wFix,vFix,uFix,wconjFix,vconjFix,uconjFix,hTop,rhoFix,stepSize);
    else
        %update top tensor, fill with lowest energy eigenvalues of H_(T-1)
       
        [V,D] = eigs(reshape(hTop,chiFix^4,chiFix^4),1);
        evals = real(diag(D)); %take real part so sorting works
        evals = sort(evals);
        V = V(:,1);
        top = V;
        
    end
    toc
    
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
    
    disp('Efix')
    %disp(log(real(vev2(hAvg,rhoFix)/E))/log(3));
    disp('de');
    dE = E-lastE;
    
    disp(dE);
    disp('M');
    disp(M);
  
    if(SIMERA)
    S2 = entanglementEntropy(rhoFix,chiFix^4);
    S1 = entanglementEntropy(singleSiteRho(rhoFix),chiFix)
    else
    S2 = entanglementEntropy(rho0,d^4);
    S1 = entanglementEntropy(singleSiteRho(rho0),d)
    end
    
    %S0 = entanglementEntropy(rho0,d*d);
    cc = 3*(S2-S1);
    %disp('c');
    %disp(cc);
    S = S2;
    dS = S-lastS;
    disp('dS');
    disp(dS);
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
    %stablist(it,:) = eachVev(w0,v0,u0,h0,rho{1})+4;
    Mlist(it) = M;
    disp(transpose(stablist(it)));
    if(abs(dE) < ansprec)
        Elist = Elist(1:it);
        Slist = Elist(1:it);
        clist = Elist(1:it);
        Mlist = Elist(1:it);
        break
    end
end
if SIMERA
    out = {Elist,Slist,stablist,Mlist,{rho0,rho,rhoFix}, ...
        {w0,w,wFix},{v0,v,vFix},{u0,u,uFix},{h0,h,hTop}};
    %disp(singleSiteScaling(wFix,wconjFix,5));
else
    out = {Elist,Slist,clist,Mlist,{rho0,rho,rhoFix}, ...
        {w0,w},{v0,v},{u0,u},{h0,h,hTop},top};
end

end


