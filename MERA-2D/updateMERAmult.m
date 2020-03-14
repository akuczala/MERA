function [neww,newv,newu,newwconj,newvconj,newuconj] = updateMERAmult(w,v,u,wconj,vconj,uconj,h,rho,steps)
    neww = w; newv = v; newu = u;
        newwconj = wconj; newvconj = vconj; newuconj = uconj;
    for j = 1:steps
           [neww,newv,newu,newwconj,newvconj,newuconj]...
               = updateMERA(neww,newv,newu,newwconj,newvconj,newuconj,h,rho);
    end
end