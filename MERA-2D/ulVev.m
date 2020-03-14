function out = ulVev(w,v,u,h,r)

tensors = {wconj,wconj,wconj,wconj,vconj,vconj,vconj,vconj,uconj,h,u,v,v,v,v,w,w,w,w,r};
links = {[1 2 3 4 5 61],[6 7 8 9 10 62],[11 12 13 14 15 63],[16 17 18 19 20 64],[21 22 23 24 4 13],[25 26 27 28 9 18],[29 30 31 32 2 6],[33 34 35 36 12 16],[37 38 39 40 31 24 36 27],[5 29 21 37 41 42 43 44],[44 38 39 40 45 46 47 48],[43 22 23 46 49 50],[25 26 48 28 51 52],[42 30 45 32 53 54],[33 34 35 47 55 56],[1 53 3 49 41 57],[54 7 8 51 10 58],[11 55 50 14 15 59],[56 17 52 19 20 60],[57 58 59 60 61 62 63 64]};

out = ncon(tensors,links,seq);
end