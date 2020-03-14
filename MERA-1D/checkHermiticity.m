function diff = checkHermiticity(A)
d= length(A);
B = reshape(A,d^2,d^2);
diff = B - transpose(conj(B));
diff = max(diff(:));
end