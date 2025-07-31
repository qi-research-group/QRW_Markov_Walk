function K = buildK(a,b,c)
% K = buildK(a,b,c)
% m = number of states
% a = oﬀ diag left
% b = diag
% c = oﬀ diag right
m = size(a,1);
K = zeros(m);
K([1 2],1)= [ b(1) -b(1)];
K([m-1 m], m) = [-b(m) b(m)];
for k = 2:(m-1)
v = [(k-1) k (k+1)] ;
K(v,k) = [a(k) b(k) c(k)];
end