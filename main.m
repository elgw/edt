M = zeros(1024,1024);
M(3,1) = 1;
M(12,1) = 1;


tic
for kk = 1:1000
    D = bwdist(M);
end
toc
