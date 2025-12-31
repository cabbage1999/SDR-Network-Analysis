function mask = FDR_Threshold(alpha, thre, adj, pval, n)

%% Extract the upper triangle
triuIdx = find(triu(true(n),1));
pvec    = pval(triuIdx);

%% Benjamini–Hochberg FDR correction
[pvecSort, idxSort] = sort(pvec);
m = numel(pvec);
qvec = zeros(size(pvec));
qvec(idxSort) = pvecSort .* m ./ (1:m)'; 
qvec = max(qvec, pvec);

qmat = zeros(n);
qmat(triuIdx) = qvec;
qmat = qmat + qmat';

%% FDR + Threshold
mask = (qmat <= alpha) & (abs(adj) >= thre);