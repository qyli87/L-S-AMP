%--------------------------------------------------------
% Yangqing Li 
% 31 Oct 2015
%--------------------------------------------------------

function rankGuess = Rank_contraction(X, tau, Rank_contract)

% Data check
X(isnan(X)) = 0;

% Compute singular values of X
d = svd(X) + eps;

% Compute ratios
dtilde = d(1:(end-1)) ./ d(2:end)/tau;

% Compare ratio
[~,p] = find(dtilde>1);
rankGuess = sum(p);

% Use singular value profile of DCT
if ~isempty(Rank_contract)
    q = find(dtilde<1);
	rankGuess = q(1);
end

% Force at least rank = 1 
if rankGuess < 1,
    rankGuess = 1;
end