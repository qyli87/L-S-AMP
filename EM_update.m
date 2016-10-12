function Updates = EM_update(y, A, State, Params)
%% Check for errors and unpack inputs

T = length(y) - 1;
N = size(State.Mu_x,1);

% Unpack State structure
MU_STATE = State.Mu_x;
V_STATE = State.V_x;
PI_OUT = State.Pi_out;
    
% Unpack Params structure
lambda = Params.lambda;

%% Start estimating the updated parameters

% Initial values for updated variables
Updates.lambda = Params.lambda;
Updates.sig2e = Params.sig2e;

% Updated noise variance
sig2e_sum = 0; M = NaN*ones(1,length(y));
for t = 1:T+1
    sig2e_sum = sig2e_sum + norm(y{t} - A{t}*MU_STATE(:,t))^2 + ...
        sum(V_STATE(:,t));
    M(t) = length(y{t});
end
Updates.sig2e = (sig2e_sum)/sum(M);

% Updated spasity level
MU_S = lambda.*prod(PI_OUT, 2) ./ (prod((1 - PI_OUT), 2).*(1 - lambda) ...
    + prod(PI_OUT, 2).*lambda);
Updates.lambda = sum(MU_S) / N;
