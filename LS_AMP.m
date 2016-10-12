%--------------------------------------------------------
% Yangqing Li 
% 31 Oct 2015
%--------------------------------------------------------

function x_hat = LS_AMP(y, A, rank_contract)
%% Initialize variables
% Initialize dimensionality of data
T = length(y);
N = size(A{1}, 2);

% Initialize variables
tol = 5e-3;             % Tolerance for early termination
resid_energy_tol = 0.975;  % Tolerance of the difference of residual energy
rank_tau = 1.5;         % Minimum ratio of singular value gap
last_iter = 0;          % Flag for early termination
eps = 1e-6;             % Small positive scalar for msg approx.
lambda = 0.1*ones(N,1); % Initialize to an arbitrary value, EM will learn it
avg_resid_energy = 0;   % The column-averaged residual energy

% Default options
inter_iters = 25;       % Max # of inter-phase iterations
min_iters = 10;          % Min # of inter-phase iterations
inner_iters = 20;       % Max # of intra-phase AMP iterations
energy_win_l = 6;       % The observed energy window length
energy_win_l = min(energy_win_l, min_iters); % Cannot larger than min_iters
verbose = true;         % Run silently or not

%Check for provided state
if nargin < 3
    rank_contract = [];
end

%% Run the M-GAMP
% Message matrix declarations and initializations (all are N-by-T dim)
PI_IN = lambda*ones(1,T);	% Matrix of messages from s to f(t)
PI_OUT = NaN*ones(N,T);     % Matrix of messages from f(t) to s
XI_IN = zeros(N,T);         % Matrix of means from theta(t) to f(t)
PSI_IN = ones(N,T);         % Matrix of vars from theta(t) to f(t)
PHI_OUT = NaN*ones(N,T);

% Create space to store the state of the GAMP variables between iterations
C_STATE = 100*ones(1,T);	% Initialize to large c values
Z_State = y;                % Initialize residual to the measurements
MU_STATE = zeros(N,T);

% Declare constants
Eq_Params.eps = eps;        % Gaussian "squelch" parameter on p(x_n)
Eq_Params.sig2e = 0.001;    % additive Gaussian noise variance
Eq_Options.iter = inner_iters;  % Number of GAMP iterations

% Initialize some variables
x_hat = cell(1,T);          % Placeholder for MMSE of x from M-GAMP
X_hat = zeros(N, T);        % Placeholder for MMSE of X from M-GAMP
StateOutAMP = cell(1,T);    % Placeholder for M-GAMP output
x_old = cell(1,T);          % Hold previous estimates for early termination
V_STATE = NaN*ones(N,T);    % Placeholder for x variances

old_resid_energy = inf*ones(1,inter_iters-min_iters+energy_win_l);
resid_energy_ratios = zeros(1,energy_win_l);

start_time = tic;           % Start stopwatch

for k = 1:inter_iters       % Inter-phase iterations
%% M-GAMP iteration loop
    stop_MGAMP = false;
    % Recoveries of x(t) for every colume of X "in parallel"
    while ~stop_MGAMP
        for t = 1:T           % Current column-step index
            Eq_Params.pi = PI_IN(:,t); Eq_Params.xi = XI_IN(:,t); 
            Eq_Params.psi = PSI_IN(:,t);
            EqStateIn.c = C_STATE(t); EqStateIn.z = Z_State{t};
            EqStateIn.mu = MU_STATE(:,t);
            
            [x_hat{t}, StateOutAMP{t}] = GAMP_frame(y{t}, A{t}, ...
                Eq_Params, Eq_Options, EqStateIn);
        end

        % Compute the column-averaged residual energy
        pre_energy = avg_resid_energy;
        avg_resid_energy = 0;
        for t = 1:T
            avg_resid_energy = avg_resid_energy + norm(y{t} - A{t}*x_hat{t})^2 / T;
        end

        % In rare cases, rerun MGAMP when it diverges in small size problem
        if avg_resid_energy < 100
            stop_MGAMP = true;
            if k > min_iters && avg_resid_energy/pre_energy > 2
                last_iter = 1;
                x_hat = x_old;
            end
        else
            if k < min_iters	% Reinitialization
                if (length(lambda) == 1)
                lambda = 0.1*ones(N,1);
                end
                PI_IN = lambda*ones(1,T);
                PI_OUT = NaN*ones(N,T);
                XI_IN = randn(N,T);
                PSI_IN = ones(N,T);
                C_STATE = 100*ones(1,T);
                Z_State = y;
                MU_STATE = zeros(N,T);
            else
                stop_MGAMP = true;
                last_iter = 1;
                x_hat = x_old;
            end
        end
    end
    
    % Compute the column-averaged residual energy ratios
    if k > min_iters-energy_win_l
        old_resid_energy(k-min_iters+energy_win_l) = avg_resid_energy;
    end
    if k > min_iters-energy_win_l+1
        resid_energy_ratios(1:energy_win_l-1) = resid_energy_ratios(2:energy_win_l);
        resid_energy_ratios(energy_win_l) = old_resid_energy(k+energy_win_l-min_iters)/old_resid_energy(k+energy_win_l-1-min_iters);
    end
    
    % Display execution info
    if verbose
        fprintf('Completed %d inter-phase iterations\n', k);
        fprintf('Total elapsed time: %f s\n', toc(start_time));
        fprintf('Column-averaged residual energy: %f\n', avg_resid_energy);
        fprintf('---------------------------------------------------------\n');
    end
    
    for t = 1:T
        PHI_OUT(:,t) = StateOutAMP{t}.phi;
        C_STATE(t) = StateOutAMP{t}.c;
        PI_OUT(:,t) = StateOutAMP{t}.pi; 
        Z_State{t} = Z_State{t} + StateOutAMP{t}.z;
        MU_STATE(:,t) = StateOutAMP{t}.mu;
        V_STATE(:,t) = StateOutAMP{t}.v;
        X_hat(:,t) = x_hat{t};
    end
    
    %Check for rank update for each inter-phase interation
    R = Rank_contraction(X_hat, rank_tau, rank_contract);

%% Run the L&SPD
    [XI_IN, PSI_IN] = LSPD(X_hat, PHI_OUT, C_STATE, R);
    
%% Back to the M-GAMP
    for t = 1:T
        if t < T
            % Update all the PI_INs for the next round
            inc_ind = [1:t-1, t+1:T];     % Timesteps in the PI_OUT prod.
            PI_IN(:,t) = lambda.*prod(PI_OUT(:,inc_ind), 2) ./ ...
                (lambda.*prod(PI_OUT(:,inc_ind), 2) + ...
                (1 - lambda).*prod(1 - PI_OUT(:,inc_ind), 2));
            for i_PI = 1:N
                if isnan(PI_IN(i_PI,t)),
                    PI_IN(i_PI,t) = 0;
                end
            end
        else
            % Update PI_IN(:,T) for the next round
            inc_ind = 1:T-1;
            PI_IN(:,t) = lambda.*prod(PI_OUT(:,inc_ind), 2) ./ ...
                (lambda.*prod(PI_OUT(:,inc_ind), 2) + ...
                (1 - lambda).*prod(1 - PI_OUT(:,inc_ind), 2));
            for i_PI = 1:N
                if isnan(PI_IN(i_PI,t)),
                    PI_IN(i_PI,t) = 0;
                end
            end
        end
    end
    
%% Run the EM-learning
    % Pass current state of factor graph to the parameter update function
    State = struct('Mu_x', []);     % Clear contents of structure
    State.Mu_x = MU_STATE;
    State.V_x = V_STATE;
    State.Pi_out = PI_OUT;
        
    % Also pass the current values of the parameters
    UpdParams.lambda = lambda;
    UpdParams.sig2e = Eq_Params.sig2e;
        
    % Now compute the updates
    Updates = EM_update(y, A, State, UpdParams);
        
    % Store the updated hyperparameters back
    lambda = Updates.lambda;
    Eq_Params.sig2e = .01*Eq_Params.sig2e + .99*Updates.sig2e;
    
    
%% Check termination
    % Check for early termination this round and terminate the inter-phase
    % interations
    if last_iter, break; end
    
    % Check for early inter-phase termination for next round
    if k > min_iters && norm([x_hat{:}] - [x_old{:}], 'fro')/norm([x_hat{:}], 'fro') < tol
        last_iter = 1;      % Set the flag for last iteration
    else
        x_old = x_hat;
        resid_energy_sum = sum(resid_energy_ratios);
    end
    
    % Check for early inter-phase termination for next round
    if k > min_iters && resid_energy_sum > resid_energy_tol*energy_win_l
        last_iter = 1;      % Set the flag for last iteration
    end
end