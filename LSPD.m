%--------------------------------------------------------
% Yangqing Li 
% 31 Oct 2015
%--------------------------------------------------------

function [theta_Uout, theta_Vout] = LSPD(Y, PHI_OUT, C_STATE, RankN)
%% Problem Setup

% Initialize dimensionality of data
[M, L] = size(Y);
N = RankN;

% Create Options object
BiGAMPopt = BiGAMPOpt();

BiGAMPopt.M = M;
BiGAMPopt.L = L;
BiGAMPopt.N = N;

% Options to control EM-learning
EMopt.learn_var = true; % learn the variance of the X entries
EMopt.learn_mean = true; % learn the mean of the X entries
EMopt.sig_dim = 'joint'; % learn a single variances for X entries ('joint')
% or a different variance per (row) 'row' or (column) 'col'

EMopt.maxTol = 1e-4; % largest allowed tolerance for a single iteration

%% Initial Setup

% Initial values
nuw = norm(Y,'fro')^2/(M*L*101); % initial noise variance estimate
meanX = zeros(N,L); % initial mean estimate of X elements
nuX = (norm(Y,'fro')^2/M/L - nuw)/N; % initial variance estimate of active X elements

% Initialize xhat and Ahat
% xhat = sqrt(nuX)*randn(N,L);
xhat = zeros(N,L); % seems to perform better
Ahat = randn(M,N);

% Set initial step size small
BiGAMPopt.step = BiGAMPopt.stepMin;

% Set them
BiGAMPopt.xhat0 = xhat;
BiGAMPopt.Ahat0 = Ahat;
BiGAMPopt.Avar0 = ones(M,N);
BiGAMPopt.xvar0 = nuX*ones(N,L);

state = [];

%% Main loop
% Initialize loop
stop = 0;

% The < condition is so that we can do an iterative EM update within L&SPD
while stop < 2

    % Prior on A
    gA = AwgnEstimIn(zeros(M,N), ones(M,N));
%     PI = repmat(0.5,M,N);
%     gA = SparseScaEstim(gA,PI);

    % Prior on X
    gX = AwgnEstimIn(meanX, nuX);

    % Output log likelihood
    gOut = AwgnEstimIn(PHI_OUT, repmat(C_STATE,M,1));

    % Run BiG-AMP
    [xhat, xvar, Ahat, Avar, EMstate, state] = BiGAMP(gX, gA, gOut, BiGAMPopt, state);

    % Estimate new X parameters
    [meanX, nuX] = BG_update(EMstate.rhatOpt, EMstate.rvarOpt, meanX, nuX, EMopt);

    % Reinitialize AMP estimates
    BiGAMPopt.xhat0 = xhat;
    BiGAMPopt.xvar0 = xvar;
    BiGAMPopt.shat0 = state.shatOpt;
    BiGAMPopt.Ahat0 = Ahat;
    BiGAMPopt.Avar0 = Avar;
    BiGAMPopt.step = BiGAMPopt.stepMin;
    
    % Count stop
    stop = stop + 1;
end

theta_Uout = Ahat*xhat;
theta_Vout = Avar*(abs(xhat).^2) + (abs(Ahat).^2)*xvar;



