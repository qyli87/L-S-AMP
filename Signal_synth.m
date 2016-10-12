%--------------------------------------------------------
% Yangqing Li 
% 31 Oct 2015
% -------------------------------------------------------

function [x_true, y, A, sig2e] = Signal_synth(SignalParams)
%% Unpack the SigGenParams structure
% Dimension of the signal matrix X: NxT
N = SignalParams.N;
T = SignalParams.T;
% # of measurements
M = SignalParams.M;                    
% Sparsity parameter Pr{s_n = 1}
lambda = SignalParams.lambda;
% Rank of X, must <= min(K,T)
Ra = SignalParams.rank;
% SNR (dB)
SNRmdB = SignalParams.SNRmdB;
%  'Gaussian' or 'Bernoulli' projection
Atype1 = SignalParams.Atype1;
%  'different' or 'common' projection
Atype2 = SignalParams.Atype2;

%% Create an original L&S matrix
% Number of non-zero coeffs
K = round(lambda*N);
% Force at least one non-zero coeff
if K == 0, K = 1; end;
% Indices of non-zero coeffs
[~, locs] = sort(rand(N,1));

% Generate factor matrices
U = orth(randn(K,Ra));
V = orth(randn(T,Ra)).';

% Generate singular values
alpha = 0.5;
svec(1:Ra) = exp(-alpha*(1:Ra));

% Build theta based on specified singular values
theta = U*diag(svec)*V;

% Normalization
norm_theta = norm(theta,'fro')/sqrt(K*T);
theta = theta/norm_theta;

% Create original synthetic L&S signal
theta_true = zeros(N,T);
theta_true(sort(locs(1:K)),:) = theta;
x_true = cell(1,T);
for t = 1:T
    x_true{t} = theta_true(:,t);
end

%% Create the measurement matrices   
A = cell(1,T);

if isequal(Atype1,'Gaussian')
    A{1} = randn(M,N); % random Gaussian matrix
else
    A{1} = 2*randi([0 1],M,N)-1; % random ¡À1 matrix
end
for n=1:N, A{1}(:,n) = A{1}(:,n)/norm(A{1}(:,n)); end;
% Compute total amount of signal power
total_signal_power = norm(A{1}*x_true{1})^2;

for t = 2:T
    if isequal(Atype2,'com')
        A{t} = A{1};
    else
        if isequal(Atype1,'Gaussian')
            A{t} = randn(M,N); % random Gaussian matrix
        else
            A{t} = 2*randi([0 1],M,N)-1; % random ¡À1 matrix
        end
    end
    for n=1:N, A{t}(:,n) = A{t}(:,n)/norm(A{t}(:,n)); end;
    % Compute total amount of signal power
    total_signal_power = total_signal_power + norm(A{t}*x_true{t})^2;
end

%% Construct noisy measurements
% Determine the AWGN variance
sig2e = (total_signal_power/M/T)*10^(-SNRmdB/10);

% Compute noisy measurements
e = cell(1,T);
y = cell(1,T);
for t = 1:T
    e{t} = sqrt(sig2e)*randn(M,1); %AWGN
    y{t} = A{t}*x_true{t} + e{t};
end

