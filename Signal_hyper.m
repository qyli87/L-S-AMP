%--------------------------------------------------------
% Full image is cropped down to the scene of 20*20 pixels.
% Full data available at 
% http://hyperspectral.ee.uh.edu/?page_id=459.
% 
% Yangqing Li 
% 31 Oct 2015
% --------------------------------------------------------

function [f_true, y, A, sig2e] = Signal_hyper(DataParams, fname, Ns)
%% Unpack the SigGenParams structure
N = DataParams.N;             % # of spectrum bands
M = DataParams.M;             % # of measurements
T = DataParams.T;             % # of pixels
SNRmdB = DataParams.SNRmdB;   % SNR (dB)
Atype = DataParams.Atype;     %  'Gaussian' or 'Bernoulli' projection

%% Placeholder initializations
f_true = cell(1,T);
A = cell(1,T);
Phi = cell(1,T);
w = cell(1,T);
y = cell(1,T);

%% load agriculture-oriented data
if isequal(fname,'agr')
    Hdata=load('HyperData.mat');
    
    %the pixel coordinates of the cropped region is (1:20, 1:20)
    Hdata=Hdata.Hyper(1:20, 1:20, :);

    [Tr, Tl, n] = size(Hdata);

    l=(Ns-1)*T+1;
    x = reshape(Hdata,Tr*Tl,n);
    F = x(l:l+T-1,:)';

    % Coarse normalization
    F = F/2500; 
end

%% load urban data
if isequal(fname,'urban')
    Hdata=load('uh_data.mat');
    
    %the pixel coordinates of the cropped region is (1:20, 1:20)
    Hdata=Hdata.uh_data(1:20, 1:20, :);
    
    [Tr, Tl, n] = size(Hdata);

    l=(Ns-1)*T+1;
    x = reshape(Hdata,Tr*Tl,n);
    F = x(l:l+T-1,:)';

    % Coarse normalization
    F = F/10000;
end

%% Create projection matrices and output noisy measurements
total_output_power = 0;    % Empirical overall power of output measurements 

for t = 1:T
    % True signal
    f_true{t} =  F(:,t);
    
    % Create projection matrices
    if isequal(Atype,'Gaussian')
        Phi{t} = randn(M,N); % Gaussian matrix
    else
        Phi{t} = 2*randi([0 1],M,N)-1; % Bernoulli (random ±1) matrix
    end
    for n=1:N, Phi{t}(:,n) = Phi{t}(:,n)/norm(Phi{t}(:,n)); end;
    % Compute empirical overall signal power
    total_output_power = total_output_power + norm(Phi{t}*f_true{t})^2;
end

% Determine the AWGN variance
sig2e = (total_output_power/M/T)*10^(-SNRmdB/10);

% Compute noisy measurements
for t = 1:T
    w{t} = sqrt(sig2e)*randn(M,1);  % AWGN
    y{t} = Phi{t}*f_true{t} + w{t};
    A{t} = Phi{t}*dctmtx(N);        % A = Phi * Psi
end
