%--------------------------------------------------------
% A simple example to resolve the CS problem of L&S model
% using L&S-AMP. 
% Author:	Yangqing Li     2016
%--------------------------------------------------------

clear all
close all
clc

addpath('gampmatlab');

% Initialize random number stream
randn('state',0); rand('state',0); %#ok<RAND>

M = 35;
CNMSE = zeros(1,length(M));

for ii = 1:length(M)
    Acc_CNMSE = 0;  % Accumulated CNMSE
    for iii = 1:1  % The number of trials for each simulation settings
        SigGenObj.N = 256;               % Dimension of the signal matrix X
        SigGenObj.T = 50;                % # of columes
        SigGenObj.M = M(ii);             % # of measurements
        SigGenObj.lambda = 1/8;          % Pr{s_n = 1} for synthetic X
        SigGenObj.rank = 3;              % Rank for synthetic X
        SigGenObj.SNRmdB = 25;           % SNR (dB)
        SigGenObj.Atype1 = 'Bernoulli';  %  'Gaussian' or 'Bernoulli'
        SigGenObj.Atype2 = 'diff';       %  'diff' or 'com' projection
        
        [x_true, y, A] = Signal_synth(SigGenObj);
        
        x_hat = LS_AMP(y, A);
        
        MMSEs = sum(abs([x_true{:}]-[x_hat{:}]).^2, 1)./sum(abs([x_true{:}]).^2, 1);
        CNMSE_temp = sum(MMSEs)/SigGenObj.T;
        Acc_CNMSE = Acc_CNMSE + CNMSE_temp;
    end
    CNMSE(ii) = Acc_CNMSE/iii;
end

disp(['CNMSE: ' num2str(10*log10(CNMSE)) 'dB']);
plot(M/SigGenObj.N,10*log10(CNMSE),'o-r');
legend('L&S-AMP')
xlabel('{\it M/N} [dB]');
ylabel('CNMSE [dB]');
grid

