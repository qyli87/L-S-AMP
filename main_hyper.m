%--------------------------------------------------------
% A demo of pixel-based compressive hyperspectral imaging
% using L&S-AMP.
% Data available at 
% http://hyperspectral.ee.uh.edu/?page_id=459.
% 
% Author:	Yangqing Li     2016
%--------------------------------------------------------

clear all
close all
clc

addpath('data'); 
addpath('gampmatlab');

% Initialize random number stream
randn('state',1); rand('state',1); %#ok<RAND>

f_true_image = [];
f_hat_image = [];
CNMSE = 0;

% Choose dataset
fname = 'urban'; %  'agr' or  'urban'

for j = 1:4  % The number of sub-scenes (block partition strategy)
    if isequal(fname,'agr')
        SigGenObj.N = 224;          % # of spectrum bands
    else
        SigGenObj.N = 144;
    end
    SigGenObj.M = 50;               % # of measurements
    SigGenObj.T = 100;              % # of pixels
    SigGenObj.SNRmdB = 25;          % Per-measurement empirical SNR (dB)
    SigGenObj.Atype = 'Bernoulli';  %  'Gaussian' or 'Bernoulli' projection
    
    [f_true, y, A, sig2e] = Signal_hyper(SigGenObj,fname,j);
    
    x_hat = LS_AMP(y, A, 'dct');
    
    % Reconstruct images
    dct_matrix = dctmtx(SigGenObj.N);
    f_hat = cell(1,SigGenObj.T);
    for t = 1:SigGenObj.T
        f_hat{t} = dct_matrix*x_hat{t};
        f_hat_image = [f_hat_image f_hat{t}']; %#ok<AGROW>
        f_true_image = [f_true_image f_true{t}']; %#ok<AGROW>
    end
    
    MMSEs = sum(abs([f_true{:}]-[f_hat{:}]).^2, 1)./sum(abs([f_true{:}]).^2, 1);
    CNMSE_temp = sum(MMSEs)/SigGenObj.T;
    disp(['CNMSE: ' num2str(10*log10(CNMSE_temp)) 'dB']);
    CNMSE = CNMSE + CNMSE_temp;
end
CNMSE = CNMSE/j;

% Visual results
figure('Position',[200 450  600 450]);

subplot(2,2,1)
f = reshape(f_true_image,SigGenObj.N,20,20);
ff = f(50,:,:);
ff = reshape(ff,20,20);
imagesc(ff,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
title('Original image')
xlabel('# 50 spectrum band')
colormap(gray)

subplot(2,2,2)
f_e = reshape(f_hat_image,SigGenObj.N,20,20);
ff_e = f_e(50,:,:);
ff_e = reshape(ff_e,20,20);
imagesc(ff_e,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
title(['L&S-AMP: CNMSE=' num2str(10*log10(CNMSE)) 'dB'])
xlabel('# 50 spectrum band')

subplot(2,2,3)
f = reshape(f_true_image,SigGenObj.N,20,20);
ff = f(100,:,:);
ff = reshape(ff,20,20);
imagesc(ff,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('# 100 spectrum band')

subplot(2,2,4)
f_e = reshape(f_hat_image,SigGenObj.N,20,20);
ff_e = f_e(100,:,:);
ff_e = reshape(ff_e,20,20);
imagesc(ff_e,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('# 100 spectrum band')
