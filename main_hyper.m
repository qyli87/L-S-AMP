%--------------------------------------------------------
% A demo of pixel-based compressive hyperspectral imaging
% using L&S-AMP. 
% using L&S-AMP. We cropped down the full image to the 
% scene of 20*20 pixels. Full data available at 
% http://hyperspectral.ee.uh.edu/?page_id=459.
% 
% Author:	Yangqing Li     2016
%--------------------------------------------------------

@@ -12,23 +15,23 @@ addpath('data');
addpath('gampmatlab');

% Initialize random number stream
randn('state',0); rand('state',0); %#ok<RAND>
randn('state',1); rand('state',1); %#ok<RAND>

f_true_image = [];
f_hat_image = [];
CNMSE = 0;

% Choose dataset
fname = 'agr'; %  'agr' or  'urban'
fname = 'urban'; %  'agr' or  'urban'

for j = 1:8  % The number of sub-scenes (block partition strategy)
for j = 1:4  % The number of sub-scenes (block partition strategy)
    if isequal(fname,'agr')
        SigGenObj.N = 224;          % # of spectrum bands
    else
        SigGenObj.N = 144;
    end
    SigGenObj.M = 40;               % # of measurements
    SigGenObj.T = 200;              % # of pixels
    SigGenObj.M = 50;               % # of measurements
    SigGenObj.T = 100;              % # of pixels
    SigGenObj.SNRmdB = 25;          % Per-measurement empirical SNR (dB)
    SigGenObj.Atype = 'Bernoulli';  %  'Gaussian' or 'Bernoulli' projection
    
@@ -56,9 +59,9 @@ CNMSE = CNMSE/j;
figure('Position',[200 450  600 450]);

subplot(2,2,1)
f = reshape(f_true_image,SigGenObj.N,40,40);
f = reshape(f_true_image,SigGenObj.N,20,20);
ff = f(50,:,:);
ff = reshape(ff,40,40);
ff = reshape(ff,20,20);
imagesc(ff,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
@@ -67,9 +70,9 @@ xlabel('# 50 spectrum band')
colormap(gray)

subplot(2,2,2)
f_e = reshape(f_hat_image,SigGenObj.N,40,40);
f_e = reshape(f_hat_image,SigGenObj.N,20,20);
ff_e = f_e(50,:,:);
ff_e = reshape(ff_e,40,40);
ff_e = reshape(ff_e,20,20);
imagesc(ff_e,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
@@ -77,18 +80,18 @@ title(['L&S-AMP: CNMSE=' num2str(10*log10(CNMSE)) 'dB'])
xlabel('# 50 spectrum band')

subplot(2,2,3)
f = reshape(f_true_image,SigGenObj.N,40,40);
f = reshape(f_true_image,SigGenObj.N,20,20);
ff = f(100,:,:);
ff = reshape(ff,40,40);
ff = reshape(ff,20,20);
imagesc(ff,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('# 100 spectrum band')

subplot(2,2,4)
f_e = reshape(f_hat_image,SigGenObj.N,40,40);
f_e = reshape(f_hat_image,SigGenObj.N,20,20);
ff_e = f_e(100,:,:);
ff_e = reshape(ff_e,40,40);
ff_e = reshape(ff_e,20,20);
imagesc(ff_e,[min(ff(:)),max(ff(:))])
set(gca,'XTick',[])
set(gca,'YTick',[])
