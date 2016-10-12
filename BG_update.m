%--------------------------------------------------------
% Yangqing Li 
% 31 Oct 2015
%--------------------------------------------------------

function [theta, phi] = BG_update(Rhat, Rvar, theta, phi, EMopt)
% Initialize problem dimensions
[N, T] = size(Rhat);

% Calculate posterior means and variances
post_var_scale = phi + Rvar + eps;
gamma_n = (Rhat.*phi+Rvar.*theta)./post_var_scale;
nu_n = Rvar.*phi./post_var_scale;

% Update mean and variance
if EMopt.learn_mean
    if strcmp(EMopt.sig_dim,'joint')
        theta = sum(sum(gamma_n))/N/T;
        theta = repmat(theta,[N T]);
    elseif strcmp(EMopt.sig_dim,'col')
        theta = sum(gamma_n)./N;
        theta = repmat(theta,[N 1]);
    elseif strcmp(EMopt.sig_dim,'row')
        theta = sum(gamma_n,2)./T;
        theta = repmat(theta,[1 T]);
    end
end

if EMopt.learn_var
    if strcmp(EMopt.sig_dim,'joint')
        phi = sum(sum((nu_n+(gamma_n-theta).^2)))/N/T;
        phi = repmat(phi,[N T]);
    elseif strcmp(EMopt.sig_dim,'col')
        phi = sum(nu_n+(gamma_n-theta).^2)./N;
        phi = repmat(phi,[N 1]);
    elseif strcmp(EMopt.sig_dim,'row')
        phi = sum((nu_n+(gamma_n-theta).^2),2)./T;
        phi = repmat(phi,[1 T]);
    end
end

return