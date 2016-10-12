% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu

function [x_hat, StateOut] = GAMP_frame(y, A, Params, Options, StateIn)
%% Variables Declaration
% Initialize dimensionality of data
[M, N] = size(A);

% Declare important algorithmic constants
maxprob = 1 - 1e-15; 	% OutMessages.pi(n) cannot exceed this
minprob = 1 - maxprob;  % OutMessages.pi(n) cannot fall below this
tol = 1e-5;             % Early termination tolerance

% Unpack 'params'
pi = Params.pi;
xi = Params.xi;
psi = Params.psi;
sig2e = Params.sig2e;

% Unpack 'Options'
iter = Options.iter;     

% Unpack 'StateIn'
c = StateIn.c;
z = StateIn.z;
Mu = StateIn.mu;

% Initialize flag for early termination
tol_flag = 0;


%% Execute GAMP algorithm

% Declare static runtime variables
gam_A = ((1 - pi)./pi);             % Constant "A" in gamma definition

oneN = ones(N,1);

for k = 1:iter      % Number of AMP updates to perform
    
    % Update x-to-g messages
    Phi = A'*z + Mu;            % AMP method for calculating Phi (Eqn. (A4))
    % Set-up variables to compute gamma
    gam_B = (c*oneN + psi);
    gam_C = c./psi;
    
    if (k == iter) || (nargin == 0) || (tol_flag == 1)
        phi = Phi;
        c = c*oneN;     % Vectorize the scalar c
        % Compute gamma vector
        gam_exp = exp(-(1/2)*(abs(phi + gam_C.*xi).^2 - ...
            gam_C.*(1 + gam_C).*(abs(xi).^2))./(gam_C.*gam_B));
        % Clip the Gamma exponent values that are too large or small
        % due to numerical precision issues
        gam_exp(gam_exp == 0) = 1/realmax;
        gam_exp(gam_exp == inf) = realmax;
        gamma = gam_A.*sqrt(gam_B./c).*gam_exp;
        % Compute mu, x_hat, and v
        gamma_redux = 1./(1 + gamma);
        mu = gamma_redux.*((phi.*psi + c.*xi) ./ gam_B);
        
        % Compute outgoing messages
        v = gamma_redux.*(c.*psi)./gam_B + gamma.*abs(mu).^2;
        % Correct Inf*0 numerical precision problem
        v(gamma == inf) = 0;
        % Calculate outgoing pi
        StateOut.pi = max(min(1./(1 + sqrt((gam_B./c)).*gam_exp), ...
            maxprob*ones(N,1)), minprob*ones(N,1));
        % Next, store outgoing conditional variances
        StateOut.v = v;
        
        % Break from "for" loop if tol_flag == 1
        if tol_flag == 1, break; end
    end

    % Update g-to-x messages
    % Calculate gamma
    Gam_exp = exp(-(1/2)*(abs(Phi + gam_C.*xi).^2 - ...
        gam_C.*(1 + gam_C).*(abs(xi).^2))./(gam_C.*gam_B));  % Eqn. xx
    % Clip the Gamma exponent values that are too large or small
    % due to numerical precision issues
    Gam_exp(Gam_exp == 0) = 1/realmax;
    Gam_exp(Gam_exp == inf) = realmax;
    Gamma = gam_A.*sqrt(gam_B./c).*Gam_exp;
    % Compute means and variances of x-to-g messages
    Gamma_redux = 1./(1 + Gamma);
    Mu_old = Mu;    % Move mu vector to old slot
    Mu = Gamma_redux.*((Phi.*psi + c.*xi) ./ gam_B);    % Eqn. (A5)
    V = Gamma_redux.*(c.*psi)./gam_B + Gamma.*abs(Mu).^2;    % Eqn. (A6)
    V(Gamma == inf) = 0;
        
    % Calculate F', c, and z
    F_prime = V./c;
    c = sig2e + (1/M)*V.'*oneN;     % Eqn. (A7)
    z = y - A*Mu + (1/M)*sum(F_prime)*z;        % Eqn. (A8)
        
    % Check for early termination
    if norm(Mu_old - Mu)^2/N < tol && k > 1
        tol_flag = 1;
    end
end

% Save final estimate of x
x_hat = mu;

% Save final states of AMP variables
StateOut.phi = phi;
StateOut.c = mean(c);

try
    StateOut.z = z;
    StateOut.mu = Mu;
catch
    StateOut.z = Z;
    StateOut.mu = NaN;  % Mu don't need initialization
end

end     % End of main GAMP function