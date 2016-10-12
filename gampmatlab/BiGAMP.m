function [xhatOpt, xvarOpt, AhatOpt,AvarOpt,...
    EMstate,state] = BiGAMP(gX, gA, gOut, opt, state)
% BiGAMP:  Bilinear Generalized Approximate Message Passing
%
% The BiG-AMP algorithm is intended for the estimation of
% random matrices A and X observed through the Markov Chain
%
%   X,A -> Z = A*X -> Y,
%
% where the components of X and A are independent and the mapping Z -> Y is
% separable. X is NxL, A is MxN, and Z,Y are consequently MxL.
%
% The function takes five arguments:
%
% gX:  An input estimator derived from the EstimIn class
%    based on the input distribution p_x_{nl}(x_nl).
% gA:  An input estimator derived from the EstimIn class
%    based on the input distribution p_a_{mn}(a_mn).
% gOut:  An output estimator derived from the EstimOut class
%    based on the output distribution p_{Y|Z}(y_ml|z_ml).
% opt:  A set of options of the class BiGAMPOpt.
% state (optional): A structure containing all the values needed to warm
%   start BiG-AMP
%
%The function returns final values for numerous estimated quantities along
%with the optional:
%
% EMstate: Several state variables useful when doing EM learning of BiG-AMP
%   parameters
% state: The values of all parameters required to warm start the algorithm


%% Setup

% Get options
nit     = opt.nit;              % number of iterations
nitMin  = opt.nitMin;           % minimum number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
stepFilter = opt.stepFilter;    % step filter setting, <1 for no effect
adaptStep = opt.adaptStep;      % adaptive step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
stepWindow = opt.stepWindow;    % step size check window size
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
pvarStep = opt.pvarStep;        % incldue step size in pvar/zvar
uniformVariance =...
    opt.uniformVariance;        % use scalar variances
compVal = adaptStep;            % only compute cost function for adaptive
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed bad steps
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures
zvarToPvarMax = opt.zvarToPvarMax;  % maximum zvar/pvar ratio

%This option was included for testing purposes and should not be modified
gainMode = opt.gainMode;
if gainMode ~= 1
    warning(...
        'Running with gainMode other than 1 is not recommended') %#ok<*WNTAG>
end

%Determine requested outputs
saveEM = true;
saveState = true;

%Get problem dimensions
M = opt.M;
L = opt.L;
N = opt.N;

%Assign Avar and xvar mins
xvarMin = opt.xvarMin;
Avarmin = opt.AvarMin;

%% Initialization

%Check for provided state
if nargin < 5
    state = [];
end

if isempty(state) %if no state is provided
    
    %Initialize Avar
    [~,Avar] = gA.estimInit();
    
    %Initialize Xvar
    [xhat,xvar] = gX.estimInit();
    
    %Handle case of scalar input distribution estimator on xhat
    if (length(xhat) == 1)
        xhat = gX.genRand([N L]);
    end
    
    
    %Handle case of scalar input distribution estimator on Ahat
    Ahat = gA.genRand([M N]);
    
    
    %Initialize valIn
    valIn = -inf;
    
    %Replace these defaults with the warm start values if provided in the
    %options object
    if ~isempty(opt.xhat0)
        xhat = opt.xhat0;
        %If warm starting, set valIn to be a negative inf to avoid problems
        valIn = -inf;
    end
    if ~isempty(opt.xvar0)
        xvar = opt.xvar0;
    end
    if ~isempty(opt.Ahat0)
        Ahat = opt.Ahat0;
    end
    if ~isempty(opt.Avar0)
        Avar = opt.Avar0;
    end
    
    %Handle uniform variances
    if (length(Avar) == 1)
        Avar = repmat(Avar,M,N);
    end
    if (length(xvar) == 1)
        xvar = repmat(xvar,N,L);
    end
    
    %Placeholder initializations- values are not used
    xhatBar = 0;
    AhatBar = 0;
    shat = 0;
    svar = 0;
    pvarOpt = 0;
    zvarOpt = 0;
        
    %Scalar variance
    if uniformVariance
        Avar = mean(Avar(:));
        xvar = mean(xvar(:));
    end
    
    %Address warm starting of shat0
    if ~isempty(opt.shat0)
        shat = opt.shat0;
    end
    
    %Init valOpt empty
    valOpt = [];
    
    %Set pvar_mean to unity initially
    pvar_mean = 1;
    
else %Use the provided state information
    
    %A variables
    Ahat = state.Ahat;
    Avar = state.Avar;
    AhatBar = state.AhatBar;
    AhatOpt = state.AhatOpt;
    AhatBarOpt = state.AhatBarOpt;
    
    %X Variables
    xhat = state.xhat;
    xvar = state.xvar;
    xhatBar = state.xhatBar;
    xhatOpt = state.xhatOpt;
    xhatBarOpt = state.xhatBarOpt;
    
    %S variables
    shat = state.shat;
    svar = state.svar;
    shatOpt = state.shatOpt;
    svarOpt = state.svarOpt;
    shatNew = state.shatNew;
    svarNew = state.svarNew;
    
    %Cost stuff
    valIn = state.valIn;
    
    %Variance momentum terms
    pvarOpt = state.pvarOpt;
    zvarOpt = state.zvarOpt;
    
    %Step
    step = state.step;
        
    %Old cost values
    valOpt = state.valOpt;
       
    %Set pvar_mean
    pvar_mean = state.pvar_mean;
end

%Specify minimum variances
pvarMin = opt.pvarMin;

%Placeholder initializations
rhat = 0;
rvar = 0;
qhat = 0;
qvar = 0;

%Cost init
val = zeros(nit,1);
zhatOpt = 0;

%% Iterations

%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;

%Handle first step
if isempty(state)
    step1 = 1;
else
    step1 = step;
end

% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
        
    % Check for final iteration
    if it >= nit
        stop = true;
    end
    
    
    if ~uniformVariance
        %Precompute squares quantities for use, these change on every
        %iteration
        Ahat2 = abs(Ahat).^2;
        xhat2 = abs(xhat).^2;
        
        %Compute pvar-
        zvar = Avar*xhat2 + Ahat2*xvar;
        
        %Compute pvar
        pvar = zvar + Avar*xvar;%R3
    else
        
        %Compute useful quantities
        mxhat2 = norm(xhat,'fro').^2/numel(xhat);
        mAhat2 = norm(Ahat,'fro').^2/numel(Ahat);
        mAvar = norm(Avar(:),1)/numel(Avar);
        mxvar = norm(xvar(:),1)/numel(xvar);
        
        %Compute variances
        zvar = (mAvar*mxhat2 + mAhat2*mxvar)*N;
        pvar = zvar + N*mAvar*mxvar;
    end
    
    %Include pvar step
    if pvarStep
        pvar = step1*pvar + (1-step1)*pvarOpt;
        zvar = step1*zvar + (1-step1)*zvarOpt;
    end
    
    %Update zhat
    zhat = Ahat * xhat;
    
    % Compute log likelihood at the output and add it the total negative
    % K-L distance at the input.
    if (compVal)
        valOut = sum(sum(gOut.loglikey(zhat,pvar)));
        val(it) = valOut + valIn;
    end
    
    % Determine if candidate passed
    if ~isempty(valOpt)
        
        %Check against worst value in last stepWindow good steps
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        
        %Check the step
        pass = (val(it) > min(valOpt(startInd:stopInd))) ||...
            ~adaptStep || (step <= stepMin);
        
    else
        pass = true;
    end

    % If pass, set the optimal values and compute a new target shat and
    % snew.
    if (pass)
        
        %Slightly inrease step size after pass if using adaptive steps
        if adaptStep
            step = stepIncr*step;
        end
        
        % Set new optimal values
        shatOpt = shat;
        svarOpt = svar;
        xhatBarOpt = xhatBar;
        xhatOpt = xhat;
        AhatBarOpt = AhatBar;
        AhatOpt = Ahat;
        pvarOpt = pvar;
        zvarOpt = zvar;
        
        %Bound pvar
        pvar = max(pvar, pvarMin);
        
        %We keep a record of only the succesful step valOpt values
        valOpt = [valOpt val(it)]; %#ok<AGROW>
        
        % Continued output step
        phat = zhat...
            - shat.*(zvar/pvar_mean);%...R4
        %+ abs(shat).^2 .* ( (Ahat .* Avar) * (xhat .* xvar)  );
        
        % Output nonlinear step
        [zhat0,zvar0] = gOut.estim(phat,pvar);%R5 R6
        
        %Compute 1/pvar
        pvarInv = pvar_mean ./ pvar;
        
        %Update the shat quantities
        shatNew = pvarInv.*(zhat0-phat);
        svarNew = pvarInv.*(1-min(zvar0./pvar,zvarToPvarMax));
        
        
        %Scalar Variance
        if uniformVariance
            svarNew = mean(svarNew(:));
        end
        
        %Enforce step size bounds
        step = min([max([step stepMin]) stepMax]);
        
    else
        
        %Check on failure count
        failCount = failCount + 1;
        if failCount > maxBadSteps
            failCount = 0;
            stepMax = max(stepMin,maxStepDecr*stepMax);
        end
        % Decrease step size
        step = max(stepMin, stepDecr*step);
        
        %Check for minimum step size
        if step < stepTol
            stop = true;
        end
    end
    
    
    % Check for convergence if step was succesful
    if pass
        if any(isnan(zhat(:))) || any(isinf(zhat(:)))
            stop = true;
        else
            testVal = norm(zhat(:) - zhatOpt(:)) / norm(zhat(:));
            if (it > 1) && ...
                    (testVal < tol)
                stop = true;
            end
        end
        
        %Set other optimal values- not actually used by iterations
        AvarOpt = Avar;
        xvarOpt = xvar;
        zhatOpt = zhat;
        
        %Save EM variables if requested
        if saveEM
            EMstate.rhatOpt = rhat;
            EMstate.rvarOpt = pvar_mean*rvar;
            EMstate.qhatOpt = qhat;
            EMstate.qvarOpt = pvar_mean*qvar;
            EMstate.zhatOpt = zhat0;
            EMstate.zvarOpt = zvar0;
            EMstate.pvarOpt = pvar;
        end
    end
    
    % Create new candidate shat
    if it > 1 || ~isempty(state)
        step1 = step;
        if stepFilter >= 1
            step1 = step1*it/(it+stepFilter);
        end
    end
    shat = (1-step1)*shatOpt + step1*shatNew;
    svar = (1-step1)*svarOpt + step1*svarNew;
    xhatBar = (1-step1)*xhatBarOpt + step1*xhatOpt;
    AhatBar = (1-step1)*AhatBarOpt + step1*AhatOpt;
    
    %Compute rvar and correct for infinite variance
    if uniformVariance
        mAhatBar2 = norm(AhatBar,'fro')^2/numel(AhatBar);
        rvar = 1/(mAhatBar2*svar*M);
    else
        rvar = 1./((abs(AhatBar).^2)'*svar);%R9
    end
    rvar(rvar > opt.varThresh) = opt.varThresh;
    
    %Update rhat
    if uniformVariance
        switch gainMode
            case 1,
                rGain = 1 - Avar/mAhatBar2;
            case 2,
                rGain = 1 - ...
                    Avar*norm(shat,'fro')^2/numel(shat)/mAhatBar2/svar;
            case 3,
                rGain = 1;
        end
    else
        switch gainMode
            case 1,
                rGain = (1 - (rvar.*(Avar'*svar)));
            case 2,
                rGain = (1 - (rvar.*(Avar'*shat.^2)));
            case 3,
                rGain = 1;
        end
    end
    rGain = min(1,max(0,rGain));
    rhat = xhatBar.*rGain + rvar.*(AhatBar'*shat);%R10
    rvar = max(rvar, xvarMin);
    
    % Input linear step for A
    if uniformVariance
        mxhatBar2 = norm(xhatBar,'fro')^2/numel(xhatBar);
        qvar = 1/(svar*mxhatBar2*L);
    else
        qvar = 1./(svar*(abs(xhatBar).^2)');%R11
    end
    qvar(qvar > opt.varThresh) = opt.varThresh;
    
    
    %Update qhat
    if uniformVariance
        switch gainMode
            case 1,
                qGain = 1 - xvar/mxhatBar2;
            case 2,
                qGain = 1 - xvar*norm(shat,'fro')^2/numel(shat)/mxhatBar2/svar;
            case 3,
                qGain = 1;
        end
    else
        switch gainMode
            case 1,
                qGain = (1 - (qvar.*(svar*xvar')));
            case 2,
                qGain = (1 - (qvar.*(shat.^2*xvar')));
            case 3,
                qGain = 1;
        end
    end
    qGain = min(1,max(0,qGain));
    qhat = AhatBar.*qGain + qvar.*(shat*xhatBar');%R12
    qvar = max(qvar,Avarmin);
    
    % Input nonlinear step
    if compVal
        [xhat,xvar,valInX] = gX.estim(rhat, rvar*pvar_mean);%R13 14
        [Ahat,Avar,valInA] = gA.estim(qhat, qvar*pvar_mean);%R15 16
    else %method may avoid computation if the vals are not needed
        [xhat,xvar] = gX.estim(rhat, rvar*pvar_mean);
        [Ahat,Avar] = gA.estim(qhat, qvar*pvar_mean);
    end
    
    %Scalar variances
    if uniformVariance
        Avar = mean(Avar(:));
        xvar = mean(xvar(:));
    end
    
    %Update valIn
    if compVal
        valIn = sum( valInX(:) ) + sum ( valInA(:) );
    end
    
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end

end

%% Save the state

if saveState
    
    %A variables
    state.Ahat = Ahat;
    state.Avar = Avar;
    state.AhatBar = AhatBar;
    state.AhatOpt = AhatOpt;
    state.AhatBarOpt = AhatBarOpt;
    
    %X Variables
    state.xhat = xhat;
    state.xvar = xvar;
    state.xhatBar = xhatBar;
    state.xhatBarOpt = xhatBarOpt;
    state.xhatOpt = xhatOpt;
    
    %s variables
    state.shat = shat;
    state.svar = svar;
    state.shatOpt = shatOpt;
    state.shatNew = shatNew;
    state.svarOpt = svarOpt;
    state.svarNew = svarNew;
    
    %Cost stuff
    state.valIn = valIn;
    
    %Variance momentum terms
    state.pvarOpt = pvarOpt;
    state.zvarOpt = zvarOpt;
    
    %Step
    state.step = step;
        
    %Old cost values
    state.valOpt = valOpt;
     
    %Problem dimensions
    state.N = N;
    
    %pvar_mean
    state.pvar_mean = pvar_mean;
end



