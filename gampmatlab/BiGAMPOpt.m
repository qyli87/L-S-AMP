classdef BiGAMPOpt
    % Options for the BiG-AMP.
    
    properties
        
        %Initialize problem dimensions.
        N = [];
        M = [];
        L = [];
        
        %Number of iterations
        nit = 100;
        
        %Minimum number of iterations
        nitMin = 30; %0 for no effect
        
        %Specify a convergence tolerance.
        tol = 1e-8;
                
        %***** Initialization
        xhat0 = [];
        Ahat0 = [];
        xvar0 = [];
        Avar0 = [];
        shat0 = [];
        
        x2hat0 = [];
        x2var0 = [];
        
        %***** Step Size Control
        
        %Logical flag to include a step size in the pvar/zvar calculation.
        %This momentum term often improves numerical performance. On by
        %defualt.
        pvarStep = true;
        
        %Initial step size, or fixed size for non-adaptive steps
        step = 0.05;
        
        % Enable adaptive step size
        adaptStep = true;
        
        % Minimum step size
        stepMin = 0.05;
        
        % Maximum step size
        stepMax = 0.5;
        
        % Multiplicative step size increase, when successful
        stepIncr = 1.1;
        
        % Multiplicative step size decrease, when unsuccessful
        stepDecr = 0.5;
        
        %Maximum number of allowed failed steps before we decrease stepMax,
        %inf for no effect
        maxBadSteps = inf;
        
        %Amount to decrease stepMax after maxBadSteps failed steps, 1 for
        %no effect
        maxStepDecr = 0.5;
        
        %Create a window for the adaptive step size test. Setting this to
        %zero causes it to have no effect. For postive integer values,
        %creats a moving window of this length when checking the step size
        %acceptance criteria. The new value is only required to be better
        %than the worst in this window, i.e. the step size criteria is not
        %required to monotonicaly increase. 
        stepWindow = 1;
        
        % Iterations are termined when the step size becomes smaller
        % than this value. Set to -1 to disable
        stepTol = -1;
        
        %This is a filter value on the steps. It slows down the early steps
        %in a smooth way. Set to less than 1 for no effect, In particular,
        %the steps are computed as step1 = step*it/(it + stepFilter)
        stepFilter = 0;
        
        
        %***** Variance Control
        
        %Switch to enable uniform variance mode. Replaces all variances
        %with scalar approximations
        uniformVariance = false;
        
        %Minimum variances. See code for details of use. 
        pvarMin = 1e-13;
        xvarMin = 0;
        AvarMin = 0;
        zvarToPvarMax = 0.99;   % prevents negative svar, should be near 1
        
        %Variance threshold for rvar and qvar, set large for no effect
        varThresh = 1e6;
        
        %This was included for testing and should not be modified
        gainMode = 1;
        
    end
    
    methods
        % Constructor with default options
        function opt = BiGAMPOpt()
        end
    end
    
end