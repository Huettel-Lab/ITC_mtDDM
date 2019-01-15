function fit_discount_k(sample)
%Fits hyperbolic discount rate
% sample - 1 (primary) or 2 (replication) data

%REVISION HISTORY:
%     brian  03.10.06 written
%     brian  03.14.06 added fallback to FMINSEARCH, multiple fit capability
%     kenway 11.29.06 added CI evaluation for FMINSEARCH (which returns
%     Hessian)
%     joe kable 03.12.09 modified to work with revised version of
%     choice_prob and discount
%     khoi 06.26.09 simplified

    dataPath=pwd; %adapt to your location
    cd(dataPath)
    if sample ==1 %Primary sample
        load('amasinoETAl_behavior.csv') %load primary sample data
        data=amasinoETAl_behavior;
    else % replication sample
        load('amasinoETAl_behavior_rep.csv')
        data=amasinoETAl_behavior_rep;
    end
    
    kvals=[];
    noise=[];
    subj=1:data(end,1);
    for i = 1:length(subj)  %loop over all subjects
        %Only fit a given subject, and exclude non-responses and the immediate $10 option
        %(sanity check to make sure they prefer $10 today to $10 tomorrow)
        ind=data(:,1)==subj(i) & ~isnan(data(:,6)) & data(:,2)~=10; 
        [info,p]=fit_discount_model(data(ind,6),data(ind,2),data(ind,4),data(ind,3),data(ind,5)); 
        if info.b(2)<.0001 && mean(data(ind,6))>.95 %Extremely patient
            kvals=[kvals; -9.5];
        elseif info.b(2)<0 %Inconsistent choice, can't properly fit
            kvals=[kvals; NaN];
        else %Fit within typical range, keep fit
            kvals=[kvals; log(info.b(2))]; %Save log(k) as output
        end
        noise=[noise; info.b(1)]; %Save noise as output
    end
    if sample==1;
        csvwrite('allLogk.csv',kvals) % write log k-vals to csv
       %csvwrite('allnoisek.csv',noise) % write temperature parameter to csv
    else
       csvwrite('allLogk_rep.csv',kvals) % write log k-vals to csv
%     %csvwrite('allnoisek_rep.csv',noise) % write temperature parameter to csv
    end
end

%----- FIT DISCOUNTING MODEL - LOGISTIC FIT OF HYPERBOLIC DISCOUNTING
%     [info,p] = fit_discount_model(choice,v1,d1,v2,d2)
%
%     Fits a probabilistic discounting model to binary choices by maximum likelihood.
%
%     Calls local functions: log-likelihood, choice probability and discount
%
%     INPUTS
%     choice    - Dependent variable. The data should be *ungrouped*,
%                   such that CHOICE is a column of 0s and 1s, where 1 indicates 
%                   a choice of OPTION 2.
%     v1        - values of option 1 (ie, sooner option)
%     d1        - delays of option 1
%     v2        - values of option 2 (ie, later option)
%     d2        - delays of option 2
%
%     OUTPUTS
%     info      - data structure with following fields:
%                     .nobs      - number of observations
%                     .nb        - number of parameters
%                     .optimizer - function minimizer used
%                     .exitflag  - see FMINSEARCH
%                     .a         - fitted parameters; note that for all the
%                                  available models, the first element of A
%                                  is a noise term for the choice
%                                  function, the remaining elements are
%                                  parameters for the selected discount
%                                  functions. eg., for dfn='exp', A(2) is
%                                  the time constant of the exponential decay.
%                     .LL        - log-likelihood evaluated at maximum
%                     .LL0       - restricted (minimal model) log-likelihood
%                     .AIC       - Akaike's Information Criterion 
%                     .BIC       - Schwartz's Bayesian Information Criterion 
%                     .r2        - pseudo r-squared
%                   This is a struct array if multiple discount functions are fit.
%     p         - Estimated choice probabilities evaluated at the values & 
%                   delays specified by the inputs v1, v2, p1, p2. This is
%                   a cell array if multiple models are fit.
%

function [info,p] = fit_discount_model(choice,v1,d1,v2,d2)

    nobs = length(choice);
    alpha = 0.05;
    b0 = [.09 .008]; %noise, kval starting points--can test different values if needed 
    % Fit model, attempting to use FMINUNC first, then falling back to FMINSEARCH
    if exist('fminunc','file')
       try
          optimizer = 'fminunc';
          OPTIONS = optimset('Display','off','LargeScale','off','MaxIter',1000);
          [b,negLL,exitflag,convg,g,H] = fminunc(@local_negLL,b0,OPTIONS,choice,v1,d1,v2,d2);
          if exitflag ~= 1 % trap occasional linesearch failures
             optimizer = 'fminsearch';
             fprintf('FMINUNC failed to converge, switching to FMINSEARCH\n');
          end         
       catch
          optimizer = 'fminsearch';
          fprintf('Problem using FMINUNC, switching to FMINSEARCH\n');
       end
    else
       optimizer = 'fminsearch';
    end

    if strcmp(optimizer,'fminsearch')
       optimizer = 'fminsearch';
       OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,... 
          'DiffMinChange',1e-4,'Maxiter',1000,'MaxFunEvals',2000);
       [b,negLL,exitflag,convg] = fminsearch(@local_negLL,b0,OPTIONS,choice,v1,d1,v2,d2);
    end

    if exitflag ~= 1
       fprintf('Optimization FAILED, #iterations = %g\n',convg.iterations);
    else
       fprintf('Optimization CONVERGED, #iterations = %g\n',convg.iterations);
    end

    % Choice probabilities (for OPTION 2)
    p = choice_prob(v1,d1,v2,d2,b);
    % Unrestricted log-likelihood
    LL = -negLL;
    % Restricted log-likelihood
    LL0 = sum((choice==1).*log(0.5) + (1 - (choice==1)).*log(0.5));

    % Confidence interval, requires Hessian from FMINUNC
    try
        invH = inv(-H);
        se = sqrt(diag(-invH));
    catch
    end

    info.nobs = nobs;
    info.nb = length(b);
    info.optimizer = optimizer;
    info.exitflag = exitflag;
    info.b = b;

    try
        info.se = se;
        info.ci = [b-se*norminv(1-alpha/2) b+se*norminv(1-alpha/2)]; % Wald confidence
        info.tstat = b./se;
    catch
    end

    info.LL = LL;
    info.LL0 = LL0;
    info.AIC = -2*LL + 2*length(b);
    info.BIC = -2*LL + length(b)*log(nobs);
    info.r2 = 1 - LL/LL0;
end

%----- LOG-LIKELIHOOD FUNCTION
function sumerr = local_negLL(beta,choice,v1,d1,v2,d2)
    p = choice_prob(v1,d1,v2,d2,beta);

    % Trap log(0)
    ind = p == 1;
    p(ind) = 0.9999;
    ind = p == 0;
    p(ind) = 0.0001;
    % Log-likelihood
    err = (choice==1).*log(p) + (1 - (choice==1)).*log(1-p);
    % Sum of -log-likelihood
    sumerr = -sum(err);
end

%----- CHOICE PROBABILITY FUNCTION - LOGIT
%     p = choice_prob(v1,d1,v2,d2,beta);
%
%     INPUTS
%     v1    - values of option 1 (ie, sooner option)
%     d1    - delays of option 1
%     v2    - values of option 2 (ie, later option)
%     d2    - delays of option 2
%     beta  - parameters, noise term (1) and discount rate (2)
%
%     OUTPUTS
%     p     - choice probabilities for the **OPTION 2**
%
%     REVISION HISTORY:
%     brian lau 03.14.06 written
%     khoi 06.26.09 simplified 

function [p,u1,u2] = choice_prob(v1,d1,v2,d2,beta)
    u1 = discount(v1,d1,beta(2:end));
    u2 = discount(v2,d2,beta(2:end));

    % logit, smaller beta = larger error
    p = 1 ./ (1 + exp(beta(1).*(u1-u2)));
end

%----- DISCOUNT FUNCTION - HYPERBOLIC
%     y = discount(v,d,beta)
%
%     INPUTS
%     v     - values
%     d     - delays
%     beta  - discount rate
%
%     OUTPUTS
%     y     - discounted values
%
%     REVISION HISTORY:
%     brian lau 03.14.06 written
%     khoi 06.26.09 simplified 

function y = discount(v,d,beta)
    y = v./(1+beta(1).*d);
end