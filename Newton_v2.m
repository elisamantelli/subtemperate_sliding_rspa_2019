function [xout, error_flag, faux] = Newton_v2(f,df,x,parameters,srchparams)
%[xout, error_flag, faux] = Newton(f,df,x,parameters,srchparams)
%Newton: Newton's method for finding zeros of function f. Works the same as
%Newton but passes only parameter x to functions f, df if no input variable
%'parameters' is provided
%
%   Input arguments:
%   f:      function handle of function whose roots are to be found
%           If Newton requires one or two output arguments
%           f must take the form
%               fout = f(x,parmaeters)
%           If Newton requires three output arguments, f must take the form
%               [fout, faux] = f(x,parameters)
%           The Newton iteration will solve for fout=0; faux can be
%           additional outputs computed within f that may be used elsewhere
%   df:     function handle of Jacobian matrix of f
%           df must take the form
%               fout = df(x,parameters)
%   x:      initial guess; must have correct dimensions for an input
%           argument of f, df
%   parameters: optional parameter structure to be passed to f, df. If
%           parameters is not passed to Newton_v2, an empty array is passed
%           to f, df instead. If parameters os passed to Newton_v2 but has
%           a field 'noparameters' set to 'true', no parameters structure
%           is passed to f, df.
%   srchparams: optional structure containing search parameters for Newton
%               iteration. This can contain the fields:
%               itmax:  maximum number of iterations
%               toldelta:   error bound on x; last Newton step must be <=
%                           0.5*toldelta
%               tolgrad:    error bound on df = 0; 0.5*norm(df,2) <=
%                           tolgrad to terminate iteration
%               verbose:    0 produces no runtime output
%                           1 displays step size and current value of
%                           function f
%                           2 also displays current iterate
%
%   The algorithm used is Newton's method as a default; intend to upgrade
%   to line search algorithm
%
%   Output:
%   xout:       root of f
%   error_flag: returns logical one if convergence to prescribed tolerance
%               has not been achieved
%   faux:       optional second output of function f, if auxiliary
%               calculations within f generate additional output
%
%Christian Schoof, May 2013

%NOTE: Check against Numerical Recipes for more efficient version; see
%below for comments regarding line search.

    %set up Newton iteration parameters
    if nargin == 5
        if isfield(srchparams,'itmax')
            itmax = srchparams.itmax;
        else
            itmax = 25;
        end
        if isfield(srchparams,'tolF')
            tolF = srchparams.tolF;
        else
            tolF = sqrt(eps);
        end
        if isfield(srchparams,'toldelta')
            toldelta = srchparams.toldelta;
        else
            toldelta = sqrt(eps);
        end
        if isfield(srchparams,'verbose')
            verbose = srchparams.verbose;
        else
            verbose = 0;
        end
    else
        itmax = 25;
        tolF = sqrt(eps);
        toldelta = sqrt(eps);
        verbose = 0;
    end

    %Deal with missing input parameters
    if nargin < 4
        parameters.noparameters = true;
    end
    
    %set error_flag to false if required
    if nargout > 1, error_flag = false; end
    
    %set up first Newton iteration step
    if nargout < 3
        if isfield(parameters,'noparameters') && parameters.noparameters
            F = f(x);
        else
            F = f(x,parameters);
        end
    else
        if isfield(parameters,'noparameters') && parameters.noparameters
            [F, faux] = f(x);
        else
            [F, faux] = f(x,parameters);
        end
    end
    if verbose, disp(strcat('Initial norm of f =',num2str(norm(F,2)))), end
    Dx = toldelta + eps;
    iter = 1;
    while iter < itmax && (norm(F) >= tolF || norm(Dx) >= toldelta)
        if isfield(parameters,'noparameters') && parameters.noparameters
            DF = df(x);
        else
            DF = df(x,parameters);
        end
        Dx = -DF\F;
        x = x + Dx;
        if nargout < 3
            if isfield(parameters,'noparameters') && parameters.noparameters
                F = f(x);
            else
                F = f(x,parameters);
            end
        else
            if isfield(parameters,'noparameters') && parameters.noparameters
                [F, faux] = f(x);
            else
                [F, faux] = f(x,parameters);
            end
        end
        if verbose, disp(strcat('Updated norm of f =',num2str(norm(F,2)))), disp(strcat('Size of step Dx =',num2str(norm(Dx,2)))), end
        iter = iter + 1;
        if verbose == 2, disp('Current iterate ='), disp(x), end
    end
    %output warning, set error_flag if convergence not achieved
    if iter >= itmax, warning('Newton did not converge'), if nargout == 2, error_flag = true; end, end
    xout = x;
end