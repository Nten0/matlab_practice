function [x,flag,relres,iter,resvec] = my_pcg(b,tol,maxit,n_in,a)

if (nargin < 2)
    error(message('MATLAB:pcg:NotEnoughInputs'));
end

n = n_in;
existM1 = 0;
existM2 = 0;
n2b = norm(b); 
x = zeros(n_in,1);

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
%r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
r = b - matrix_mv_5853(a,x);
normr = norm(r);                   % Norm of residual
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('pcg',tol,maxit,0,flag,iter,relres);
    end
    return
end

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1,:) = normr;               % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of minimum residual
rho = 1;
stag = 0;                          % stagnation of the method
moresteps = 0;
maxmsteps = min([floor(n/50),5,n-maxit]);
maxstagsteps = 3;

% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    if existM1
        y = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
        if ~all(isfinite(y))
            flag = 2;
            break
        end
    else % no preconditioner
        y = r;
    end
    
    if existM2
        z = iterapp('mldivide',m2fun,m2type,m2fcnstr,y,varargin{:});
        if ~all(isfinite(z))
            flag = 2;
            break
        end
    else % no preconditioner
        z = y;
    end
    
    rho1 = rho;
    rho = r' * z;
    if ((rho == 0) || isinf(rho))
        flag = 4;
        break
    end
    if (ii == 1)
        p = z;
    else
        beta = rho / rho1;
        if ((beta == 0) || isinf(beta))
            flag = 4;
            break
        end
        p = z + beta * p;
    end
    %q = iterapp('mtimes',afun,atype,afcnstr,p,varargin{:});
    q = matrix_mv_5853(a,p);
    pq = p' * q;
    if ((pq <= 0) || isinf(pq))
        flag = 4;
        break
    else
        alpha = rho / pq;
    end
    if isinf(alpha)
        flag = 4;
        break
    end
    
    % Check for stagnation of the method    
    if (norm(p)*abs(alpha) < eps*norm(x))
        stag = stag + 1;
    else
        stag = 0;
    end
    
    x = x + alpha * p;             % form new iterate
    r = r - alpha * q;
    normr = norm(r);
    normr_act = normr;
    resvec(ii+1,1) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
        %r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        r = b - matrix_mv_5853(a,x);
        normr_act = norm(r);
        resvec(ii+1,1) = normr_act;
        if (normr_act <= tolb)
            flag = 0;
            iter = ii;
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:pcg:tooSmallTolerance'));
                end
                flag = 3;
                iter = ii;
                break;
            end
        end
    end
    if (normr_act < normrmin)      % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    if stag >= maxstagsteps
        flag = 3;
        break;
    end
end                                % for ii = 1 : maxit

% returned solution is first with minimal residual
if (flag == 0)
    relres = normr_act / n2b;
else
    %r_comp = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    r_comp = b - matrix_mv_5853(a,xmin);
    if norm(r_comp) <= normr_act
        x = xmin;
        iter = imin;
        relres = norm(r_comp) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% truncate the zeros from resvec
if ((flag <= 1) || (flag == 3))
    resvec = resvec(1:ii+1,:);
else
    resvec = resvec(1:ii,:);
end

% only display a message if the output flag is not used

