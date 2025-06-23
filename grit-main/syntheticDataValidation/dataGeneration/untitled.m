function [scdata,Tgrid,GT] = linData_custom(simOpts, A_custom, b_custom, x0_custom, Tgrid)
% Modified version of linData that accepts custom A, b, x0 parameters
% 
% Inputs:
%   simOpts - simulation options structure
%   A_custom - custom A matrix (n_genes x n_genes)
%   b_custom - custom b vector (n_genes x 1)  
%   x0_custom - custom initial state (n_genes x 1)
%   Tgrid - optional custom time grid

% Use custom parameters instead of getLinSystem
A = A_custom;
b = b_custom;
x0 = x0_custom;
n_genes = length(x0);

% Create ground truth network from A matrix
GT = abs(A) > 0.01;
GT = GT - diag(diag(GT));  % Remove diagonal

dt = simOpts.dt;

clear('tim')
if exist('Tgrid', 'var') && ~isempty(Tgrid)
    tim(1,:) = Tgrid;
    tim(2,:) = (0:length(Tgrid)-1)*simOpts.nc + 1;
    tim(3,:) = (1:length(Tgrid))*simOpts.nc;
else
    tim(1,:) = simOpts.d*(0:simOpts.T-1);
    tim(2,:) = (0:simOpts.T-1)*simOpts.nc + 1;
    tim(3,:) = (1:simOpts.T)*simOpts.nc;
end

%Initial distribution covariance
rng(simOpts.seed)
P = randn(n_genes, n_genes);
P = P*P';
dp = diag(P).^-.5;
P = (x0.^.5.*((dp.*P).*dp')).*x0'.^.5+.3*eye(n_genes);
CP = chol(P)';

%Steady-state calculation for noise intensity scaling
% Calculate steady-state: A*xss + b = 0, so xss = -A\b
if det(A) ~= 0
    xss = -A\b;
else
    xss = x0;  % Fallback to initial state if A is singular
end

% Ensure steady-state is positive
xss = max(abs(xss), 0.1);

Dnoise = (-2*diag(A).*xss).^.5;
% Handle complex numbers (if diagonal elements are positive)
Dnoise = real(Dnoise);
Dnoise(Dnoise <= 0) = 0.1;  % Set minimum noise level

if isfield(simOpts,'noisetype') && strcmp(simOpts.noisetype,'uniform')
    Dnoise = mean(Dnoise)*ones(n_genes,1);
end

%Initial distribution
Y = x0 + simOpts.Ks*CP*randn(n_genes,tim(3,end));
Y = max(Y, 0);  % Ensure non-negative values

for jtim = 2:size(tim,2)
    
    if simOpts.cont
        %Discretisation for continuous time
        nt = round(tim(1,jtim)/dt);
        for jt = 1:nt
            for jcell = tim(2,jtim):tim(3,jtim)
                dx = A*Y(:,jcell) + b;
                Y(:,jcell) = Y(:,jcell) + dt*dx + simOpts.Kn*dt^.5/2^.5*Dnoise.*sum(randn(n_genes,2),2);
                Y(:,jcell) = max(Y(:,jcell),0);
            end
        end
    else
        %Discrete-time system
        for jt = 2:jtim
            dt_step = tim(1,jt)-tim(1,jt-1);
            Y(:,tim(2,jtim):tim(3,jtim)) = (eye(size(A,1)) + dt_step*A) * Y(:,tim(2,jtim):tim(3,jtim)) + dt_step*b + simOpts.Kn*dt_step^.5*Dnoise.*randn(size(A,1),tim(3,jtim)-tim(2,jtim)+1);
            Y(:,tim(2,jtim):tim(3,jtim)) = max(Y(:,tim(2,jtim):tim(3,jtim)), 0);
        end
    end
    
end

%Put data into required format
Tgrid = tim(1,:);
clear('scdata')
for jt = 1:size(tim,2)
    scdata{jt} = Y(:,tim(2,jt):tim(3,jt));
end

end