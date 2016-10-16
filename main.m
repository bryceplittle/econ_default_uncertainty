%  Conventions:
%  Index variables are permanently mapped to certain state variables.
%       i:  income today
%       j:  income tomorrow
%       l:  debt today
%       k:  debt tomorrow
%       m:  volatility today
%       n:  volatility tomorrow
%       s:  value functions

clear, clc;

%% Parameters

beta    = 0.98;
gamma   = 2;
psi     = 0.75;
phi     = 0.969;
r       = 0.017;
theta   = 0.282;


%% Tauchen Discretization for Volatility

n_sig   	= 5;
rho_sig     = 0.95;
std_sig     = 0.06;
sig_bar     = exp(0);

[logsig_grid, P_sig] = tauchen(n_sig, (1-rho_sig)*sig_bar, rho_sig, std_sig, 1.5);

sig_grid = exp(logsig_grid);


%% Modified Tauchen Discretization for Income

n_y     	= 21;
mu_y        = 0;
rho_y       = 0.945;
std_y       = 0.025;
P_y         = zeros(n_y, n_y, n_sig);
logy_grid   = linspace(-3*sqrt(std_y^2/(1 - rho_y^2)), 3*sqrt(std_y^2/(1 - rho_y^2)), n_y)';

for m = 1:n_sig
    P_y(:,:,m) = tauchen_mod(logy_grid, (1-rho_y)*mu_y, rho_y, std_y*sig_grid(m));
end

y_grid      = exp(logy_grid);
pi_y        = zeros(n_y, 1, n_sig);
mu_y        = zeros(1, 1, n_sig);
y_def_grid  = zeros(n_y, 1, n_sig);

for m = 1:n_sig
    temp            = P_y(:,:,m)^1e10;
    pi_y(:,:,m)     = temp(1,:)';
    mu_y(:,:,m)     = pi_y(:,:,m)'*y_grid;
    y_def_grid(y_grid >  phi*mu_y(:,:,m),1,m)   = phi*mu_y(:,:,m);
    y_def_grid(y_grid <= phi*mu_y(:,:,m),1,m)   = y_grid(y_grid <= phi*mu_y(:,:,m),:);
end


%% Asset Grid

b_min       = -0.7;                                                           % lower bound of debt
b_max       =  0.2;                                                           % upper bound of debt
b_grid      = b_min:0.02:b_max;                                               % equally spaced grid of debt states
n_b     	= size(b_grid,2);                                                 % number of points in grid for debt states


%% Guess Price Function

q 		= ones(n_b, n_y, n_sig)*(1/(1+r));
def 	= zeros(n_b, n_y, n_sig);


%% Guess Value Functions

v       = zeros(n_b, n_y, n_sig);
v_def   = zeros(n_y, 1, n_sig);
for m = 1:n_sig
    for i = 1:n_y
        v(:,i,m)        = (b_grid + y_grid(i)).^(1-1/psi)';
        v_def(i,1,m)    = (y_def_grid(i,1,m)).^(1-1/psi);
    end
end


%% Value Function Iteration

count_q     = 1;
q_diff 		= 1;
tol 		= 1e-1;
count_max 	= 200;

while (count_q < 100  && q_diff > tol)
    
    q_old = q;
    
    [v,v_def,b_pol] = solve_valfun(v,v_def,y_grid,y_def_grid,b_grid,q_old,beta,psi,gamma,theta,n_b,n_y,n_sig,P_sig,P_y,tol,count_max);
    
	% optional: use a c++ mex function to solve value function iteration
    % [v,v_def,b_pol] = solve_valfun_mex(v,v_def,y_grid,y_def_grid,b_grid,q_old,beta,psi,gamma,theta,n_b,n_y,n_sig,P_sig,P_y,tol,count_max);
    
    for m = 1:n_sig
        def(:,:,m) = v(:,:,m) < repmat(v_def(:,:,m)', n_b, 1);
    end
    
    [q] = solve_px(def,r,n_b,n_y,n_sig,P_sig,P_y);
    
    count_q = count_q + 1;
    q_diff  = max(max(max(abs(q - q_old))));
	
	fprintf('px iteration: %s \n', convertnum(count_q));
    fprintf('px fn error: %s \n', convertnum((q_diff/tol)*1e2)); 
end