function [Zprob] = tauchen_mod(y_grid, mu, rho, sigma)
% This code modifies Martin Floden's code for the Tauchen method
% to implement my method of discretizing an AR(1) process with stochastic
% volatility.

if iscolumn(y_grid) == 0
    y_grid = y_grid';
end 

Z     = y_grid;
N     = size(y_grid, 1);
Zprob = zeros(N,N);
a     = (1-rho)*mu;
zstep = (Z(N) - Z(1)) / (N - 1);

for i=2:(N-1)
    Z(i) = Z(1) + zstep * (i - 1);
end 

Z = Z + a / (1-rho);

for j = 1:N
    for k = 1:N
        if k == 1
            Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma);
        elseif k == N
            Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma);
        else
            Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - ...
                         cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
        end
    end
end

end

function c = cdf_normal(x)

    c = 0.5 * erfc(-x/sqrt(2));
	
end

