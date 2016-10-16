function [q] = solve_px(def,r,n_b,n_y,n_sig,P_sig,P_y)
q = ones(n_b, n_y, n_sig)*(1/(1+r));
for l = 1:n_b
    for i = 1:n_y
        for m = 1:n_sig
            for n = 1:n_sig
                for j = 1:n_y
                    q(l,i,m) = q(l,i,m) - P_sig(m,n)*P_y(i,j,n)*def(l,j,n)/(1+r);
                end
            end
        end
    end
end
end

