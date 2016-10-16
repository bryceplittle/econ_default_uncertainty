function [v,v_def,b_pol] = solve_valfun(v,v_def,y_grid,y_def_grid,b_grid,q_old,beta,psi,gamma,theta,n_b,n_y,n_sig,P_sig,P_y,tol,count_max)

v_diff      = 1;
count       = 1;
Ev_def      = complex(zeros(n_y, n_sig));
b_pol       = zeros(n_b, n_y, n_sig);
b_def_node  = find(b_grid == 0);

while (count < count_max && v_diff > tol)
    
    v_old       = v;
    v_def_old   = v_def;
    
    % compute value of not defaulting
    parfor m = 1:n_sig
        for i = 1:n_y
            for k = 1:n_b
                c   = complex(zeros(n_b,1));
                u   = complex(zeros(n_b,1));
                Ev  = complex(zeros(n_b,1));
                for l = 1:n_b
                    for n = 1:n_sig
                        c(l) = complex(y_grid(i) + b_grid(k) - q_old(l,i,m)*b_grid(l));
                        if c(l) > 0
                            u(l)  = c(l)^(1-1/psi);
                        else
                            u(l)  = -Inf;
                        end
                        Ev(l) = beta*(P_sig(m,n)*P_y(i,:,n)*bsxfun(@max, v_old(l,:,n).^(1-gamma)', v_def_old(:,1,n).^(1-gamma))).^((1-1/psi)/(1-gamma));
                    end
                end
                [v(k,i,m), b_pol(k,i,m)] = max((u + Ev).^(1-1/psi));
            end
        end
    end
    
    % compute value of defaulting
    u_def   = (y_def_grid).^(1-1/psi);
    for m = 1:n_sig
        for i = 1:n_y
            for n = 1:n_sig
                Ev_def(i,m) = beta*(P_sig(m,n)*P_y(i,:,n)*(theta*bsxfun(@max, v_old(b_def_node,:,n).^(1-gamma)', v_def_old(:,1,n).^(1-gamma)) + (1-theta)*v_def_old(:,1,n).^(1-gamma))).^((1-1/psi)/(1-gamma));
            end
            v_def(i,m) = (u_def(i,m) + Ev_def(i,m)).^(1-1/psi);
        end;
    end
    
    count   = count + 1;
    v_diff  = max(max(max(max(abs(v - v_old)))), max(max(max(abs(v_def - v_def_old)))));
    
end
end

