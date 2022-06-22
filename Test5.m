clear 
close

% First load the dataset variables
load dataset3.mat

% Now start by stating some definitions
% Motion model noise variance
Q = diag([v_var; w_var]);    % Diagonal 6x6 matrix

% Measurement model noise variance
R = diag(y_var);    % Diagonal 4x4 matrix

Q_k = eye(6,6);

% Set transformation matrix T_cv
T_cv = eye(4);
T_cv(1:3,1:3) = C_c_v;
T_cv(1:3,4) = -C_c_v * rho_v_c_v;

% D traspose matrix creation
DT = [1,0,0,0;0,1,0,0;0,0,1,0];

% Limits
k1 = 1;
k2 = 200;%1714;

for ite = 1:2
    
    store = 1;
    for i = k1:k2
    
        if ite ==1
            % Motion model equations to estimate T_kpre
            if i == k1
                TK = 0;
                % Translation vector
                r_k0 = r_i_vk_i(:,i);
                % Rotation matrix
                psi_k0 = theta_vk_i(:,i);
                C_k0 = vec2rot(psi_k0);
                
                r_k(1:3,store) = r_k0;
                
                % Rotation matrix
                C_k{store,1} = vec2rot(psi_k0);
                
                % Transformation matrix
                T = eye(4);
                T(1:3,1:3) = C_k{store,1};
                T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
                T_kpre{store,1} = T;
                
            else
                Tk = t(i)-t(i-1);   % Variable period at each timestep
                % Translation vector
                d_k = v_vk_vk_i(1:3,i-1) * Tk;
                r_k(1:3,store) = r_k(1:3,store-1) + C_k{store-1,1}.' * d_k;
                
                % Rotation matrix
                psi_k(1:3,store) = w_vk_vk_i(1:3,i-1) * Tk;
                C_k{store,1} = vec2rot(psi_k(1:3,store)) * C_k{store-1,1};
                
                % Transformation matrix
                T = eye(4);
                T(1:3,1:3) = C_k{store,1};
                T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
                T_kpre{store,1} = T;
                
                % W matrix append
                Q_k = blkdiag(Q_k,(Tk^2) * Q);
            end
        else
            
            if i == k1
                TK = 0;
                T_eps{store,1} = vec2tran(eps(store*6-5:store*6,1));
                NewT_op{store,1} =  T_eps{store,1} * T_kpre{store,1};
                T2 = NewT_op{store,1};
                C_k{store,1} = T2(1:3,1:3);
                r = inv(NewT_op{store,1}) * [0;0;0;1];
                r_k(1:3,store) = r(1:3,1);
                
                % Transformation matrix
                T = eye(4);
                T(1:3,1:3) = C_k{store,1};
                T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
                T_kpre{store,1} = T;
                
                
            else
                Tk = t(i)-t(i-1);   % Variable period at each timestep
                % Translation vector
                d_k = v_vk_vk_i(1:3,i-1) * Tk;
                r_k(1:3,store) = r_k(1:3,store-1) + C_k{store-1,1}.' * d_k;
                
                % Rotation matrix
                psi_k(1:3,store) = w_vk_vk_i(1:3,i-1) * Tk;
                C_k{store,1} = vec2rot(psi_k(1:3,store)) * C_k{store-1,1};
                
                
                % Transformation matrix
                T = eye(4);
                T(1:3,1:3) = C_k{store,1};
                T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
                T_kpre{store,1} = T;
                
                T_eps{store,1} = vec2tran(eps(store*6-5:store*6,1));
                NewT_op{store,1} =  T_eps{store,1} * T_kpre{store,1};
                T2 = NewT_op{store,1};
                C_k{store,1} = T2(1:3,1:3);
                r = inv(NewT_op{store,1}) * [0;0;0;1];
                r_k(1:3,store) = r(1:3,1);
                
                T = eye(4);
                T(1:3,1:3) = C_k{store,1};
                T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
                T_kpre{store,1} = T;
            end
        
        end
        
        
        % Input error term and F_kpre
        if i == k1
            
            e_vk(1:6,store) = tran2vec(T_kpre{store,1} / T_kpre{store,1});
            
        else
            
            e_vk(1:6,store) = tran2vec(T_kpre{store-1,1} / T_kpre{store,1});
            F_kpre{store-1,1} = tranAd(T_kpre{store,1} / T_kpre{store-1,1});
            
        end
        
        % Measurement model
        % Error term e_jk = y_jk - g(arg)
        av_points = y_k_j(:,i,:);
        av_pos = rho_i_pj_i;
        colsav = any(av_points ~= -1,1);
        num(store) = sum(colsav);
        
        if num(store) > 0
            [y_jk, g_k, G_cam, P_j ] = measModel(av_points, av_pos, DT, T_cv, T_kpre{store,1}, fu, fv, cu, cv, b);
            
            % Compute stacked error term and G_jk
            % G_jk = G_cam * DT * T_cv * point2fs(T_kpre{store,1} * P_j)
            error_k = [];
            G_instk = [];
            R_k = [];
            for j = 1:num(store)
                
                error_k(j*4 - 3:j*4, 1) = y_jk(1:4,j) - g_k(1:4,j);
                G_instk(j*4 - 3:j*4, :) = G_cam{j,1} * DT * T_cv * point2fs(T_kpre{store,1} * P_j(1:4,j));
                % Now compute and store R_K matrices
                
                if j ==1
                    R_k = R;
                else
                    R_k = blkdiag(R_k,R);
                end
                
            end
            
            e_jk{store,1} = error_k;
            G_jk{store,1} = G_instk;
            R_F{store,1} = R_k;
            
        end
        
        store = store + 1;
    end
    
    [H] = matHcompute(F_kpre, G_jk);

    % Stack e_jk errors
    E_JK = [];
    R_Mat = [];
    SizeObs = size(e_jk);
    for i = 1: SizeObs(1,1)
        
        if isempty(e_jk{i,1}) == 0
            if i ==1
                E_JK = e_jk{i,1};
                
            else 
                E_JK = [E_JK;e_jk{i,1}];
            end  
        else
            if i ==1
                E_JK = zeros(4,1);
                
            else 
                E_JK = [E_JK;zeros(4,1)];
            end
        end
        
        
        
        if isempty(R_F{i,1}) == 0
            if i == 1
                R_final = R_F{i,1};
            else
                R_final = blkdiag(R_final, R_F{i,1});
            end
        else
            if i == 1
                R_final = R;
            else
                R_final = blkdiag(R_final, R);
            end
        end
    end
    
    if ite == 1
        % Create matrix W
        W = blkdiag(Q_k,R_final);
    end
    
    % Create Error vector
    Error_vec = [e_vk(:); E_JK];
    % Setup solver
    A = (transpose(H) / W)  * H;
    B = (transpose(H) /W) * Error_vec;
    %eps = linsolve(A,B);
    AA = chol(A);
    eps = AA\(AA'\B);
    
end