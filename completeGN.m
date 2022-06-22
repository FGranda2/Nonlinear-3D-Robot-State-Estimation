clear 
close

% First load the dataset variables
load dataset3.mat

% Now start by stating some definitions
% Motion model noise variance
Q = diag([v_var; w_var]);    % Diagonal 6x6 matrix

% Measurement model noise variance
R = diag(y_var);    % Diagonal 4x4 matrix

%Variance Initial Condition
Q_k = eye(6,6)*0.5;

% Set transformation matrix T_cv
T_cv = eye(4);
T_cv(1:3,1:3) = C_c_v;
T_cv(1:3,4) = -C_c_v * rho_v_c_v;

% Camera matrix
M = [fu,0,cu,0;...
    0,fv,cv,0;...
    fu,0,cu,-fu*b;...
    0,fv,cv,0];

% D traspose matrix creation
DT = [1,0,0,0;0,1,0,0;0,0,1,0];

% Limits
k1 = 1215;
k2 = 1714;
store = 1;
for i = k1:k2
    % Motion model equations to estimate T_kpre
    if i == k1
        Tk = 0;
        % Translation vector
        r_k0 = r_i_vk_i(:,i);
        % Rotation matrix
        psi_k0 = -theta_vk_i(:,i);
        C_k0 = vec2rot(psi_k0);

        r_k(1:3,store) = r_k0;

        % Rotation matrix
        C_k{store,1} = C_k0;

        % Transformation matrix
        T = eye(4);
        T(1:3,1:3) = C_k{store,1};
        T(1:3,4) = -C_k{store,1} * r_k(1:3,store);
        T_kpre{store,1} = T;
        T_0 = T;

    else
        Tk = t(i)-t(i-1);   % Variable period at each timestep
        
        Win{store,1} = [v_vk_vk_i(1:3,i); w_vk_vk_i(1:3,i)];
        T = vec2tran(Tk * -Win{store,1}) * T_kpre{store-1,1};
        T_kpre{store,1} = T;
        r = inv(T_kpre{store,1}) * [0;0;0;1];
        r_k(1:3,store) = r(1:3,1);
        
        
        % W matrix append
        Q_k = blkdiag(Q_k, (Tk^2) * Q);
    end
    store=store+1;
end

T_op = T_kpre;
% Start GN update
stop = 1;
ite = 1;
while stop == 1
    store = 1;
    for i = k1:k2
       if i == k1
          
           Tk = 0;
           % Error
           e_vk(1:6,store) = tran2vec(T_0 * inv(T_op{store,1}));
           
       else
           % Get evaluation of T_opk-1
           % Translation vector
           Tk = t(i)-t(i-1);
           T = [];
           T = vec2tran(Tk * -Win{store,1}) * T_op{store-1,1};
           
           % Error values
           e_vk(1:6,store) = tran2vec(T * inv(T_op{store,1}));
           F_kpre{store-1,1} = tranAd(T_op{store,1} * ...
               inv(T_op{store-1,1}));
       end
       % Measurement model
       % Error term e_jk = y_jk - g(arg)
       av_points = y_k_j(:,i,:);
       av_pos = rho_i_pj_i;
       colsav = any(av_points ~= -1,1);
       num(store) = sum(colsav);
       if num(store) > 0
           [y_jk, g_k, G_cam, P_j ] = measModel2(av_points, ...
               av_pos, DT, T_cv, T_op{store,1}, M);
           % Compute stacked error term and G_jk
           % G_jk = G_cam * T_cv * point2fs(T_kpre{store,1} * P_j)
           error_k = [];
           G_instk = [];
           R_k = [];
           
           for j = 1:num(store)
                
               error_k(j*4 - 3:j*4, 1) = y_jk(1:4,j) - g_k(1:4,j);
               G_instk(j*4 - 3:j*4, :) = G_cam{j,1} * ...
                   T_cv * point2fs(T_op{store,1} * P_j(1:4,j));
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
       else
           e_jk{store,1} = [];
           G_jk{store,1} = [];
           R_F{store,1} = [];
           
       end    
       store=store+1;
    end
    [H] = matHcompute2(F_kpre, G_jk);
    % Stack e_jk errors
    E_JK = [];
    R_final = [];
    SizeObs = size(e_jk);
    condition1 = 0;
    condition2 = 0;
    for i = 1: SizeObs(1,1)
        
        if isempty(e_jk{i,1}) == 0
            if i ==1 || isempty(E_JK) == 1
                E_JK = e_jk{i,1};
                condition1 = 0;

            else 
                E_JK = [E_JK;e_jk{i,1}];
            end
        else
            if i ==1 || condition1 == 1
                E_JK = [];
                condition1 =1;
            end
        end
        
        if isempty(R_F{i,1}) == 0
            if ite == 1 
                if i == 1 || isempty(R_final) == 1
                    R_final = R_F{i,1};
                    condition2 = 0;
                else
                    R_final = blkdiag(R_final, R_F{i,1});
                end
            end
        else
            if ite ==1
                if i == 1 || condition2 ==1
                    R_final = [];
                    condition2 = 1;
                    
                end
            end
        end
    end
    
    if ite == 1
        % Create matrix W
        W = blkdiag(Q_k,R_final);
        W_inv = inv(W);
    end
    % Create Error vector
    Error_vec = [e_vk(:); E_JK];
    % Setup solver
    
    A = H.' * W_inv * H;
    B = H.' * W_inv * Error_vec;
    AA = chol(A);
    eps = AA\(AA'\B);
    
    %Update T_op with eps
    for i = 1:store-1
        T_eps{i,1} = vec2tran(real(0.5*eps(i*6-5:i*6,1)));
        T_op{i,1} = T_eps{i,1} * T_op{i,1};
        r = inv(T_op{i,1}) * [0;0;0;1];
        r_k(1:3,i) = r(1:3,1);
    end
    
    if norm(eps) < 1e-4
        stop =0;
    end
    ite = ite+1;
end