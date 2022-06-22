clear 
close

% First load the dataset variables
load dataset3.mat

% Now start by stating some definitions
% Motion model noise variance
Q = diag([v_var; w_var]);    % Diagonal 6x6 matrix

% Measurement model noise variance
R = diag(y_var);    % Diagonal 4x4 matrix

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
kinit = 1215;
kend = 1714;
ksize = 10;
r_init = r_i_vk_i(:,kinit);
theta_init = -theta_vk_i(:,kinit);
finalStore = 1;

for window = kinit:kend
    
    k1 = window;
    k2 = window + ksize;

    [T_FINAL{finalStore,1},A{finalStore,1}] = windowGN(k1,k2,r_init,...
        theta_init,t,Q_k,Q,y_k_j,...
        rho_i_pj_i,DT, T_cv,M,R,v_vk_vk_i,w_vk_vk_i);
    
    T = T_FINAL{finalStore,1};
    C = T(1:3,1:3);
    theta_init = rot2vec(C);
    r = inv(T_FINAL{finalStore,1}) * [0;0;0;1];
    r_init = r(1:3,1);
    r_k(1:3,finalStore) = r(1:3,1);
    finalStore = finalStore+1;
    
end