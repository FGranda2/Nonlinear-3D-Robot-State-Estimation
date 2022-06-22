function [y_jk, g_k, G_cam, P_j] = measModel(av_points, av_pos, DT, T_cv, T_op, fu, fv, cu, cv, b, M)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    colsav = any(av_points ~= -1,1);
    av_points = av_points(:, colsav);
    av_pos = av_pos(:,colsav);

    sizeTotal = size(av_pos);
    sizeTotal = sizeTotal(1,2);
    %av_pos = [av_pos; ones(1,sizeTotal)];
    
    r_vi =inv(T_op) * [0;0;0;1];
    r_vi = DT*r_vi;
    C_vi = T_op(1:3,1:3);
    
    r_cv =inv(T_cv) * [0;0;0;1];
    r_cv = DT*r_cv;
    C_cv = T_cv(1:3,1:3);

    for j =  1:sizeTotal
        % Available argument for g(.)
        %av_args(1:3, j) = DT * T_cv * T_op * av_pos(1:4, j);
        av_args(1:3, j) = C_cv * (C_vi * (av_pos(1:3, j)-r_vi)-r_cv);
        
        x = av_args(1, j);
        y = av_args(2, j);
        z = av_args(3, j);

        % Image Pixel values for each available point
        g_k(1:4, j) = (1/z) * [fu*x;fv*y;fu*(x-b);fv*y]+[cu;cv;cu;cv];

        % Compute camera model Jacobian
        %G_cam{j,1} = [fu/z,0,-fu*x/z^2;0,fv/z,-fv*y/z^2;...
        %    fu/z,0,fu*(b-x)/z^2;0,fv/z,-fv*y/z^2];

    end
  
    for j =  1:sizeTotal
        G_cam{j,1} = camJac(M,[av_args(1:3, j);1]);
    end
    
    y_jk = av_points;
    av_pos = [av_pos; ones(1,sizeTotal)];
    P_j = av_pos;
    
end

