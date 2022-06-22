function [y_jk, g_k, G_cam, P_j] = measModel2(av_points, ...
    av_pos, DT, T_cv, T_op, M)

    colsav = any(av_points ~= -1,1);
    av_points = av_points(:, colsav);
    av_pos = av_pos(:,colsav);

    sizeTotal = size(av_pos);
    sizeTotal = sizeTotal(1,2);
    av_pos = [av_pos; ones(1,sizeTotal)];

    for j =  1:sizeTotal
        
        g_k(1:4,j) = cam(M,T_cv*T_op*av_pos(1:4, j));
        G_cam{j,1} = camJac(M,T_cv*T_op*av_pos(1:4, j));
    end
    
    y_jk = av_points;
    P_j = av_pos;
    
end