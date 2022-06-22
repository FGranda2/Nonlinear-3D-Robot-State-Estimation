function [H_final] = matHcompute(F_kpre, G_k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Size = size(F_kpre);
Size = Size(1,1);
%oneMat = ones(6,6);
oneMat = eye(6);

for i = 1:Size
   
    if i == 1
        H = blkdiag(oneMat,oneMat);
        H(7:12,1:6) = -F_kpre{i,1};
        Total = H;
    else
        Total = blkdiag(Total,oneMat);
        Total((6*i)+1:6*(i+1),(6*i)-5:6*i) = -F_kpre{i,1};
    end
    
end

Size = size(G_k);
Size = Size(1,1);

for i = 1:Size
   if isempty(G_k{i,1}) == 0
       
        if i == 1
           G_mat =  G_k{i,1};
        else
           G_mat = blkdiag(G_mat, G_k{i,1});
        end
   else
        if i == 1
            G_mat = zeros(4,6);
        else
            newSize = size(G_mat);
            G_mat = [G_mat, zeros(newSize(1,1),6)];
        end
    end 
end
H_final = [Total;G_mat];
end

