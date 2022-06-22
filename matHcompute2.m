function [H_final] = matHcompute2(F_kpre, G_k)

Size = size(F_kpre);
Size = Size(1,1);
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
condition = 0;
adjust = 0;

for i = 1:Size
   if isempty(G_k{i,1}) == 0
       
        if i == 1
           G_mat =  G_k{i,1};
           condition = 0;
           
        elseif isempty(G_mat) == 1
           sizeGK = size(G_k{i,1});
           G_mat(1:sizeGK(1,1),(6*adjust) + 1 :(6*adjust)+6) = G_k{i,1};
           condition = 0;
        else
           G_mat = blkdiag(G_mat, G_k{i,1});
        end
   else
        if i == 1 || condition ==1
            G_mat = [];
            adjust = adjust + 1;
            condition = 1;
        else
            newSize = size(G_mat);
            G_mat = [G_mat, zeros(newSize(1,1),6)];
        end
    end 
end
H_final = [Total;G_mat];
end