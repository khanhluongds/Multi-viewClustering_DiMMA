% cho Matran R
% R = cell(3,3);
% m = 3;
% ntype = zeros(m,1);
% ntype(1) = 3; ntype(2) = 3; ntype(3) = 4;
% R{1,2} = [1 2 1; 2 2 0; 4 0 3];
% R{1,3} = [0 1 3 1; 1 2 3 0; 4 1 1 2];
% R{2,3} = [2 1 0 3; 4 0 0 1; 1 2 3 2];
% R{2,1} = R{1,2}'; 
% R{3,1} = R{1,3}';
% R{3,2} = R{2,3}'; 
% options.p = 1;
% for i = 1:m
%     for j = 1:m
%         temp = zeros(ntype(i), ntype(j));
%         R{i,j} = temp;
%     end
% end 
function [T,Q] = constructTQ(R, options, ntype, m)
% khoi tao Rnew
for i = 1:m
    for j = 1:m
        temp = zeros(ntype(i), ntype(j));
        Rnew{i,j} = temp;
    end
end 
% construcR de co Rnew


for i = 1:m
    for j=1:m
        if i~=j
            Rnew{i,j} = constructR(R{i,j}, options);
            Rnew{i,j} = full(Rnew{i,j});
        end
    end
end

% construct T

Rnew_sum1 = cell(m,m);
Rnew_sum2 = cell(m,m);
for i = 1:m
    for j = 1:m
%         if i~=j
            Rnew_sum1{i,j} = sum(Rnew{i,j},2);
            Rnew_sum2{i,j} = sum(Rnew{i,j},1);
            Rnew_sum1{i,j} = diag(Rnew_sum1{i,j});
            Rnew_sum2{i,j} = diag(Rnew_sum2{i,j});
%         end
    end
end

for i = 1:m
    temp = zeros(ntype(i), ntype(i));
%     temp = [];
    for j = 1:m
        if i~=j
            temp = temp + Rnew_sum1{i,j} + Rnew_sum2{j,i}';
        end
    end 
    T{i} = temp;
end

% construct Q
% khoi tao Q
for i = 1:m
    for j = 1:m
        temp = zeros(ntype(i), ntype(j));
        Q{i,j} = temp;
    end
end 
for i = 1:m
%     temp = zeros(ntype(i), ntype(i));
%     temp = [];
    for j = 1:m
        if i<j
            Q{i,j} = Q{i,j} + Rnew{i,j} + Rnew{j,i}';
        end
    end 
%     T{i} = temp;
end




