%initialize multi-type, m is the number of types
function [G] = initializeMTRD(R, nClass, m)
labeltemp = cell(m,1);
rand('twister',5489);
labeltemp{1,1} = litekmeans(R{1,2}, nClass); 
for i = 2:m
    rand('twister',5489);
    labeltemp{i,1} = litekmeans(R{i,1}, nClass); 
end
for h = 1:m
    for i = 1:nClass  % nCla == mCla
        G{h,1}(labeltemp{h,1} ==i, i) = 1;
    end
end
for i = 1:m
    G{i,1} = G{i,1}+0.2;
end    
