clear;
addpath(genpath(pwd));
datasetname = 'Movie'
load(datasetname);

m = nviews + 1;

% File name
nClass = length(unique(gnd)); %document cluster number

%Normalize each data vector to have L2-norm equal to 1  
 
%Clustering in the original space using K-means
rand('twister',5489);
label = litekmeans(fea{1},nClass,'Replicates',20);
NMI_Kmeans = MutualInfo(gnd,label);
disp(['Clustering in the original space. NMI: ',num2str(NMI_Kmeans)]);

%DiMMA learning
options = [];
options.WeightMode = 'Binary';  
options.maxIter = 100;

% initialise R is a matrix of mxm matrix storing inter relationships
% between object types
R = cell(m,m);
% number of elements (objects) in each object types
ntype = zeros(m,1); %ntype{1} stores the number of objects in object type 1
% assign the number of object to ntype based on dataset
ntype(1) = size(fea{1},1);
for i = 2:m
    ntype(i) = size(fea{i-1},2);
end

for i = 1:m
    for j = 1:m
        temp = zeros(ntype(i),ntype(j));
        R{i,j} = temp;
    end 
end 
for i = 2:m
    R{1,i} = fea{i-1}; 
    R{i,1} = R{1,i}';
end

% initialise and assign intra-type relationship
W = cell(m,1);
tempW = cell(m,1);
alpha = 1;
options.alpha = alpha;
options.normW = 1;
for i = 1:m
    tempW{i} = [];
    for j = 1:m
        if (i~=j) && ~isequal(R{i,j},zeros(ntype(i),ntype(j)))
            tempW{i} = [tempW{i}, R{i,j}];
        end
    end
    W{i} = constructW(tempW{i},options);
    L{i} = constructL(W{i}, alpha, options);
end
 
G = cell(m,1);
% Initialize for G
label = cell(m,1);
rand('twister',5489);
label{1} = litekmeans(R{1,2}, nClass); 
for i = 2:m
    label{i} = litekmeans(R{i,1}, nClass); 
end
for h = 1:m
    for i = 1:nClass  
        G{h,1}(label{h} ==i, i) = 1;
    end
end
for i = 1:m
    G{i,1} = G{i,1}+0.2;
end
%% Using both inter and intra manifold leanring
options.alpha = 1;
options.beta = 1;
options.delta = 0.1*options.alpha;
options.p = 5;
Rnew = cell(m,m);    
[T, Q] = constructTQ(R, options, ntype, m);
%%
tic
[G_final, nIter_final, S_final, objhistory_final, nIteration] = DiMMA_function(gnd, m, R, G, L, nClass, ntype, options, Q, T);
time = toc;

rand('twister',5489);
label = litekmeans(G_final{1,1},nClass,'Replicates',20);
MIhat = MutualInfo(gnd,label); 
d = [gnd, label];
fscore = FScr(d);
gnd1 = gnd;
labelnew = bestMap(gnd1, label);
AC = length(find(gnd == labelnew))/length(gnd);

disp(['p =--------------------------------- ',num2str(options.p)]);
disp(['alpha = ',num2str(options.alpha),', delta = ',num2str(options.delta),', beta = ',num2str(options.beta)]);
disp(['Final MIhat = ',num2str(MIhat)]);
disp(['Final F_score = ',num2str(fscore)]);
disp(['Final Accuracy = ',num2str(AC)]);
disp(['Running time (seconds) = ',num2str(time)]);


