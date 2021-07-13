function [G_final, nIter_final, S_final, objhistory_final, nIteration] = DiMMA_function(gnd, m, R, G, L, nClass, ntype, options, Q, T);
 

if ~isfield(options,'error')
    options.error = 1e-6;
end
if ~isfield(options, 'maxIter')
    options.maxIter = [];
end
 
if ~isfield(options,'nRepeat')
    options.nRepeat = 10;
end
 
if ~isfield(options,'minIter')
    options.minIter = 30;
end
 
if ~isfield(options,'meanFitRatio')
    options.meanFitRatio = 0.1;
end
 
if ~isfield(options,'alpha')
    options.alpha = 100;
end

if ~isfield(options,'Optimization')
    options.Optimization = 'Multiplicative';
end
 
if ~exist('G1','var') %k
    G1 = [];
    G2 = [];
    S = [];
end
 
differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;
 
alpha = options.alpha;
beta = options.beta; %kk910
delta = options.delta;
G_final = cell(m,1);
for i = 1:m
    G_final{i,1} = G{i,1};
end 
Norm = 2;
NormV = 0;
 
%construct L based on u
 
u = ones(1, 1);
u = [1];
q = size(u,2);%number of W
sumu = sum(u);
for i = 1:size(u,2)
    u(:,i) = u(:,i)/sumu;
end
 
u_final = u;
Lnorm = cell(m,1);
for i = 1:m
    Lnorm{i} = zeros(size(L{i},1));
    for j = 1:size(u_final,2)
        Lnorm{i} = Lnorm{i} + u_final(1,j)*L{i};
    end
end
% connect T and L
for i = 1:m
    Lnorm{i} = Lnorm{i} + delta*T{i};
end
% end connect T and L
Lnorm1 = cell(m,1);
Lnorm0 = cell(m,1);
for i = 1:m
    Lnorm1{i} = zeros(size(Lnorm{i},1),size(Lnorm{i},2));
    Lnorm0{i} = zeros(size(Lnorm{i},1),size(Lnorm{i},2));
    Lnorm1{i} = (abs(Lnorm{i}) + Lnorm{i})*0.5;
    Lnorm0{i} = (abs(Lnorm{i}) - Lnorm{i})*0.5;   
end
 
%end chuan bi L
 
selectInit = 1;
Rd_t = R{1,2};
nSmp = size(Rd_t,1);
mFea = size(Rd_t,2);
% Initialize the data and feature matrices 
G = initializeMTRD(R, nClass,m);

% initialise matrix S
S = cell(m,m);
for i=1:m
    for j=1:m
        S{i,j}=zeros(ntype(i),ntype(j));
    end
end
tryNo = 0;
nIter = 0;
nIteration = 0;
while tryNo < nRepeat   
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror)
       nIteration = nIteration + 1;
        
         % ===================== update S ========================
         S = updateS(S,G,R,m);
         
         % ===================== update G ========================   
         G = updateG(G,S,R,m, ntype, nClass, alpha, beta, delta, Lnorm0, Lnorm1, Q, T) ;    
         G = l1_norm(G,m);
         
    
        
        nIter = nIter + 1;
        
        
%         When U, V run nIter times
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(R,S, G, Lnorm, options,m, Q, T); %kk
                maxErr = 0;
            else
                if isempty(maxIter)
                    newobj = CalculateObj(R,S, G, Lnorm, options,m, Q, T);
                    objhistory = [objhistory newobj]; %#ok<AGROW>
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge') && options.Converge
                        newobj = CalculateObj(R,S, G, Lnorm, options,m, Q, T);
                        objhistory = [objhistory newobj]; %#ok<AGROW>
                        meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;%k
                        maxErr = (meanFit-newobj)/meanFit;% k
                        
                    end
%                     maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = 0;
                        end
                    end
                end
            end
        end
    end
    
%     When nIter achieves minIter, run the following code segment
    if tryNo == 1
        for i = 1:m
            G_final{i,1} = G{i,1};
        end
%         for i=1:m
%             for j=1:m
%                 if j > i
%                     S_final{i,j} = S{i,j};
%                 end
%             end
%         end
        S_final = S;
        nIter_final = nIter;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end)
            for i = 1:m
                G_final{i,1} = G{i,1};
            end
%             for i=1:m
%                 for j=1:m
%                     if j > i
%                         S_final{i,j} = S{i,j};
%                     end
%                 end
%             end
            S_final = S;
           nIter_final = nIter;
           objhistory_final = objhistory;
       end
    end
 
    if selectInit
        if tryNo < nRepeat
            %re-start
            G = initializeMTRD(R, nClass,m);
% % % % % % %             rand('twister',5489);
% % % % % % %             labeld = litekmeans(Rd_t, nClass);  % doc
% % % % % % %             rand('twister',5489);
% % % % % % %             labelt = litekmeans(Rd_t', nClass);  % term
% % % % % % %             % labelc = litekmeans(Rd_c', nClass);  %concept
% % % % % % %             for i = 1:nClass  % nCla == mCla
% % % % % % %                G{1}(labeld ==i, i) = 1;
% % % % % % %                G{2}(labelt ==i, i) = 1;
% % % % % % %             %    G{3}(labelc ==i, i) = 1;
% % % % % % %             end
% % % % % % %             for i = 1:m
% % % % % % %                 G{i} = G{i}+0.2;
% % % % % % %             end
%             end re-start
 
            nIter = 0;
        else
            tryNo = tryNo - 1;
            nIter = minIter+1;
            selectInit = 0;
            for i = 1:m
                G{i,1} = G_final{i,1};
            end
%             for i=1:m
%                 for j=1:m
%                     if j > i
%                         S{i,j} = S_final{i,j};
%                     end
%                 end
%             end
            S = S_final;
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end
end
 
[G{2},G{1}] = NormalizeUV(G_final{2,1}, G_final{1,1}, NormV, Norm);
 
 
%==========================================================================
function objhistory_final = CalculateObj(R,S, G, Lnorm, options, m, Q, T)

        alpha = options.alpha;
        beta = options.beta;
        delta = options.delta;
        
        obj_NMF = 0;
        obj_Poss = 0;
        for i=1:m
            for j=1:m
                if j~=i
                    obj_NMF = obj_NMF + norm(R{i,j} - G{i}*S{i,j}*G{j}','fro');
                    obj_Poss = obj_Poss + trace(G{i}'*Q{i,j}*G{j});
                end
            end
        end
        obj_Lap = 0; 
        obj_norm = 0;
        for i = 1:m
            obj_Lap = obj_Lap + trace(G{i}'*Lnorm{i}*G{i});
%             obj_norm = obj_norm + norm(G{i},'fro');
%             obj_norm = obj_norm + trace(G{i}'*T{i}*G{i});
        end
        objhistory_final = obj_NMF + 2*alpha*obj_Lap - 2*delta*obj_Poss; % + beta*obj_norm;     
 
    
function G = l1_norm(G,m)
for p = 1:m      
    for i = 1:size(G{p},1)
                if sum(G{p}(i,:))~= 0
                    G{p}(i,:) = G{p}(i,:)/sum(G{p}(i,:));
                else
                    for j = 1:size(G{p},2)
                        G{p}(i,j) = 1/(size(G{p},2));
                    end
                end
    end
end 
function [U, V] = NormalizeUV(U, V, NormV, Norm)
    K = size(U,2);
    if Norm == 2
        if NormV
            norms = max(1e-15,sqrt(sum(V.^2,1)))';
            V = V*spdiags(norms.^-1,0,K,K);
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sqrt(sum(U.^2,1)))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = V*spdiags(norms,0,K,K);
        end
    else
        if NormV
            norms = max(1e-15,sum(abs(V),1))';
            V = V*spdiags(norms.^-1,0,K,K);
            U = U*spdiags(norms,0,K,K);
        else
            norms = max(1e-15,sum(abs(U),1))';
            U = U*spdiags(norms.^-1,0,K,K);
            V = V*spdiags(norms,0,K,K);
        end
    end
function S = updateS(S,G,R,m)
    for i=1:m
          for j=1:m
                if j > i
%                     VV1 = (G{i}'*G{i})^-1;
%                     VV2 = (G{j}'*G{j})^-1;
                    VV1 = pinv(G{i}'*G{i});
                    VV2 = pinv(G{j}'*G{j});
                    S{i,j} = (((VV1*G{i}')*R{i,j})*G{j})*VV2;
                    S{j,i} = S{i,j}';
                end
           end
    end
function G = updateG(G,S,R,m, ntype, nClass, alpha, beta, delta, Lnorm0, Lnorm1, Q, T)        
        A = cell(m,1);
        B = cell(m,1);
        Aplus = cell(m,1); 
        Aminus = cell(m,1);
        Bplus = cell(m,1);
        Bminus = cell(m,1);
%         khoi tao A{i}
        for i = 1:m
            A{i} = zeros(ntype(i),nClass);
        end
%         kiem tra Lplus and minus
%         tinh A{i}
        for i=1:m
            for j=1:m
                if j~=i
                    A{i} = A{i} + R{i,j}*G{j}*S{i,j}' + delta*Q{i,j}*G{j};
                end
            end
        end
%         Tinh A+, A-
        Aplus = A;
        Aminus = A;
        for i = 1:m 
            Aplus{i} = (abs(A{i}) + A{i})*0.5;
            Aminus{i} = (abs(A{i}) - A{i})*0.5;
        end
        for i = 1:m
            B{i} = zeros(nClass,nClass);
        end
%         Tinh B{i}
        for i=1:m
            for j=1:m
                if j~=i
                    B{i} = B{i} + S{i,j}*G{j}'*G{j}*S{i,j}';
                end
            end
        end
%         Tinh B+, B-
        
        Bplus = B;
        Bminus = B;
        for i = 1:m 
            Bplus{i} = (abs(B{i}) + B{i})*0.5;
            Bminus{i} = (abs(B{i}) - B{i})*0.5;
        end
        
        for p = 1:m
            tempup = alpha*Lnorm0{p}*G{p} + Aplus{p} + G{p}*Bminus{p}; % + beta*T{p}*G{p};
            tempun = alpha*Lnorm1{p}*G{p} + Aminus{p} + G{p}*Bplus{p}; % + beta*T{p}*G{p};
            for j = 1:size(G{p},2)
                for i = 1:size(G{p},1)
                    if tempun(i,j)~=0
                        G{p}(i,j) = G{p}(i,j)*(tempup(i,j)/tempun(i,j))^(0.5);
                        
                    else
                        G{p}(i,j) = 0;
                    end
                end
            end
        end
