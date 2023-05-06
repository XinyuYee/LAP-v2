%% LAP: Local affine preserving matching
clear;
close all
numNeigh1 = 25; %20
K = 10; % 10
alpha = 0.5;
lambda1 = 0.7; 
final_num = alpha *  nchoosek(K,3);  %60;  %10
%% 
fn_l = '.\127_r.JPG';
fn_r = '.\127_l.JPG';
Ia = imread(fn_l);
Ib = imread(fn_r);
load  .\127.mat;
Ia = imresize(Ia, [size(Ib,1), size(Ib,2)]);
if size(Ia,3) == 1
    Ia = repmat(Ia,[1,1,3]);
    Ib = repmat(Ib,[1,1,3]);
end
%% 
X_ori = X; Y_ori = Y;
idxUnique = removeRepeat(X,Y);
X=X_ori(idxUnique,:);
Y=Y_ori(idxUnique,:);

tic;
[N,~] = size(X);
Xt = X';Yt = Y';
vec=Yt-Xt;
d2=vec(1,:).^2+vec(2,:).^2;
vx = vec(1, :); vy = vec(2, :);

kdtreeX = vl_kdtreebuild(Xt);
[neighborX, ~] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', numNeigh1+1) ;
index = neighborX(2:numNeigh1+1, :);
vxi = vx(index); vyi = vy(index);
d2i = d2(index);

cos_sita = (vxi.*repmat(vx,size(vxi,1),1) + vyi.*repmat(vy,size(vyi,1),1)) ./ sqrt(d2i.*repmat(d2,size(d2i,1),1));
ratio = min(d2i, repmat(d2,size(d2i,1),1)) ./ max(d2i, repmat(d2,size(d2i,1),1));
%c2i = cos_sita.*ratio;
c2i = 0.5*cos_sita + ratio;
[~,I] = sort(c2i,1,'descend'); 

row = 0:N-1;
I = I + repmat(row.*numNeigh1,numNeigh1,1);
Index = index(I);
sort_neibor = Index(1:K,:);
%sort_neibor = index(1:K,:);
%% STEP 2
x1 = X(:,1); x2 = X(:,2); 
y1 = Y(:,1); y2 = Y(:,2); 
pum = nchoosek(1:K,3); pum_1 =  pum(:,1); pum_2 =  pum(:,2); pum_3 =  pum(:,3);
pum_size = size(pum_1,1);
sort_neibor_1 = sort_neibor(repmat(pum_1,1,N)+ repmat(row.*K,pum_size,1));
sort_neibor_2 = sort_neibor(repmat(pum_2,1,N)+ repmat(row.*K,pum_size,1));
sort_neibor_3 = sort_neibor(repmat(pum_3,1,N)+ repmat(row.*K,pum_size,1));

% pum_rat = nchoosek(1:pum_size,2); pum_rat_1 = pum_rat(:,1); pum_rat_2 = pum_rat(:,2); 
% pum_rat_size = size(pum_rat,1);
% s1_index = repmat(pum_rat_1,1,N)+ repmat(row.*pum_size,pum_rat_size,1);
% s2_index = repmat(pum_rat_2,1,N)+ repmat(row.*pum_size,pum_rat_size,1);
%% 计算原图像的比率向量
%第一个点
M11 = x1(sort_neibor_1) - repmat(x1',pum_size,1);
M12 = x2(sort_neibor_1) - repmat(x2',pum_size,1);
%第二个点
M21 = x1(sort_neibor_2) - repmat(x1',pum_size,1);
M22 = x2(sort_neibor_2) - repmat(x2',pum_size,1);
%第三个点
M31 = x1(sort_neibor_3) - repmat(x1',pum_size,1);
M32 = x2(sort_neibor_3) - repmat(x2',pum_size,1);
%计算面积
S_11 = abs(M11.*M22- M12.*M21);
S_12 = abs(M31.*M22- M32.*M21);
S_13 = abs(M11.*M32- M12.*M31);
%Sx_rat = (1 - exp(-S_11./S_12)) + (1 - exp(-S_11./S_13)) + (1 - exp(-S_12./S_13));

%% 计算目标图像的比率向量
%第一个点
M11 = y1(sort_neibor_1) - repmat(y1',pum_size,1);
M12 = y2(sort_neibor_1) - repmat(y2',pum_size,1);
%第二个点
M21 = y1(sort_neibor_2) - repmat(y1',pum_size,1);
M22 = y2(sort_neibor_2) - repmat(y2',pum_size,1);
%第三个点
M31 = y1(sort_neibor_3) - repmat(y1',pum_size,1);
M32 = y2(sort_neibor_3) - repmat(y2',pum_size,1);
%计算面积
S_21 = abs(M11.*M22- M12.*M21);
S_22 = abs(M31.*M22- M32.*M21);
S_23 = abs(M11.*M32- M12.*M31);
%Sy_rat =(1 - exp(-S_21./S_22)) + (1 - exp(-S_21./S_23)) + (1 - exp(-S_22./S_23));

S_rat = (1-exp(-abs(S_11./S_12 - S_21./S_22))) + (1-exp(-abs(S_11./S_13 - S_21./S_23))) + (1-exp(-abs(S_12./S_13 - S_22./S_23)));
S_rat = abs(S_rat);
%S_rat = abs(Sy_rat - Sx_rat);
S_rat_sort = sort(S_rat,1); 
S_rat_sort_no = S_rat_sort(1:round(final_num),1:N);
%S_rat_sort_final = 1-exp(-S_rat_sort_no);
S_rat_sort_final = S_rat_sort_no;
differ = sum(S_rat_sort_final,1)/final_num;
differ_X = differ;% 1 - exp(-differ/0.85); %

%%
X_temp = X;
Y_temp = Y;
Y = X_temp;
X = Y_temp;

Xt = X';Yt = Y';
vec=Yt-Xt;
d2=vec(1,:).^2+vec(2,:).^2;
vx = vec(1, :); vy = vec(2, :);

kdtreeX = vl_kdtreebuild(Xt);
[neighborX, ~] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', numNeigh1+1) ;
index = neighborX(2:numNeigh1+1, :);
vxi = vx(index); vyi = vy(index);
d2i = d2(index);

cos_sita = (vxi.*repmat(vx,size(vxi,1),1) + vyi.*repmat(vy,size(vyi,1),1)) ./ sqrt(d2i.*repmat(d2,size(d2i,1),1));
ratio = min(d2i, repmat(d2,size(d2i,1),1)) ./ max(d2i, repmat(d2,size(d2i,1),1));
%c2i = cos_sita.*ratio;
c2i = 0.5*cos_sita + ratio;
[~,I] = sort(c2i,1,'descend'); 

row = 0:N-1;
I = I + repmat(row.*numNeigh1,numNeigh1,1);
Index = index(I);
sort_neibor = Index(1:K,:);
% sort_neibor = index(1:K,:);
%% STEP 2
x1 = X(:,1); x2 = X(:,2); 
y1 = Y(:,1); y2 = Y(:,2); 
pum = nchoosek(1:K,3); pum_1 =  pum(:,1); pum_2 =  pum(:,2); pum_3 =  pum(:,3);
pum_size = size(pum_1,1);
sort_neibor_1 = sort_neibor(repmat(pum_1,1,N)+ repmat(row.*K,pum_size,1));
sort_neibor_2 = sort_neibor(repmat(pum_2,1,N)+ repmat(row.*K,pum_size,1));
sort_neibor_3 = sort_neibor(repmat(pum_3,1,N)+ repmat(row.*K,pum_size,1));

% pum_rat = nchoosek(1:pum_size,2); pum_rat_1 = pum_rat(:,1); pum_rat_2 = pum_rat(:,2); 
% pum_rat_size = size(pum_rat,1);
% s1_index = repmat(pum_rat_1,1,N)+ repmat(row.*pum_size,pum_rat_size,1);
% s2_index = repmat(pum_rat_2,1,N)+ repmat(row.*pum_size,pum_rat_size,1);
%% 计算原图像的比率向量
%第一个点
M11 = x1(sort_neibor_1) - repmat(x1',pum_size,1);
M12 = x2(sort_neibor_1) - repmat(x2',pum_size,1);
%第二个点
M21 = x1(sort_neibor_2) - repmat(x1',pum_size,1);
M22 = x2(sort_neibor_2) - repmat(x2',pum_size,1);
%第三个点
M31 = x1(sort_neibor_3) - repmat(x1',pum_size,1);
M32 = x2(sort_neibor_3) - repmat(x2',pum_size,1);
%计算面积
S_11 = abs(M11.*M22- M12.*M21);
S_12 = abs(M31.*M22- M32.*M21);
S_13 = abs(M11.*M32- M12.*M31);
%Sx_rat = (1 - exp(-S_11./S_12)) + (1 - exp(-S_11./S_13)) + (1 - exp(-S_12./S_13));

%% 计算目标图像的比率向量
%第一个点
M11 = y1(sort_neibor_1) - repmat(y1',pum_size,1);
M12 = y2(sort_neibor_1) - repmat(y2',pum_size,1);
%第二个点
M21 = y1(sort_neibor_2) - repmat(y1',pum_size,1);
M22 = y2(sort_neibor_2) - repmat(y2',pum_size,1);
%第三个点
M31 = y1(sort_neibor_3) - repmat(y1',pum_size,1);
M32 = y2(sort_neibor_3) - repmat(y2',pum_size,1);
%计算面积
S_21 = abs(M11.*M22- M12.*M21);
S_22 = abs(M31.*M22- M32.*M21);
S_23 = abs(M11.*M32- M12.*M31);
%Sy_rat =(1 - exp(-S_21./S_22)) + (1 - exp(-S_21./S_23)) + (1 - exp(-S_22./S_23));

S_rat = (1-exp(-abs(S_11./S_12 - S_21./S_22))) + (1-exp(-abs(S_11./S_13 - S_21./S_23))) + (1-exp(-abs(S_12./S_13 - S_22./S_23)));
S_rat = abs(S_rat);
%S_rat = abs(Sy_rat - Sx_rat);
S_rat_sort = sort(S_rat,1); 
S_rat_sort_no = S_rat_sort(1:round(final_num),1:N);
%S_rat_sort_final = 1-exp(-S_rat_sort_no);
S_rat_sort_final = S_rat_sort_no;
differ = sum(S_rat_sort_final,1)/final_num;
differ_Y = differ;% 1 - exp(-differ/0.85); %
%%
differ = (differ_X+differ_Y)/2;
p1 = differ <= lambda1.*ones(1,N);
ind_pro = find( p1 == 1 );
ind = ind_pro;
toc
[inlier_num,inlierRate,precision_rate,Recall_rate,F1_score] = evaluatePR(X,CorrectIndex,idxUnique(ind))
%CorrectIndex = intersect(CorrectIndex,idxUnique);


