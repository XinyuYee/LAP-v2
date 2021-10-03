function [inlier_num,inlierRate,precision_rate,Recall_rate,F1_score]=evaluatePR(X,CorrectIndex,ind)
 N=size(X,1);
inlier_num = length(CorrectIndex);
% outlier_num= N-inlier_num;
inlierRate = inlier_num./N;

Correct_num=sum(double(ismember(ind,CorrectIndex)));
% Error_num=length(ind)-Correct_num;
% unRecall_num=length(CorrectIndex)-Correct_num;
precision_rate=Correct_num/length(ind);
if length(ind)==0
    precision_rate=inlierRate;
end
Recall_rate=Correct_num/length(CorrectIndex);
if precision_rate == 0 && Recall_rate == 0
    F1_score = 0;
else
    F1_score = (2 * precision_rate * Recall_rate)/(precision_rate + Recall_rate);   
end