function [idxUnique,ID_Bin] = removeRepeat(X,Y)
% CC = UniqueSearch(X);
% EE = UniqueSearch(Y);
% ID_Bin = (CC&EE);
% idxUnique=find(ID_Bin==1);

CC= UniqueSearch_revision(X);
EE = UniqueSearch_revision(Y);
if length(CC)==2
    CC = CC';
end
if length(EE)==2
    EE = EE';
end
double_rep = intersect(CC,EE,'rows');
single_rep = setdiff([CC;EE],double_rep,'rows');
if size(single_rep,1)==1
   single_rep = single_rep';
end
if size(double_rep,1)==1
   double_rep = double_rep';
end
remove_ps = [unique(single_rep);double_rep(:,1)];
idxUnique = setdiff(1:size(X,1),remove_ps);

function indNew = UniqueSearch(X)
ind0 = false(size(X,1),1);
[aa,bb] = sortrows(X,1);
cc0=logical(sum(abs(diff(aa)),2)); % cc0 = (sum(abs(diff(aa)),2)) > 2
cc1 = [true;cc0];
cc2 = [cc0;true];
idx=cc1&cc2;
ind0(bb(idx))= true;
indNew=ind0;
end

function indNew = UniqueSearch_revision(X)
[aa,bb] = sortrows(X,1);
cc0 = (sum(abs(diff(aa)),2)) > 2;
rep_ps = find(cc0==0);
indNew(:,1) = rep_ps;
indNew(:,2) = rep_ps+1;
indNew = bb(indNew);
indNew = sort(indNew,2);
end


end