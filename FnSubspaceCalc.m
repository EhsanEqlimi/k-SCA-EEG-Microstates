function [Sub,OrthSub]=FnSubspaceCalc(X,Lables,k,c)
for i=1:max(Lables)

   [row,col,v]= find(Lables==i);
    [Sub(:,:,i),OrthSub(:,:,i)]=FnSubspaceCalcofInleiers(X,col,k);
end
    %