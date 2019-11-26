function EvolDirect = coefValueQuery(diff)
% 查询差异最大的维度及其系数值
% increaseCoef: 该值的数量增加了
% decreaseCoef: 该值的数量减少了

absDiff = abs(diff);
[~, ind] = sort(absDiff, 'descend');
maxInd = ind( 1:round(0.02 * length(diff)) );
ind = maxInd>274;
a = zeros(size(maxInd));
a(ind) = 274;
dims = maxInd-a;        % dim,维度分量
increaseDims = dims( diff(maxInd)>0 );
decreaseDims = dims( diff(maxInd)<0 );

increaseCoef = dimToCoef(increaseDims);
decreaseCoef = dimToCoef(decreaseDims);
symbioticMatrix = [];
if(nnz(dims>168))
    symbioticMatrix = blockCorrelation( dims(dims>168), diff );
end
EvolDirect = struct('increaseCoef', increaseCoef, 'decreaseCoef', decreaseCoef,...
    'symbioticMatrix', symbioticMatrix);
end

function coefs=dimToCoef(dims)
% 维度隐射到系数
% 值为零对应的dims
% 非零值表示绝对值

zero =  [17,28,39,50,61, 112:120];
one =   [16,27,38,49,60,18,29,40,51,62, 121:129, 103:111];
two =   [15,26,37,48,59,19,30,41,52,63, 130:130+8, 94:102];
three = [14,25,36,47,58,20,31,42,53,64, 139:139+8, 85:93];
four =  [13,24,35,46,57,21,32,43,54,65, 76:84, 148:156];
five =  [12,23,34,45,56,22,33,44,55,66, 67:75, 157:165];

% ----------------------------------------------------------
countZero = 0;countOne = 0;countTwo = 0;
countThree = 0;countFour = 0;countFive = 0;
for i=1:length(dims)
    if( ismember( dims(i),zero ) )
        countZero=countZero+1;
    elseif(ismember( dims(i),one ))
        countOne=countOne+1;
    elseif(ismember( dims(i),two ))
        countTwo=countTwo+1;
    elseif(ismember( dims(i),three ))
        countThree=countThree+1;
    elseif(ismember( dims(i),four ))
        countFour=countFour+1;
    elseif(ismember( dims(i),five ))
        countFive=countFive+1;
    end
end
coefs = [0,countZero;
         1,countOne;
         2,countTwo;
         3,countThree;
         4,countFour;
         5,countFive];
[~,ind] = sort( coefs(:,end), 'descend' );
coefs = coefs(ind,:);
end

function symbioticMatrix = blockCorrelation(dims,diff)
% 块间相关性隐射
dims = dims(dims<194);
PN = zeros(size(dims));
PN( diff(dims)>0 ) = 1;    % 大于0或小于0
PN( diff(dims)<0 ) = -1;

dims = dims-169;
% 取商
r = 1+fix(dims/5)-3;
% 取余
c = 1+mod(dims,5)-3;
% 共生矩阵
symbioticMatrix = [r',c',PN'];
% (a,b,出现次数)
end
