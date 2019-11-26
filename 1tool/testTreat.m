function tredObj = testTreat(cover, changeRates, tredPath, correlaRata, Dests)
% 免疫测试, 修改方案
% cover可以是路径, jobj结构体
% tredPath, Dests是可选项
% 输出tredObj
% changeRates:        被修改的系数的比例
% correlaRata         相关性比例
% -------------------------------------------------------------------------
seed=99;rng('default');rng(seed);
if(ischar(cover))
    [jobj, quality, QTable, Cb, Cr] = getQDCT(cover);
    QDCT=jobj.coef_arrays{1};
else
    QDCT = cover.coef_arrays{1};
    jobj = cover;
end
% 4 1;3 2;2 3;1 4
mode=[
    2 1;1 2;    
   3 1;2 2;1 3;
  4 1;3 2;2 3;1 4      % Q=90
  ];
o=QDCT;

marks=zeros(size(QDCT));
for i=1:size(mode,1)   
    marks(mode(i,1):8:end, mode(i,2):8:end) = true;
end
% 低频系数线性索引
LFInd = find(marks);

%--------------------------------------------------------------------------
nChange = length(changeRates);
Srcs = 0:nChange-1;
if(~exist('Dsts', 'var') || length(Dests)~=length(Srcs))
    Dests = Srcs+1;
end
for i=1:nChange
    % 指定系数下标
    embedMarks = ( abs(QDCT)==Srcs(i) );
    embedMarks(1:8:end,1:8:end) = false;
    specificInd = find(embedMarks);
    linearInd = intersect(LFInd, specificInd);
    specificInd = linearInd(randperm(length(linearInd)));   % 指定系数位置
    if(changeRates(i)>0)
        QDCT = moverValue(QDCT, specificInd, Srcs(i), Dests(i), changeRates(i));
    end
end

% 增强块间相关性
if(exist('correlaRata', 'var') && correlaRata>0 && correlaRata<1)
    QDCT=VCO_QDCT(QDCT, correlaRata);
end


% QDCT=HCO_QDCT(QDCT, 0.025);

%{
lena图像测试
% 将20%的零系数移至[1,2,3]
QDCT = moverValue(QDCT, lfInd, 0, [1,2,3],0.2);
QDCT=VCO_QDCT(QDCT, 0.5);

man图像测试
QDCT = moverValue(QDCT, lfInd, 0, [4,6,6],0.5);
QDCT=VCO_QDCT(QDCT, 0.5);
%}
% 保存免疫载体图像
tredObj = jobj;
tredObj.coef_arrays{1} = QDCT;
if(exist('tredPath', 'var') && ischar(tredPath))
    jpeg_write(tredObj, tredPath);
end
end

%% 指定值的系数集-0
function QDCT=moverValue(QDCT, specificInd, Rst, Dsts, changeRate)
% 移动一类值为Rst的系数到Dsts，如将20%的零值系数移至[1,2,3]

T = zeros(2*length(Dsts), 1);
T(1:2:end)=Dsts; T(2:2:end)=-1*Dsts;
Dsts=T;

linearInd = specificInd( 1:round(length(specificInd)* changeRate) );
% 将linearInd分nSections段,一行一段，每段移动到dst的一个值
nSections = length(Dsts);
% 剔除多余的部分
rest = mod(length(linearInd), nSections);
linearInd(1:rest)=[];
linearInd = reshape(linearInd, nSections,[]);

% 移动
if(Rst==0)
    for i=1:nSections
        QDCT( linearInd(i,:) ) = Dsts(i);
    end
else
    for i=1:nSections
        character = sign(QDCT( linearInd(i,:)));        
        QDCT( linearInd(i,:)) = character .* abs(Dsts(i));
    end
end
end

%% 25为共生矩阵,增强垂直方向块间相关性
function QDCT=VCO_QDCT(QDCT, changeRate)
% 
seed=99;rng('default');rng(seed);
QDCT_O = QDCT;
DC = QDCT(9:end,:); DC = DC(1:8:end, 1:8:end);
QDCT(1:8:end,1:8:end) = 1e8;
% vertical part
mBlockV = reshape(QDCT(1:end-8,:), [],1);
nBlockV = reshape(QDCT(9:end,:), [],1);
% Mv:7*7， Mv(i,j)表示第k块值为i，第k+1块值为j的频数，i,j范围为[-3,3]

% 第k块值=1，第k+1块值=0
ind = [];
BUF = find( (abs(mBlockV)==1) & (nBlockV==0) );
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];

BUF = find((mBlockV==2) & (nBlockV==1));
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];

BUF = find((mBlockV==-2) & (nBlockV==-1));
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];
% 修改
nBlockV(ind) = mBlockV(ind);
nBlockV=reshape(nBlockV, size(QDCT(9:end,:)) );
nBlockV(1:8:end, 1:8:end) = DC;
QDCT = [QDCT_O(1:8,:);nBlockV];
end

%% 25为共生矩阵,增强水平方向块间相关性
function QDCT=HCO_QDCT(QDCT, changeRate)
% 
QDCT_O = QDCT;  DC = QDCT(:,9:end);
DC = DC(1:8:end, 1:8:end);
QDCT(1:8:end,1:8:end) = 1e8;
% 水平方向
mBlockV = reshape(QDCT(:, 1:end-8), [],1);
nBlockV = reshape(QDCT(:, 9:end), [],1);
% Mv:7*7， Mv(i,j)表示第k块值为i，第k+1块值为j的频数，i,j范围为[-3,3]

% 第k块值=1，第k+1块值=0
ind = [];
BUF = find( (abs(mBlockV)==1) & (nBlockV==0) );
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];

BUF = find((mBlockV==2) & (nBlockV==1));
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];

BUF = find((mBlockV==-2) & (nBlockV==-1));
BUF = BUF(randperm(length(BUF)));
BUF = BUF( 1: round(length(BUF)*changeRate) );
ind = [ind;BUF];
% 修改
nBlockV(ind) = mBlockV(ind);
nBlockV=reshape(nBlockV, size(QDCT(:,9:end)) );
nBlockV(1:8:end, 1:8:end) = DC;
QDCT = [QDCT_O(:,1:8), nBlockV];
end