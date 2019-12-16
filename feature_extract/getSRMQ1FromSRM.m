function SRMQ1 = getSRMQ1FromSRM(SRM)
% 从 SRM 特征中提取 SRMQ1 子特征
%%
indLower= [1,6501,11376,21126,25026,28926,29602,30616,31630,33658]';
indUpper= [3250,8125,14625,22425,26325,29263,29939,30953,32305,33995]';
SRMQ1= zeros(size(SRM,1),12753,'single');
start= 1;
for i=1:length(indLower)
  L = indUpper(i)-indLower(i);
  SRMQ1(:,start:start+L) = SRM(:, indLower(i):indUpper(i));
  start = start+L+1;
end
end