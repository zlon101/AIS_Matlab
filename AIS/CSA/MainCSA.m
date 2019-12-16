function MainCSA(coverRoot, startInd, endInd,saveRoot)
% startInd='1'; endInd='2';
coverRoot = 'E:\astego\Images\BOSS_ALL\';
payload = single(0.4);
if(~exist('saveRoot','var'))
  saveRoot = 'E:\astego\CSA\';
end
% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']);
num = length(coverDirs);
if(exist([saveRoot,'bestAbs.mat'],'file'))
  load([saveRoot,'bestAbs.mat']);
else
  bestAbs = cell(num,2);
end
if(exist('startInd','var') && str2double(startInd)>0)
  startInd = single(str2double(startInd));
else
  startInd = single(getStart(bestAbs));
end
if(exist('endInd','var'))
  endInd=single(str2double(endInd));
else 
  endInd=single(num);
end
clear getStart;
fprintf('# start\n#count:%d - %d\n',startInd,endInd);

%% 
old=''; t0=datetime('now');
for i = startInd:endInd
  if(~isempty(bestAbs{i,1}))
    continue;
  end
  cPath = [coverRoot, coverDirs(i).name];
  save([saveRoot,'coverDirs.mat'],'coverDirs'); clear coverDirs;
  save([saveRoot,'bestAbs.mat'],'bestAbs'); clear bestAbs;
  clear bestFits TAbs
  t0=tic;
  [bestFits,TAbs] = CSA(cPath,payload);
  disp(toc(t0));
  load([saveRoot,'coverDirs.mat']); load([saveRoot,'bestAbs.mat']);
  [vmin,~] = min(bestFits); 
  inds = (bestFits==vmin);
  TAbs = TAbs(inds,:); [~,ind] = min(TAbs(:,1));
  Ab = TAbs(ind,:);
  bestAbs{i,1} = coverDirs(i).name;
  bestAbs{i,2} = Ab;
   
  if(mod(i,10)==0)
    %save([saveRoot,'bestAbs.mat'], 'bestAbs');
    clear functions mex global;
  end
  % 打印
  msg=sprintf('- count: %3d/%d',i,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
fprintf('\n耗时: '); disp(datetime('now')-t0);
save([saveRoot,'bestAbs.mat'], 'bestAbs');
fprintf('\n# end\n');
end

function start=getStart(Abs)
for i=1:size(Abs,1)
  if(isempty(Abs{i,1}))
    break;
  end
end
start = i;
end