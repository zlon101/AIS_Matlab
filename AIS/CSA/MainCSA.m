function MainCSA(coverRoot, startInd, endInd)
% startInd='1';endInd='10';
% coverRoot = 'E:\astego\Images\BOSS_ALL\';
payload = single(0.4);
saveRoot = 'E:\astego\CSA\';

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']);
num = length(coverDirs);
if(exist('endInd','var'))
  endInd=int8(str2double(endInd));
else 
  endInd=num;
end
if(exist([saveRoot,'bestAbs.mat'],'file'))
  load([saveRoot,'bestAbs.mat']);
  startInd = getStart(bestAbs);
else
  bestAbs = cell(num,2);
  if(exist('startInd','var'))
    startInd=int8(str2double(startInd));
  else
    startInd=1;
  end
end
clear getStart;
fprintf('# start\n#count:%d - %d\n',startInd,endInd);

%% 
old=''; % t0=datetime('now');
for i = startInd:endInd
  cPath = [coverRoot, coverDirs(i).name];
  save([saveRoot,'coverDirs.mat'],'coverDirs'); clear coverDirs;
  save([saveRoot,'bestAbs.mat'],'bestAbs'); clear bestAbs;
  
  [bestFits,TAbs] = CSA(cPath,payload);
  load([saveRoot,'coverDirs.mat']); 
  load([saveRoot,'bestAbs.mat']);
  [vmin,~] = min(bestFits); 
  inds = (bestFits==vmin);
  Ab = min(TAbs(inds));
  bestAbs{i,1} = coverDirs(i).name;
  bestAbs{i,2} = Ab;
   
  if(mod(i,10)==0)
    save([saveRoot,'bestAbs.mat'], 'bestAbs');
    clear functions mex global;
  end
  % 打印
  msg=sprintf('- count: %3d/%d',i,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
fprintf('\n耗时: '); disp(datetime('now'));
save([saveRoot,'bestAbs.mat'], 'bestAbs');
fprintf('\n# end');
end

function start=getStart(Abs)
for i=1:size(Abs,1)
  if(isempty(Abs{i,1}))
    break;
  end
end
start = i;
end