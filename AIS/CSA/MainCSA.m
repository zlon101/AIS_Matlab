function MainCSA(coverRoot)
fprintf('# start\n');
coverRoot = 'E:\astego\Images\BOSS_ALL\';
saveRoot = 'E:\astego\CSA\';
payload = single(0.4);

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']); % coverDirs(1)=[];coverDirs(1)=[];
num = length(coverDirs);
bestAbs = cell(num,2);

old='';t0=datetime('now');
for i = 1:num
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
    clear functions;
  end
  % 打印
  msg=sprintf('- count: %3d/%d',i,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
  %fprintf('\n耗时: '); disp(datetime('now')-t0);
end
fprintf('\n耗时: '); disp(datetime('now')-t0);
save('E:\astego\CSA\bestAbs.mat', 'bestAbs');
fprintf('\n# end');
end