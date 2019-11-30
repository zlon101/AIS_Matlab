function MainCSA(coverRoot)
fprintf('# start\n');
coverRoot = 'E:\astego\Images\BOSS_ALL\';
payload = single(0.4);

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']); % coverDirs(1)=[];coverDirs(1)=[];
num = length(coverDirs);
bestAbs = cell(num,2);
% save('coverDirs','coverDirs');

old='';t0=datetime('now');
for i = 1:3
    % load coverDirs;
    cPath = [coverRoot, coverDirs(i).name];
    bestAbs{i,1} = coverDirs(i).name;
    % save('bestAbs','bestAbs'); clear bestAbs;
    % clear coverDirs bestAbs bestFits TAbs;
    
    [bestFits,TAbs,~,~] = CSA(cPath,payload);
	% clear CSA; clear mex;
    TAbs = cell2mat(TAbs);
    [vmin,~] = min(bestFits); 
    inds = (bestFits==vmin);
    Ab = min(TAbs(inds));
    
    % load('bestAbs.mat');
    bestAbs{i,2} = Ab;
    % 打印
    msg=sprintf('- count: %3d/%d',i,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
    %fprintf('\n耗时: '); disp(datetime('now')-t0);
    if(mod(i,10)==0)
        save([pwd,'/bestAbs.mat'], 'bestAbs');
    end
end
fprintf('\n耗时: '); disp(datetime('now')-t0);
save([pwd,'\bestAbs.mat'], 'bestAbs');
end