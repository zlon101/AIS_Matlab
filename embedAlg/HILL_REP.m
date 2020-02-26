function stego=HILL_REP(cover, payload)
% 结合HILL和repairPixel
cover= single(imread(cover));

[optP1,optM1] = repairPixel(cover); % 选中的为1，其他为0
[rhoP1,rhoM1] = CostHILL(cover);
vmin= min(min(rhoP1(:)), min(rhoM1(:)));

% 修改的代价
rhoP1(optP1)= vmin*0.1;
rhoM1(optM1)= vmin*0.1;

cover(optP1)= cover(optP1)+1;
cover(optM1)= cover(optM1)-1;

stego = EmbeddingSimulator(cover, rhoP1, rhoM1, payload*numel(cover), false);
end