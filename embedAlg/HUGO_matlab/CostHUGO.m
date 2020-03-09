function [rhoP1,rhoM1]= CostHUGO(cover)
% HUGO 代价函数
% 返回+1 -1 的代价
%% 
params.gamma= single(1);
params.sigma= single(1);
cover= single(cover);
wetCost= 10^8;
responseP1 = [0; 0; -1; +1; 0; 0];

% create mirror padded cover image
padSize = double(3);
coverPadded = padarray(cover, [padSize,padSize], 'symmetric');
% create residuals
C_Rez_H = coverPadded(:, 1:end-1) - coverPadded(:, 2:end);
C_Rez_V = coverPadded(1:end-1, :) - coverPadded(2:end, :);
C_Rez_Diag = coverPadded(1:end-1, 1:end-1) - coverPadded(2:end, 2:end);
C_Rez_MDiag = coverPadded(1:end-1, 2:end) - coverPadded(2:end, 1:end-1);

rhoM1 = zeros(size(cover),'single');  % declare cost of -1 change           
rhoP1 = zeros(size(cover),'single');  % declare cost of +1 change
for row=1:size(cover, 1)
  for col=1:size(cover, 2)
    D_P1 = 0; D_M1 = 0;
    % Horizontal
    c_rez_sub = C_Rez_H(row+3, col:col+5)';
    stego_sub_P1 = c_rez_sub + responseP1;
    stego_sub_M1 = c_rez_sub - responseP1;
    D_M1 = D_M1 + GetLocalDistortion(c_rez_sub, stego_sub_M1, params);
    D_P1 = D_P1 + GetLocalDistortion(c_rez_sub, stego_sub_P1, params);
    % Vertical
    c_rez_sub = C_Rez_V(row:row+5, col+3);
    stego_sub_P1 = c_rez_sub + responseP1;
    stego_sub_M1 = c_rez_sub - responseP1;
    D_M1 = D_M1 + GetLocalDistortion(c_rez_sub, stego_sub_M1, params);
    D_P1 = D_P1 + GetLocalDistortion(c_rez_sub, stego_sub_P1, params);
    % Diagonal
    c_rez_sub = [C_Rez_Diag(row, col); C_Rez_Diag(row+1, col+1); C_Rez_Diag(row+2, col+2); C_Rez_Diag(row+3, col+3); C_Rez_Diag(row+4, col+4); C_Rez_Diag(row+5, col+5)];
    stego_sub_P1 = c_rez_sub + responseP1;
    stego_sub_M1 = c_rez_sub - responseP1;
    D_M1 = D_M1 + GetLocalDistortion(c_rez_sub, stego_sub_M1, params);
    D_P1 = D_P1 + GetLocalDistortion(c_rez_sub, stego_sub_P1, params);
    % Minor Diagonal
    c_rez_sub = [C_Rez_MDiag(row, col+5); C_Rez_MDiag(row+1, col+4); C_Rez_MDiag(row+2, col+3); C_Rez_MDiag(row+3, col+2); C_Rez_MDiag(row+4, col+1); C_Rez_MDiag(row+5, col)];
    stego_sub_P1 = c_rez_sub + responseP1;
    stego_sub_M1 = c_rez_sub - responseP1;
    D_M1 = D_M1 + GetLocalDistortion(c_rez_sub, stego_sub_M1, params);
    D_P1 = D_P1 + GetLocalDistortion(c_rez_sub, stego_sub_P1, params);
    
    rhoM1(row, col) = D_M1;
    rhoP1(row, col) = D_P1;
  end
end
% clear S_Rez_H S_Rez_V S_Rez_Diag S_Rez_MDiag C_Rez_H C_Rez_V C_Rez_Diag;
rhoM1(rhoM1>wetCost) = wetCost;
rhoP1(rhoP1>wetCost) = wetCost;
rhoP1(cover == 255) = wetCost;
rhoM1(cover == 0) = wetCost;

%% 嵌套函数
  function D = GetLocalDistortion(C_resVect, S_resVect, params)
    % C_resVect: cover image 的残差矩阵
    D = single(0);    
    for i=1:4
     D= D+ GetLocalPotential(C_resVect(i:i+2), S_resVect(i:i+2), params);
    end
    % D= D+ GetLocalPotential(C_resVect, S_resVect, params);
  end

  function Vc = GetLocalPotential(c_res, s_res, params)
    % sqrt 改为 abs
    c_w = (params.sigma + norm(c_res)) .^ (-params.gamma);
    s_w = (params.sigma + norm(s_res)) .^ (-params.gamma);
    Vc= c_w + s_w;
  end
end
