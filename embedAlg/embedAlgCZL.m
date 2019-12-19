function [stego,rhoP1,rhoM1]=embedAlgCZL(I,payload)
% 失真函数设计
%%
HFilter = [-1, 2, -2, 2,-1;
           2,-6,  8,-6, 2;
          -2, 8,-12, 8,-2;
           2,-6,  8,-6, 2;
          -1, 2, -2, 2,-1];
L1 = ones(3);
L2 = ones(15);
Y1 = imfilter(I,HFilter,'symmetric','same');
Y2 = imfilter(abs(Y1),L1,'symmetric','same');
Cost = imfilter(Y2.^-1,L2,'symmetric','same');
Cost = Cost./sum(L1(:))./sum(L2(:));
% adjust embedding costs
wetCost = 10^8;
Cost(Cost > wetCost) = wetCost;
Cost(isnan(Cost)) = wetCost;
rhoP1 = Cost;  % +1 的代价
rhoM1 = Cost;
rhoP1(I==255) = wetCost;
rhoM1(I==0) = wetCost;

%%
stego = EmbeddingSimulator(I, rhoP1, rhoM1, payload*numel(I), false);

  function [y] = EmbeddingSimulator(x, rhoP1, rhoM1, m, fixEmbeddingChanges)
    n = numel(x);   
    lambda = calc_lambda(rhoP1, rhoM1, m, n);
    pChangeP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
    pChangeM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
    if fixEmbeddingChanges == 1
      RandStream.setGlobalStream(RandStream('mt19937ar','seed',139187));
    else
      RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
    end
    randChange = rand(size(x));
    y = x;
    y(randChange < pChangeP1) = y(randChange < pChangeP1) + 1;
    y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;

    function lambda = calc_lambda(rhoP1, rhoM1, message_length, n)
      L3 = 1e+3;
      m3 = double(message_length + 1);
      iterations = 0;
      while m3 > message_length
        L3 = L3 * 2;
        pP1 = (exp(-L3 .* rhoP1))./(1 + exp(-L3 .* rhoP1) + exp(-L3 .* rhoM1));
        pM1 = (exp(-L3 .* rhoM1))./(1 + exp(-L3 .* rhoP1) + exp(-L3 .* rhoM1));
        m3 = ternary_entropyf(pP1, pM1);
        iterations = iterations + 1;
        if (iterations > 10)
          lambda = L3;
          return;
        end
      end        

      l1 = 0; 
      m1 = double(n);        
      lambda = 0;
      alpha = double(message_length)/n;
      % limit search to 30 iterations
      % and require that relative payload embedded is roughly within 1/1000 of the required relative payload        
      while  (double(m1-m3)/n > alpha/1000.0 ) && (iterations<30)
        lambda = l1+(L3-l1)/2; 
        pP1 = (exp(-lambda .* rhoP1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
        pM1 = (exp(-lambda .* rhoM1))./(1 + exp(-lambda .* rhoP1) + exp(-lambda .* rhoM1));
        m2 = ternary_entropyf(pP1, pM1);
        if m2 < message_length
          L3 = lambda;
          m3 = m2;
        else
          l1 = lambda;
          m1 = m2;
        end
        iterations = iterations + 1;
      end
    end

    function Ht = ternary_entropyf(pP1, pM1)
      p0 = 1-pP1-pM1;
      P = [p0(:); pP1(:); pM1(:)];
      H = -((P).*log2(P));
      H((P<eps) | (P > 1-eps)) = 0;
      Ht = sum(H);
    end
  end
end