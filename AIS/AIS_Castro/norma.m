% Function normalizes matrix over [0,1]
function [Dn] = norma(D);
% Dn  -> normalized vector over [0,1]
[np,ni] = size(D);
if ni == 1,
   Dn = (D - min(D))./(max(D)-min(D));
else,
   vmaxD = max(D); vminD = min(D);
   for i=1:ni,
      Dn(:,i) = (D(:,i) - vminD(i))./(vmaxD(i)-vminD(i));
   end;
end;
% End Function NORMA
