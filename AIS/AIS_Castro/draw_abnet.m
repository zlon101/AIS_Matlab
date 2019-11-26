function draw_abnet(w);

[N,L] = size(w);
r = [0.2.*ones(1,N);1:1:N]';
c = [1.2.*ones(1,L);1:1:L]';
if N > L,
   y = N;
else, 
   y = L;
end;
plot(r(:,1),r(:,2),'sk'); hold on;
plot(c(:,1),c(:,2),'or');
axis([0 1.4 0 y+1]);
for i = 1:N,
   for j = 1:L,
      if w(i,j) ~= 0,
         line([r(i,1),c(j,1)],[r(i,2),c(j,2)]);
      end;
   end;
end;
title('ABNET (0-weights ommitted)');
hold off;