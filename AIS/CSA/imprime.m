function imprime(PRINT,vx,vy,vz,x,y,fx,it,mit)
% x,fx				-> actual values
% vxplot, vplot	-> original (base) function
if PRINT == 1
   if rem(it,mit) == 0
      mesh(vx,vy,vz); hold on; axis([-1 1 -1 1 -1 2.5]);
      xlabel('x'); ylabel('y'); zlabel('f(x,y)');
      plot3(x,y,fx,'k*'); drawnow; hold off;
   end
end
end