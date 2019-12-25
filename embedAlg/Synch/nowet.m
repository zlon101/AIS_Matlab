function y = nowet(y)

l = sort( unique(y), 'descend' );
M = l( find(l < 10^7, 1, 'first') );

y( y > 10^7 ) = M;

end