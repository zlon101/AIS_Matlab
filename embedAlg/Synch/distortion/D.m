function y = D( s, rho )
y = 0;
indices = [ 1 2; 2 1; 2 3; 3 2 ];
for i = 1:4
    s1 = s( 2, 2 );
    s2 = s( indices(i,1), indices(i,2) );
    r1 = rho( 2, 2 );
    r2 = rho( indices(i,1), indices(i,2) );
    A = ( r1 + r2 ) / 2;
    
    if     s1 == -1
        if     s2 == -1
            y = y + 0;
        elseif s2 ==  0
            y = y + A;
        else
            y = y + 5 * A;
        end
    elseif s1 ==  0
        if     s2 == -1
            y = y + A;
        elseif s2 ==  0
            y = y + 0;
        else
            y = y + A;
        end
    else
        if     s2 == -1
            y = y + 5 * A;
        elseif s2 ==  0
            y = y + A;
        else
            y = y + 0;
        end
    end
end
end