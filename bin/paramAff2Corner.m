function L = paramAff2Corner( A, sz )

    xl = 1 - sz(1)/2;
    xr = sz(1)/2;
    yd = 1 - sz(2)/2;
    yu = sz(2)/2;
    
    corners = ([A(3) A(4) A(1); A(5) A(6) A(2)] * [[xl;yu;1], [xr;yu;1], [xr;yd;1], [xl;yd;1]]);
    L = corners(:)';

end
