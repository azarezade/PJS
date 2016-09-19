function  corners = paramLoc2Corner( L )
    center = [L(1), L(2)];
    width = L(3);
    height = L(4);
    corners = [L(1)-(width/2), L(1)+(width/2), L(1)+(width/2), L(1)-(width/2);
    L(2)-(height/2), L(2)-(height/2), L(2)+(height/2), L(2)+(height/2)];
    t=L(5); A=[cos(t) -sin(t); sin(t) cos(t)];
     corners = corners - repmat([center(1);center(2)],1,4);
     corners = A * corners;
     corners = reshape(corners, 1, 8);
     corners = corners + repmat(center,1,4);
end
