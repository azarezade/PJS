function G = paramLoc2Geom( L, objWidth)
    % G(1,2) = center     G(3) = width/32,  G(4) = ratoation, G(5) =
    % height /width  G(6) = skew = 0
    G = [L(1), L(2), L(3)/objWidth, L(5), 0, 0];
    G(5) = L(4)/(G(3)*objWidth);
end

