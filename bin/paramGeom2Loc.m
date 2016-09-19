function L = paramGeom2Loc( G, objWidth )
    L = [G(1), G(2), G(3)*objWidth,  G(5)*G(3)*objWidth, G(4)];
end

