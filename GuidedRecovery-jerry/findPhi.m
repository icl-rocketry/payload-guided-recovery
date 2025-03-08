function phi = findPhi(l,d,dl)   
    theta = acos(0.5*d/l); 
    c = sqrt(d^2 + dl^2 - 2*d*dl*cos(theta));
    phi = acos( (l^2 - c^2 - (l-dl)^2) / (-2*c*(l-dl)) ) - theta;


end