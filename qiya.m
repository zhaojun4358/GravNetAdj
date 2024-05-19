
function deltga = qiya(P,H)

Pn = 1.01325*10^3*(1-(0.0065*H/288.15)^5.2559);
deltga = 0.3*(P-Pn);