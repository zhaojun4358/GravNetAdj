function tide_v = load_tide_revise(geodeticfai,ymdh,fanweijiao,s)

y = ymdh(1);
m = ymdh(2);
d = ymdh(3);
t = ymdh(4);
Mi = pi/180;
T0 = Julian(y,m,d);
T = (T0  - 2415020.0 + (t-8)/24)/36525;
geodeticlamda = geodeticlamda*Mi;
geodeticfai =  geodeticfai *Mi;
%***************************************************

s = 270.43659 + 481267.89057*T + 0.00198*T^2 + 0.000002*T^3;
h = 279.69668 + 36000.76892*T + 0.0003*T^2;
p = 334.32956 + 4069.03403*T - 0.01032*T^2 - 0.00001*T^3;
N = 259.18328 - 1934.14201*T + 0.00208*T^2 + 0.000002*T^3;
ps = 281.22083 + 1.71902*T + 0.00045*T^2 + 0.000003*T^3;
eslon = 23.45229 - 0.01301*T - 0.000002*T^2;   
s = s*Mi;
h = h*Mi;
p = p*Mi;
N = N*Mi;
ps = ps*Mi;
eslon = eslon*Mi;
%****************************************************

lamda = s - 0.0032*sin(h-ps) - 0.001*sin(2*h-2*p) + 0.001*sin(s-3*h+p+ps )+0.0222*sin(s-2*h+p)+0.0007*sin(s-h-p+ps )-0.0006*sin(s-h) + 0.1098*sin(s-p) - 0.0005*sin(s+h-p-ps)+0.0008*sin(2*s-3*h+ps) +0.0115*sin(2*s-2*h) + 0.0037*sin(2*s-2*p)-0.0020*sin(2*s-2*N) +0.0009*sin(3*s-2*h-p);

lamda_s = h + 0.0335*sin(h-ps)+0.0004*sin(2*(h-ps));

beta = -0.0048*sin(p-N )-0.0008*sin(2*h-p-N) +0.003*sin(s-2*h+N) + 0.0895*sin(s-N) + 0.001*sin(2*s-2*h+p-N) + 0.0049*sin(2*s-p-N)+0.0006*sin(3*s-2*h-N);

hexingshi =(t-8)*15*Mi + h + geodeticlamda- pi;
%**************************************************

Cm_rm = 1 + 0.0545*cos(s-p ) + 0.0030*cos(2*(s-p))+0.01*cos(s-2*h+p) + 0.0082*cos(2*(s-h) )+0.0006*cos(2*s-3*h+ps)+0.0009*cos(3*s -2*h -p);

Cs_rs = 1+0.0168*cos(h-ps) + 0.0003*cos(2*h-2*ps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cm = 384000;             
Cs = 1.496*10^8;         
Mm = 7.349*10^22;        
Ms = 1.9891*10^30;       
reath_radius  = 6371.2 ; 
G = 6.67*10^(-11);       

Dm = 3/4*G*Mm*(earth_radius*10^3)^2/(Cm*10^3)^3;
Ds = 3/4*G*Ms*(earth_radius*10^3)^2/(Cs*10^3)^3;
%*****************************************************
sin_delt_m = sin(eslon) *sin(lamda)*cos(beta) + cos(eslon)*sin(beta);       
cos_delt_m_coa_t_m = cos(beta)*cos(lamda)*cos(hexingshi) + sin(hexingshi)*(cos(eslon)*cos(beta)*sin(lamda) - sin(eslon)*sin(beta)); %cos(delt_m)*cos(t_m)
sin_delt_s = sin(eslon)*sin(lamda_s);
cso_delt_s_cost_t_s = cos(lamda_s)*cos(hexingshi) + sin(hexingshi)*sin(lamda_s);
%*****************************************************

sphericalfai = geodeticfai/Mi - 0.193296*sin(2*geodeticfai);

cos_Zm = sin(sphericalfai)*sin_delt_m + cos(sphericalfai)*cos_delt_m_coa_t_m;
cos_Zs = sin(sphericalfai)*sin_delt_s + cos(sphericalfai)*cso_delt_s_cost_t_s;
sin_Zm = sqrt(1-cos_Zm^2);
sin_Zs = sqrt(1-cos_Zs^2);
cos_Am = (sin_delt_m *cos(sphericalfai) - sin(sphericalfai)*cos_delt_m_coa_t_m)/sin_Zm;
cos_As = (sin_delt_s *cos(sphericalfai) - sin(sphericalfai)*cos_delt_m_coa_t_m)/sin_Zs;
%*****************************************************
g = 9.80665; 
gR = g*earth_radius;
theta_m = 2*Dm/gR*(Cm_rm)^3*2*sin_Zm*cos_Zm + 2*Dm/(g*Cm)*(Cm_rm)^3*(5*cos_Zm^2 - 1)*sin_Zm;
theta_s = 2*Ds/gR*(Cs_rs)^3*2*sin_Zs*cos_Zs;
cos_Am_A = cos_Am*cos(A*Mi) - sqrt(1-cos_Am^2)*sin(A*Mi);
cos_As_A = cos_As*cos(A*Mi) - sqrt(1-cos_As^2)*sin(A*Mi);
%********************************************************
gama = 0.68;
tide_v = (theta_m *cos_Am_A + theta_s *cos_As_A)*gama*s;







