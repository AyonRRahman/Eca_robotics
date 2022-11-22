function out = PilotDepth_LQR(in, memory, parameters)
% Depth piloting
% coder.extrinsic('lqrd');

Q = diag([1e-6,1,100,10,0.1]);
R = 1e9;

%% Variables intermédiaires
un_demi_rho_S = 0.5*1026*0.385;
m_kg = 1974.26 + 1845.9; Det33 = 1/m_kg;
I_kgm2= 4173.5 + 3388.76; Det55 = 1/I_kgm2;
Lref = 5;
CZuw = -3.1;
CZuq = -2;
CMuw = 1;
CMuq = -0.9;
W = 9.81*1974.26;
zg = 3e-2;
CZ = -0.41;
CM = -0.17;
%% Etats Equilibres
CM0 = -0.002;
CZ0 = -0.01;
theta0 = (CM*CZ0 - CM0*CZ)*in.u_ms*2*un_demi_rho_S*Lref*abs(in.u_ms)...
    /((- 2.*W*zg + CMuw*in.u_ms^2*2*un_demi_rho_S*Lref)*CZ - in.u_ms^2*2*un_demi_rho_S*Lref*CM*CZuw);
BAR0 = 2*((-un_demi_rho_S*CMuw*in.u_ms^2*sin(theta0)*Lref/cos(theta0) ...
    - (- 1.*W*zg)*sin(theta0))/(2*un_demi_rho_S*in.u_ms*Lref*abs(in.u_ms)) - 0.5*CM0)/CM;
w0 = in.u_ms*tan(theta0);
q0 = 0;

%% A matrix
A = zeros(5,5);

A(1,1) = un_demi_rho_S*Det33*CZuw*in.u_ms;
A(1,2) = Det33*((un_demi_rho_S*CZuq*Lref + m_kg)*in.u_ms + 2*m_kg*q0*zg);

A(2,1) = Det55*(un_demi_rho_S*CMuw*in.u_ms*Lref - 1.*m_kg*q0*zg);
A(2,2) = Det55*(un_demi_rho_S*CMuq*in.u_ms*Lref - 1.*m_kg*w0*zg);
A(2,4)= -Det55*W*zg*cos(theta0);

A(3,1) = 1;
A(4,2) = 1;
A(5,3) = 1;

%% B matrix      
B = zeros(5,1);
B(1,1) = un_demi_rho_S*Det33*CZ*in.u_ms*abs(in.u_ms);
B(2,1) = un_demi_rho_S*CM*in.u_ms*Lref*abs(in.u_ms)*Det55;

C = eye(5);

%%
delta_z = in.zc_m - in.z_m;
delta_z = EcaF_Saturate(delta_z,-parameters.delta_z_sat_m,parameters.delta_z_sat_m);

delta_w = w0 - in.w_ms;
delta_q = q0 - in.q_rads;
delta_theta = DiffAngle(theta0,in.theta_rad);

Qy = C'*Q*C;

K = zeros(1,5);
K = lqrd(A,B,Qy,R,in.delta_time_s);

%in.delta_time_s

BARc_rad = BAR0 - delta_w.*K(1) - delta_q.*K(2) - delta_z.*K(3)...
    - delta_theta.*K(4) - memory.int_z.*K(5);

%% Output saturation
out.BARc = EcaF_Saturate(BARc_rad,parameters.BAR_sat_rad(1),parameters.BAR_sat_rad(2));

end
