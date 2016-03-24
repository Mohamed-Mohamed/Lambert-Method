function [ h, mag_h, i, omega, e_vector, mag_e, w, theta1, rp, zp, epslon ] = Lambert ( r1, r2, dt, z0, state, muo, R )
% this function is about Lambert Method of determination the orbit by 2 position vectors and the time between them
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% r1   :  first position vector (1x3)
% r2   :  second position vector (1x3)
% dt   :  the time between r1 and r2 
% zo   :  initial value of z
% state : 0 for prograde trajectory and 1 for retrograde trajectory
% muo :  Gravitational Parameter
%% OUTPUTS:
% h              : specific angular momentum vector
% mag_h     : specific angular momentum magnitude
% i               : inclination angle in degree
% omega     : right ascension of the ascending node in degree
% e              : eccentricity vector
% mag_e     : eccentricity magnitude
% w             : argument of perigee in degree
% theta1      : true anomaly of r1 in degree
% epslon      :  specific energy
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
mag_r1=norm(r1);
mag_r2=norm(r2);
r1CROSSr2=cross(r1,r2);
if state == 0 && r1CROSSr2(3) >= 0
    dtheta=acosd(dot(r1,r2)/mag_r1/mag_r2);
elseif state == 0 && r1CROSSr2(3) < 0
    dtheta=360-acosd(dot(r1,r2)/mag_r1/mag_r2);
elseif state == 1 && r1CROSSr2(3) < 0
    dtheta=acosd(dot(r1,r2)/mag_r1/mag_r2);
elseif state == 1 && r1CROSSr2(3) >= 0
    dtheta=360-acosd(dot(r1,r2)/mag_r1/mag_r2);
end
A=sind(dtheta)*sqrt(mag_r1*mag_r2/(1-cosd(dtheta)));
% syms  n m z
% S(n)=1/6-n/120+n^2/5040-n^3/362880+n^4/39916800-n^5/6227020800+n^6/1.307674368e+12
% C(m)=1/2-m/24+m^2/720-m^3/40320+m^4/3628800-m^5/479001600+m^6/8.71782912e+10

syms n m z;
if z0>0 % ellipse
    S(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
    C(m)=(1-cos(sqrt(m)))/m;
elseif z0<0 % hyperbola
    S(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
    C(m)=(cosh(sqrt(-m))-1)/-m;
elseif z0==0 % parabola
    S(n)=sym(1/6);
    C(m)=sym(1/2);
end

Y(z)=mag_r1+mag_r2+A*(z*S(z)-1)/sqrt(C(z));
F(z)=(Y(z)/C(z))^(3/2)*S(z)+A*sqrt(Y(z))-sqrt(muo)*dt;
if z0 == 0
    F_dash(z)=sqrt(2)/40*(Y(0))^(3/2)+A/8*(sqrt(Y(0))+A*sqrt(1/2/Y(0)));
else
    F_dash(z)=(Y(z)/C(z))^(3/2)*(1/2/z*(C(z)-3/2*S(z)/C(z))+3/4*(S(z))^2/C(z))+A/8*(3*S(z)/C(z)*sqrt(Y(z))+A*sqrt(C(z)/Y(z)));
end
sol=z0;
ratio=double(((F(z0))/(F_dash(z0))));
while double(abs(ratio))>= 1e-8;
if sol>0 % ellipse
    S(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
    C(m)=(1-cos(sqrt(m)))/m;
elseif sol<0 % hyperbola
    S(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
    C(m)=(cosh(sqrt(-m))-1)/-m;
elseif sol==0 % parabola
    S(n)=1/6;
    C(m)=1/2;
end
Y(z)=mag_r1+mag_r2+A*(z*S(z)-1)/sqrt(C(z));
F(z)=(Y(z)/C(z))^(3/2)*S(z)+A*sqrt(Y(z))-sqrt(muo)*dt;
    if sol == 0
        F_dash(z)=double(sqrt(2)/40*(Y(0))^(3/2)+A/8*(sqrt(Y(0))+A*sqrt(1/2/Y(0))));
    else
        F_dash(z)=(Y(z)/C(z))^(3/2)*(1/2/z*(C(z)-3/2*S(z)/C(z))+3/4*(S(z))^2/C(z))+A/8*(3*S(z)/C(z)*sqrt(Y(z))+A*sqrt(C(z)/Y(z)));
    end
    sol=double(sol-ratio);
    ratio=double((F(sol))/(F_dash(sol)));
end
Z=sol;
y=double(Y(Z));
f=1-y/mag_r1;
g=A*sqrt(y/muo);
g_dot=1-y/mag_r2;
v1_vector=1/g*(r2-f*r1);
v2_vector=1/g*(g_dot*r2-r1);
v1_mag=norm(v1_vector);
v2_mag=norm(v2_vector);
% orbital element procdure
vr1=dot(r1,v1_vector)/mag_r1;
h=cross(r1,v1_vector);
mag_h=norm(h);
i=acosd(h(3)/mag_h);
N=cross([0,0,1],h);
mag_N=norm(N);
if N(2) >= 0
    omega=acosd(N(1)/mag_N);
elseif N(2) < 0
    omega=360-acosd(N(1)/mag_N);
end
e_vector=((v1_mag^2-muo/mag_r1)*r1-mag_r1*vr1*v1_vector)/muo;
mag_e=norm(e_vector);
if e_vector(3) >= 0
    w=acosd(dot(N,e_vector)/mag_N/mag_e);
elseif e_vector(3) < 0
    w=360-acosd(dot(N,e_vector)/mag_N/mag_e);
end
if vr1 >= 0
    theta1=acosd(dot(e_vector,r1)/mag_r1/mag_e);
elseif vr1 < 0
    theta1=360-acosd(dot(e_vector,r1)/mag_r1/mag_e);
end
rp=mag_h^2/muo/(1+mag_e*cosd(0));
zp=rp-R;
epslon=-1/2*muo^2/mag_h^2*(1-mag_e^2);
end