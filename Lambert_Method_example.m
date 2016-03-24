%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg






% Lambert Method example
close all; clear all; clc;
%% inputs
r1=[5e3, 1e4, 2100];   %  first position vector (1x3)
r2=[-14600, 2500, 7e3];   %  second position vector (1x3)
dt=3600;   %  the time between r1 and r2 
z0=1.5;   %  initial value of z
state=0; % 0 for prograde trajectory and 1 for retrograde trajectory
muo=398600; %  Gravitational Parameter
Re=6378;    % raduis of plant
%% Lambert Solution
[ h, mag_h, i, omega, e_vector, mag_e, w, theta, rp, zp, epslon ] = Lambert ( r1, r2, dt, z0, state, muo, Re );
%% RK4 parameter of orbital element 
order=6;
X0=[mag_h;mag_e;i;omega;w;theta];
B=[0;0;0;0;0;0;];
sol(1:6,1)=X0;
dt=1000;
t_initial=0;
t_final=2e7;
%% solution of orbital element 
for n=1:length(t_initial:dt:t_final)
    A=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      muo^2/sol(1,n)^4*(1+sol(2,n)*cosd(sol(6,n)))^2,0,0,0,0,0];    
    [ XX, t ] = RK4( A,B,sol(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol(1:6,n+1)=XX(1:6,2);
end
for m=1:length(sol(1,:))
    [ r, v ] = OrbitalElements2rvGeo( sol(1,m), muo, sol(2,m), sol(3,m), sol(4,m), sol(5,m), sol(6,m) );
    r_XYZ(1:3,m)=r;
    v_XYZ(1:3,m)=v;
end
%% RK4 parameter of orbital element with oblateness
order1=6;
X01=[mag_h;mag_e;i;omega;w;theta];
B1=[0;0;0;0;0;0;];
sol1(1:6,1)=X0;
J2=1.08263e-3;
%% solution of orbital element with oblateness
for n=1:length(t_initial:dt:t_final)
    A1=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*cosd(sol1(3,n)),0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*(5/2*(sind(sol1(3,n)))^2-2),0,0,0,0,0;...
      muo^2/sol1(1,n)^4*(1+sol1(2,n)*cosd(sol1(6,n)))^2,0,0,0,0,0];    
    [ XX1, t1 ] = RK4( A1,B1,sol1(1:6,n),dt,n*dt,(n+1)*dt,order1 );
    sol1(1:6,n+1)=XX1(1:6,2);
end
for m=1:length(sol1(1,:))
    [ r1, v1 ] = OrbitalElements2rvGeo( sol1(1,m), muo, sol1(2,m), sol1(3,m), sol1(4,m), sol1(5,m), sol1(6,m) );
    r_XYZ1(1:3,m)=r1;
    v_XYZ1(1:3,m)=v1;
end
%% plotting
figure(1);
view(3);
set(gcf,'color','w');
subplot(1,2,1)
plot3(r_XYZ(1,:),r_XYZ(2,:),r_XYZ(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution without oblateness','fontsize',18);
subplot(1,2,2)
plot3(r_XYZ1(1,:),r_XYZ1(2,:),r_XYZ1(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution with oblateness','fontsize',18);
%% error 
for k=1:length(r_XYZ1(2,:))
    E1(1,k)=norm([r_XYZ(1:3,k)]);
    E2(1,k)=norm([r_XYZ1(1:3,k)]);
end
[ E,Max_e,std_e, mean_e, RMS_e ] = ERROR ( E1(1,:),E2 );
figure(2);
set(gcf,'color','w');
plot(t_initial:dt:t_final+dt,E)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _o_b_l-Sol_w_i_t_h _o_b_l','fontsize',18);
title('Error','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e) ' Std(Error) = ' num2str(std_e) ' mean(Error) = ' num2str(mean_e) ' RMS(Error) = ' num2str(RMS_e) ]);
grid on;