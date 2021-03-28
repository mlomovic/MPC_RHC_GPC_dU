clear,clc;
%Uzimamo vec diskretizovan sistem sa periodom 1s
Ap=[1.1 2;0 0.95]
Bp=[0;0.0787]
Cp=[-1 1]
Dp=0
[Ac,Bc,Cc,Dc] = d2cm(Ap,Bp,Cp,Dp,0.1)

Np=4;
Nc=4;
r_omega=0.1;

N_sim=50; % Broj trenutaka odabiranja


sis = ss(Ap,Bp,Cp,Dp)

[Phi,F,Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e] = mpcPojacanja(Ap,Bp,Cp,Nc,Np);
[n,n_in]=size(B_e);
xm=[0;0];
Xf=zeros(n,1);

r=ones(N_sim,1);
u=0; % u(k-1) =0
y=0;

for kk=1:N_sim;
DeltaU = inv(Phi_Phi+r_omega*eye(Nc,Nc)) * (Phi_R*r(kk)-Phi_F*Xf);
deltau=DeltaU(1,1);
u=u+deltau;
u1(kk)=u;
y1(kk)=y;
xm_old=xm;
xm=Ap*xm+Bp*u;
y=Cp*xm;
Xf=[xm-xm_old;y];
end

k=0:(N_sim-1);

figure
s(1)=subplot(211)
plot(k,y1,'-b',k,r,'--b'),grid
title(s(1),['Upravljanje MPC-om sa parametrima: N_c = ',num2str(Nc),', N_p = ',num2str(Np),', r_\omega = ',num2str(r_omega)]);
xlabel('Trenutci odabiranja')
legend('Izlaz','Zeljena vrednost')
s(2)=subplot(212)
plot(k,u1),grid
xlabel('Trenutci odabiranja')
legend('Upravljanje')
