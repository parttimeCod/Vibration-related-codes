k=12.35;%Stiffness of system
c=0.274;%Damping constant calculated previously
m=0.2;%System mass given in problem
h=0.001;%Time step maximum time step allowed is 0.08sec
Rt=15;%Elastic yielding force for tension
Rc=-15;%Elastic yielding force for compression
yt=Rt/k;
yc=Rc/k;
ymin=1000;%Time is any random varriable out of the range of response used as flag varriable to identify minimum plastic response point
restime=5; %response to be determined for this time length
%These are definations of forcing functions
t1=0:0.001:0.45;
f1=(20/0.45)*t1;
t2=0.001:0.001:0.65;
f2=(20/0.65)*(0.65-t2);
t3=0.001:0.001:0.1;
f3=(-10/0.1)*(t3);
t4=0.001:0.001:0.2;
f4=(-10/0.2)*(0.2-t4);
f5=zeros(1,round((restime-1.4)/h));
F=[f1,f2,f3,f4,f5];
%Varriable decleration
u=zeros(1+restime/h,1); %Displacement
v=zeros(1+restime/h,1); %Velocity
a=zeros(1+restime/h,1); %Acceleration
R=zeros(1+round(restime/h),1);
A=6*m/h +3*c;
B=3*m+h*c/2;
a(1)=F(1)/m;
%Flang 0 to indicate elastic response and flag 1 for plastic
flag=0; %Starting with elastic response
%This for loop is for calculation of elasto-plastic response
for i=1:length(F)-1
DF=F(i+1)-F(i);
if flag==0
   kp=k;
else
    kp=0;
end
Kt=kp+A/h;
dy=(DF+A*v(i)+B*a(i))/Kt;
u(i+1)=u(i)+dy;
v(i+1)=(3/h)*dy -h*a(i)/2 - 2*v(i);
if v(i+1)>0&&u(i+1)<yt
    flag=0;
    R(i+1)=k*u(i+1);
end
if v(i)>0&&v(i+1)<0
    ymax=u(i);
end
if u(i+1)>yt&&v(i+1)>0&&i<2000
    flag=1;
    R(i+1)=Rt;
end
if v(i+1)<0&&u(i+1)>(ymax-2*yt)
    flag=0;
    R(i+1)=Rt-k*(ymax-u(i+1));
end
if v(i+1)<0&&u(i+1)<(ymax-2*yt)
    flag=-1;
    R(i+1)=Rc;
end
if v(i)<0&&v(i+1)>0&&i<2000
    ymin=u(i);
end
if ymin<1000
    flag=0;
    R(i+1)=(u(i+1)-ymin)*k+Rc;
end
if v(i+1)>0&&u(i+1)>(ymin+2*yt)
    flag=1;
    R(i+1)=Rt;
end
a(i+1)=(F(i+1)-v(i+1)*c-R(i+1))/m; %Equilibrium condition checked for i+1th timestep
end
%%Following code is for elastic response of system
u1=zeros(1+round(restime/h),1);
v1=zeros(1+round(restime/h),1);
a1=zeros(1+round(restime/h),1);
Kt=k+A/h;
for i=1:round(restime/h)-1
  DF=F(i+1)-F(i);
dy=(DF+A*v1(i)+B*a1(i))/Kt;
u1(i+1)=u1(i)+dy;
v1(i+1)=(3/h)*dy -h*a1(i)/2 - 2*v1(i);
a1(i+1)=(F(i+1)-v1(i+1)*c-k*u1(i+1))/m;
end
%Code below this is for plotting of responses in suitable format
t=0:0.001:5;
figure(1)
plot(t,u, LineWidth=2);
hold on
plot(t,u1,LineWidth=1);
xlabel("Time in sec",FontSize=15);
ylabel("Displacement in inches", Fontsize=15);
title("Linear and non-linear response of system",FontSize=20);
figure(2)
plot(u,R,LineWidth=3);
xlabel("Displacement in inches",FontSize=15);
ylabel("Elastic force in Kips", Fontsize=15);
title("Loading Hysterisis loop",FontSize=20);
