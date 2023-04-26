h=0.005;%Timestep for successive computations
k=2800;%Stiffness value in kip/ft
W=100;%Weight of vibrating mass tank
x=0.05;%damping factor 
m=W/32.174;%mass of vibrating structure g=32.17ft/sec^2
c=2*m*x*(k/m)^0.5;%Value of damping coefficient in kip-sec/ft
F=100;%Value of Fo given in question
td=0.05;%Duration of ramp load
t=0.5;%Time for which response has to be calculated
%Varriable declaration for computation
n=round(t/h);
n1=round(td/h);
rem=zeros(1,n-n1);
u=zeros(1,n);
v=zeros(1,n);
u(1)=0;
v(1)=0;
t1=0:h:td;
p=F*(td-t1)/td;% force equation (Ramp) equation
P=[p,rem];%assigning general force matrix for each calculation points
u(2)=u(1)+F(1)/m *0.5*h*h;%Calculation for displacement at 2nd time step 
%Coefficients used in recurrance formula
M=(m/h^2 +c/(2*h));
A=(2*m/h^2 -k)/M;
B=(c/(2*h)-m/h^2)/M;
C=1/M;
D=1/(2*h);
%Main calculation of each time step displacement and velocity
for i=2:n
u(i+1)=u(i)*A+u(i-1)*B-P(i)*C;
v(i)=(u(i+1)*D-u(i-1)*D);
end
t=0:h:t;
%Code below this line is just for plotting of displacement and velocity
%time history multiplied by 12 to convert to inch
figure(1);
plot(t,u*12);
xlabel("Time in sec");ylabel("Displacement in inch");title("Time History Plot of Displacement");
v(n+1)=0;
figure(2);
plot(t,v*12);
xlabel("Time in sec");ylabel("Velocity in in/sec");title("Time History Plot of Velocity");