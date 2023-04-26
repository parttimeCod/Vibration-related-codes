function [nod_disp] = beamsolver(bl,E,I,nel,P,x)
%This function gives the nodal displacemenent as output by taking beam
%parameters as input
l=bl/nel;
k=(E*I/l^3)*[12,6*l,-12,6*l;6*l,4*l*l,-6*l,2*l*l;-12,-6*l,12,-6*l;6*l,2*l*l,-6*l,4*l*l];
K_Str=zeros((nel+1)*2,(nel+1)*2);
F=zeros((nel+1)*2,1);
for i=1:nel
    K_Str(2*i-1:2*i+2,2*i-1:2*i+2)=K_Str(2*i-1:2*i+2,2*i-1:2*i+2)+k;
end
for i=1:nel
    if x>=(i-1)*l && x<i*l
        F(2*i-1,1)=F(2*i-1,1)-P*(i*l-x)/l;
        F(2*i,1)=F(2*i,1)-P*((i*l-x)^2)*(x-(i-1)*l)/l^2;
        F(2*i+1,1)=F(2*i+1,1)-P*((x-(i-1)*l))/l;
        F(2*i+2,1)=F(2*i+2,1)+P*((i*l-x))*(((x-(i-1)*l))^2)/l^2;
    end
end
%Applying boundary condition as simply supported
k_cp=K_Str;
k_cp(nel*2+1,:)=[];
k_cp(:,nel*2+1)=[];
k_cp(1,:)=[];
k_cp(:,1)=[];
f_final=F;
f_final(nel*2+1,:)=[];
f_final(1,:)=[];
nod_disp=k_cp\f_final;

end