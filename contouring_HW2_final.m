clear; close all; clc
nex = 10;%input('Give the number of elements in x direction: ');
ney = 10;%input('Give the number of elements in y direction: ');
theta1 = 90;
theta2 = 60;
theta3 = 30;

nnx = 2*nex+1;
nny = 2*ney+1;

%% Results for different angle theta

%Theta = 90 deg
[xpt_90 ypt_90 u_90 sk_90 r1_90]=FEM_2d_HW2_final(nex,ney,theta1);
[xi_90, yi_90] = meshgrid(linspace(min(xpt_90),max(xpt_90),length(xpt_90)),linspace(min(ypt_90),max(ypt_90),length(ypt_90)));
zi_90 = griddata(xpt_90,ypt_90,u_90,xi_90,yi_90);
%Theta = 60 deg
[xpt_60 ypt_60 u_60 sk_60 r1_60]=FEM_2d_HW2_final(nex,ney,theta2);
[xi_60, yi_60] = meshgrid(linspace(min(xpt_60),max(xpt_60),length(xpt_60)),linspace(min(ypt_60),max(ypt_60),length(ypt_60)));
zi_60 = griddata(xpt_60,ypt_60,u_60,xi_60,yi_60);
%Theta = 30 deg
[xpt_30 ypt_30 u_30 sk_30 r1_30]=FEM_2d_HW2_final(nex,ney,theta3);
[xi_30, yi_30] = meshgrid(linspace(min(xpt_30),max(xpt_30),length(xpt_30)),linspace(min(ypt_30),max(ypt_30),length(ypt_30)));
zi_30 = griddata(xpt_30,ypt_30,u_30,xi_30,yi_30);

%% Ploting the results
n1=nnx;
n2=nny;
l1=1;
l2=1;
k1=1;
k2=1;
for i = 1:nnx
    
    x_90(1:nnx,k1)=xpt_90(l1:n1,1);
    x_60(1:nnx,k1)=xpt_60(l1:n1,1);
    x_30(1:nnx,k1)=xpt_30(l1:n1,1);
    n1=n1+nny;
    k1=k1+1;
    l1=l1+nny;
    
end

for j=1:nny
    
    y_90(1:nny,k2)=ypt_90(l2:n2,1);
    y_60(1:nny,k2)=ypt_60(l2:n2,1);
    y_30(1:nny,k2)=ypt_30(l2:n2,1);
    n2=n2+nny;
    k2=k2+1;
    l2=l2+nny;
    
end

z_90=zeros(nnx,nny);
z_60=zeros(nnx,nny);
z_30=zeros(nnx,nny);

%% Theta = 90 deg
figure(1)
mesh(x_90,y_90,z_90,'edgecolor','k')
title('2D-Mesh (90 deg)')
xlabel('x')
ylabel('y')

figure (2)
contourf(xi_90,yi_90,zi_90)
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Contour plot Temperature Distribution (90deg)') 

figure (3)
surf(xi_90,yi_90,zi_90,'edgecolor','none')
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Surf plot Temperature Distribution (90deg)') 

%% Theta = 60 deg
figure(4)
mesh(x_60,y_60,z_60,'edgecolor','k')
title('2D-Mesh (60deg)')
xlabel('x')
ylabel('y')

figure (5)
contourf(xi_60,yi_60,zi_60)
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Contour plot Temperature Distribution (60deg)') 

figure (6)
surf(xi_60,yi_60,zi_60,'edgecolor','none')
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Surf plot Temperature Distribution (60deg)') 

%% Theta = 30 deg
figure(7)
mesh(x_30,y_30,z_30,'edgecolor','k')
title('2D-Mesh (30deg)')
xlabel('x')
ylabel('y')

figure (8)
contourf(xi_30,yi_30,zi_30)
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Contour plot Temperature Distribution (30deg)') 

figure (9)
surf(xi_30,yi_30,zi_30,'edgecolor','none')
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('T - temperature (C)')
title('Surf plot Temperature Distribution (30deg)') 


