  clear,clc,close all
  a=1; c=5; theta_slope(1)=1; phi_slope(1)=0; %intitial condition  
  %(theta_slope,phi_slope)=(1,0) gives equator and (0,1) gives meridian as Geodesic
  phi(1)=0; theta(1)=0;  %initial point.  
  h = 0.01;% step size   
 
  %Solving using Runge-Kutta Method
    for i=1:1:628   % calculation loop
        phi(i+1) = phi_slope(i)*h + phi(i); theta(i+1) = theta_slope(i)*h + theta(i);
        G = @(p,tslope,pslope) ((2*a*sin(p)*tslope*pslope)/(c+a*cos(p)));
        F = @(p,tslope) ((-1/a)*(sin(p))*(c+(a*cos(p)))*(tslope^2));
    l_1 = h*G(phi(i),theta_slope(i),phi_slope(i));
    k_1 = h*F(phi(i),theta_slope(i));
    l_2 = h*G(phi(i),theta_slope(i)+(0.5*l_1),phi_slope(i)+(0.5*k_1));
    k_2 = h*F(phi(i),theta_slope(i)+(0.5*l_1));
    l_3 = h*G(phi(i),theta_slope(i)+(0.5*l_2),phi_slope(i)+(0.5*k_2));
    k_3 = h*F(phi(i),theta_slope(i)+(0.5*l_2));
    k_4 = h*F(phi(i),theta_slope(i)+l_3);
    l_4 = h*G(phi(i),theta_slope(i)+l_3,phi_slope(i)+k_3);
    theta_slope(i+1) = theta_slope(i) + (1/6)*(l_1+2*l_2+2*l_3+l_4)*h;% main equation
    phi_slope(i+1) = phi_slope(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    end
%Plotting the Geodesic    
x1 = (c+a*cos(phi)).*cos(theta) ;
y1 = (c+a*cos(phi)).*sin(theta) ;
z1 = a*sin(phi) ;
z2 = z1';
X1 = reshape(x1,[],1);
Y1 = reshape(y1,[],1);
Z1 = reshape(z1,[],1);
plot3(X1,Y1,Z1);
hold on
grid on
%Plotting Torus
u=linspace(0,2*pi);
v1=linspace(0,2*pi);
[U,V] = meshgrid(u,v1) ;
X2 = (c+a*cos(V)).*cos(U) ;
Y2 = (c+a*cos(V)).*sin(U) ;
R = c+a*cos(V);
Z2 = a*sin(V) ;
surf(X2,Y2,Z2);
shading interp
