clc ;
clf ; 
result = load('var_u_5000.dat') ;
nx = 100;
ny = 11 ; 
dx = 0.3 ;
dy = 0.1 ;
y = (0:1:ny-1)*dy ;
% x = (0:1:nx-1)*dx ; 
x = 1 ;
x_0  = get_velocity_x(result,1,ny);
x_4 = get_velocity_x(result,5,ny);
x_70 = get_velocity_x(result,71,ny);
hold on
plot(y,x_0(:),'k','LineWidth',1)
plot(y,x_4(:),'--k','LineWidth',1)
plot(y,x_70(:),'r','LineWidth',1)
plot_ideal()
hold off
legend('at x = 0','at x = 1.5','at x = 21','Ideal','fontsize',14,'interpreter','Latex')
ylabel('Velocity\hspace{0.1cm}in\hspace{0.1cm}x\hspace{0.1cm}axis\hspace{0.1cm}(u)','fontsize',14,'interpreter','Latex')
xlabel('Channel\hspace{0.1cm}Height','fontsize',14,'interpreter','Latex')
title('Velocity\hspace{0.1cm}Profile','fontsize',14,'interpreter','Latex')
grid;
function target = get_velocity_x(result,x,ny)
    if (x == 1)
        target = zeros(ny,1);
        target(2:ny-1,1) = result(3:ny,x);
        target(1,1) = 1 ; 
        target(end,1) = 1; 
    else
    target = zeros(ny,1);
    target(2:ny-1,1) = result(3:ny,x);
    end
end
function plot_ideal() 
    a = 1; %distance between plate
    neu = 1; %dynamic viscosity
    p1 = 100 ; %pressure at inlet 
    p2 = 96.45 ; %pressure at outlet 
    dp = p2-p1;

    x1 = 0;
    x2 = 0.3;
    dx = x2-x1; %plate length

    y = linspace(0,a,100);

    u = a^2/(2*neu)*(dp/dx)*((y/a).*(y/a)-(y/a)) ;
    plot(y,u,'b','LineWidth',2);
end
