clc ;
result = load('var_u_5000.dat') ;
nx = 100;
ny = 11 ; 
dx = 0.3 ;
dy = 0.1 ;
y = (0:1:ny-1)*dy ;
x = (0:1:nx-1)*dx ; 
thickness = zeros(length(x)) ;
X = [0:1:nx]*dx ; 
for i = 1:length(x) 
    x_p = x(i);  
    target = get_velocity_x(result,i,ny);
    f = spline(y,target);
    dfdx = fnder(f);
    thickness(i) = get_y_lowest_dev(dfdx,ny,dy) ;
end
plot(X(1:30),thickness(1:30),'b','LineWidth',2)
hold on
yline((ny-1)*dy/2,'--',{'Maximum\hspace{0.1cm}thickness'},'fontsize',14,'interpreter','Latex','LineWidth',2)
xline(5.7,'--',{'$L_{f}$ = 5.7'},'fontsize',14,'interpreter','Latex','LineWidth',2)
ylim([0,0.6])
ylabel('Boundary\hspace{0.1cm}layer\hspace{0.1cm}thickness\hspace{0.1cm}$(\delta)$','fontsize',14,'interpreter','Latex')
xlabel('Channel\hspace{0.1cm}Distance\hspace{0.1cm}(x)','fontsize',14,'interpreter','Latex')
title('Velocity\hspace{0.1cm}boundary\hspace{0.1cm}layer\hspace{0.1cm}thickness\hspace{0.1cm}as\hspace{0.1cm}a\hspace{0.1cm}function\hspace{0.1cm}of\hspace{0.1cm}x'...
        ,'fontsize',14,'interpreter','Latex')
hold off
function target = get_velocity_x(result,x_p,ny) 
    target = zeros(ny,1);
    target(1:ny,1) = result(2:ny+1,x_p);
end
function y_lowest = get_y_lowest_dev(dfdx,ny,dy)
    Range = [0:0.01:ny]*dy ;
    dev = ppval(dfdx,Range) ;
    idx=find(dev<0.01&dev>-0.01) ;
    target = idx(1) ; 
    y_lowest = Range(target);
end
function plot_fn(result,x_p,ny,y) 
    target = get_velocity_x(result,x_p,ny);
    f = spline(y,target);
    dfdx = fnder(f);
    fnplt(dfdx) ;
end
