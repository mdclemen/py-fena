%% Example 4.4
% Linearization

%% Clean up
close all;
clear y true_solution lin_trap full_trap H N h yl yf;
clear STEPS LIN_ERR FULL_ERR;

%% Obtain true solution
y = dsolve('Dy + y*(1-y) = 0','y(0)=1/2')
true_solution = subs(y,'t',1);

%% Setup schemes
% Linearized Trapeziodal method
lin_trap = inline('y + h*y*(y-1)/(1-h*(y-1/2))','y','h');
full_trap = inline('(2/h+1-sqrt((2/h+1)^2-4*(2/h*y+y*(y-1))))/2','y','h');

%% Collect Data
H = round(logspace(0,2.7,10));
for i = 1:length(H)
    N = H(i);
    h = 1/N;
    yl = 1/2;
    yf = 1/2;
    for j = 1:N
        yl = lin_trap(yl,h);
        yf = full_trap(yf,h);
    end
    STEPS(i) = N;
    LIN_ERR(i) = abs(yl - true_solution);
    FULL_ERR(i) = abs(yf - true_solution);
end

%% Plot
figure(1); clf;
loglog(STEPS,LIN_ERR,'o-')
hold on
loglog(STEPS,FULL_ERR,'s-')
grid on
grid minor
xlabel('N -- Number of Steps','FontSize',14)
ylabel('Error','FontSize',14)
legend('Trapeziodal','Linearized Trapeziodal')
