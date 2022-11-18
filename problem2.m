function problem2(tpr, N)
% First order Godunov, solving the burger's equation:
%
%     p_t + (p^2/2)_x = 0
%
% with periodic boundary conditions.

% Number of grid cells, dx and dt:
x = linspace(-pi/2,pi/2,N); dx=abs(x(2)-x(1)); dt=dx/4;
% problem # for choosing initial conditions
% problem = 2;
% Printing times:

npr=size(tpr,2);

% Number of time intervals in between printing times:
Ndt=round((tpr(2:npr)-tpr(1:npr-1))/dt);
L = [1 1:N-1]; M = 1:N; R = [2:N N];
% Auxiliary vectors for looking to the left, middle and right:
% Choose initial conditions and indices based on problem #
%         L = 1:N-2; M = 2:N-1; R = 3:N; 
initialCond = -1;
% Initial data:
p = -sin(x);
        
% Colors for the output:
% col=['k' 'b' 'r' 'g' 'y' 'c']; nc=length(col);

% Adjusted printing times, since the original ones are typically not multiples of dt.
for jpr=2:npr
  tpr(jpr)=tpr(jpr-1)+Ndt(jpr-1)*dt;
end



% Plot initial data.
% figure(1)
% plot(x,p,col(1)); hold on

% Main loop over printing times.
for jpr=2:npr

% Loop over time steps.
  for j = 1:Ndt(jpr-1)
  p(1) = -1*initialCond;p(end) = 1*initialCond;
% Characteristic speed.
    c = p;

% Fluxes Q: from the right, from the left, 
% or at x/t=0 within a rarefaction, yielding c=0 => p=1/2 => Q=1/4.
    Fl = p.^2/2;
    
    Q = (c(L)>0).*(c(M)>0).*Fl(L) + ...
        (c(L)<=0).*(c(M)<=0).*Fl(M) + ...
        (c(L)>0).*(c(M)<=0).*(((c(L)+c(M))>0).*Fl(L)+((c(L)+c(M))<=0).*Fl(M)) +...
        (c(L)<=0).*(c(M)>0)*0.25;

% (The third line corresponds to a shock moving with either positive or
% negative speed.)

% Update of p using conservation.
    p = p + dt*(Q(M)-Q(R))/dx;
%     p(1) = -1*initialCond;p(end) = 1*initialCond;
  end

% Plot at printing time jpr.
  plot(x,p,'o', 'DisplayName',"t = "+num2str(tpr(jpr)));hold on;

end
legend
xlabel('x'); ylabel('u'); 
title("Snapshots of a numerical solution to the Burger's equations")
hold off
%%



