function problem1(tpr, N,a)
% First order Godunov, solving the traffic flow equation:
%
%     p_t + (p(1-p))_x = 0
%
% with periodic boundary conditions.


% N=200;a =4;
dx=a/N; dt=dx/10;

% Printing times:
npr=size(tpr,2);

% Number of time intervals in between printing times:
Ndt=round((tpr(2:npr)-tpr(1:npr-1))/dt);

% Auxiliary vectors for looking to the left, middle and right:
L = [1 1:2*N-1]; M = 1:2*N; R = [2:2*N 2*N]; % (Notice the periodicity here.)
% Adjusted printing times, since the original ones are typically not multiples of dt.
for jpr=2:npr
  tpr(jpr)=tpr(jpr-1)+Ndt(jpr-1)*dt;
end

% Position at the middle of each cell:
x1=dx/2:dx:a;
x2= -a:dx:-dx/2;
% Initial data:
p1 = x1;
p1(x1<-1) = -1;
p1(x1>1) = 1;
p2 = x2;
p2(x2<-1) = -1;
p2(x2>1) = 1;

% Plot initial data.
figure(1)
x = [x1 x2];p=[p1 p2];

% Main loop over printing times.
for jpr=2:npr

% Loop over time steps.
  for j = 1:Ndt(jpr-1)
      
% Characteristic speed.
    c = p;

% Fluxes Q: from the right, from the left, 
% or at x/t=0 within a rarefaction, yielding c=0 => p=1/2 => Q=1/4.
    Fl = p.^2./2;

    Q = (c(L)>0).*(c(M)>0).*Fl(L) + ...
        (c(L)<=0).*(c(M)<=0).*Fl(M) + ...
        (c(L)>0).*(c(M)<=0).*(((c(L)+c(M))>0).*Fl(L)+((c(L)+c(M))<=0).*Fl(M)) +...
        (c(L)<=0).*(c(M)>0)*0.25;
    
% (The third line corresponds to a shock moving with either positive or
% negative speed.)

% Update of p using conservation.
    p = p + dt*(Q(M)-Q(R))/dx;

  end

% Plot at printing time jpr.
  plot(x(1:end-1),p(1:end-1),'o','DisplayName',"t = "+num2str(tpr(jpr))); hold on;

end
legend
xlabel('x'); ylabel('u'); 
title("Snapshots of a numerical solution to the Burger's equations")
hold off


