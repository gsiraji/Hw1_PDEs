% Main 

% choose problem to show the corresponding figure for
problem  = 1; % 1 or 2
% choose printing time and grid size,...
% endpoints for problem1 
tprint = [1 2 3 4];Ngrid =200;endpoint=5;


switch problem
    case 1
        problem1(tprint,Ngrid,endpoint)
    case 2
        problem2(tprint,Ngrid)
end