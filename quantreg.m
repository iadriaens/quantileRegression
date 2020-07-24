function [ bhat,yhat,res ] = quantreg( y, x, tau, order, optie)
%   bhat are the estimates
%   y is a vector of outcomes
%   x is a vector of explanatory variables
%   tau is a scalar for choosing the conditional quantile to be estimated
%   order is the order of the polynomial

if size(x,2)~=1
    x = x';     % set x as a column vector if row
end

% construct X so order = order
x = [ones(length(x),1) x];
if order > 1
    for i = 1:order
        x = [x x(:,2).^i];
    end
end


n=size(x,1);        % number of rows of X
m=size(x,2);        % number of columns of X = polynomial grade

if optie == 1
    % vectors and matrices for linprog
    f=[tau*ones(n,1);(1-tau)*ones(n,1);zeros(m,1)];
elseif optie == 2
    % vectors and matrices for linprog - more weight in beginning for
    % negative res
    f=[0.10*ones(15,1);tau*ones(n-15,1);(1-0.10)*ones(15,1);(1-tau)*ones(n-15,1);zeros(m,1)];
end

Aeq=[eye(n),-eye(n),x];
beq=y;
lb=[zeros(n,1);zeros(n,1);-inf*ones(m,1)];
ub=inf*ones(m+2*n,1);

% Solve the linear programme
opts = optimoptions('linprog','Display','off');
[bhat,~,~]=linprog(f,[],[],Aeq,beq,lb,ub,opts);

% Pick out betas from (u,v,beta)-vector.
bhat=bhat(end-m+1:end);

% Calculate the yhat
yhat = x*bhat;

% Calculate the residuals
res = y-yhat;

end