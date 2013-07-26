
T = 10;
N = 10000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 28000;
rho = .004; %log(1+rate)
beta = 9; %.003*3000
g = .04;

% Set up the parameters for the deer simulation
alpha = 0.03;
gamma = 0.05;
P     = 277500.0;

rtilde=r1-h;
ftilde=(rtilde/r1)*F;
a=(alpha^2)-rtilde;      
b=(.5*(alpha^2))-a;
g0= 1-(rtilde/b);

% Variables for the bond fund
m0 = (beta*ftilde - P)/(rho-g);
rho = .004; %log(1+rate)
beta = 9; %.003*3000
g = .04;
m = m0;

dt=sqrt(T/N);

clf;

%Deer Population on 1
subplot(2, 1, 1);
axis([0 T ftilde*0.9 ftilde*1.1]);
hold on;

%Bond funds on 2
subplot(2, 1, 2);
axis([0 T 0 2*m0]);
hold on;



z = 1.0;
X = ftilde;
W = 0.0;
Q = 0.0;
for j=0:N,
    dW = dt*randn(1,1);
    s = j*dt;
    
    % Update the deer population
    Q = Q + exp(b*s+alpha*W).*dW;   %integral part of analytic solution
    Z = (rtilde/b) + g0*exp((-b*s)-(alpha*W)) + exp((-b*s)-(alpha* ...
                                                      W)).*((-alpha*rtilde)/b).*Q; %analytic solution of Z
    prevX = X;
    X = (ftilde./Z); % transform back to X
    
    % Update the bond fund
    m = m + dt*(rho*m - beta*prevX+P) - dW*gamma*prevX - .5*gamma*alpha*prevX*(dW*dW-dt);

    % Plot the updates
    subplot(2, 1, 1);
    plot(s,X,'rx');
    subplot(2, 1, 2);
    plot(s,m,'bo');


    W = W + dW;
    
    drawnow();
    
end
