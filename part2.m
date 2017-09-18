rIn=1;      % The inner radius of the cylinder in centimeters
rOut=2;     % Outer radius of the cylinder in centimeters
uLiq=450;   % Temperature of the fluid
uAmb=20;    % Ambient temperature
n=25;       % N=25 (according to project description)
margin=0.1;   % Assumption: the difference between the temperature of two runs
rCount=0;   % Count loops
uInitLoop=0;
uOuter=100;

r=linspace(rIn, rOut, n);
h=(rOut-rIn)/n;  % Discretization according to the finite differential method

% Boundary conditions
col=zeros(n,1);
col(1,1)=450;
col(n,1)=100;

% Matrix
mat=zeros(n,n);
mat(1,1)=1;
for i=2:n-1
    mat(i,i+1)=(r(i)/h^2)+1/(2*h); 
    mat(i,i)=-2*r(i)/h^2; 
    mat(i,i-1)=r(i)/h^2-1/(2*h);
end
mat(n,n)=1;

% Solve matrix
u=mat\col;
   
plot(r,u);
grid on;
axis([1 2 20 uLiq]);
xlabel('Radius of cylinder (cm)');
ylabel('Temperature of the metall');

% Matrix equation for row = n
k= (2*r(n)/h^2) * (1/(1+2*r(n)/h)) * (u(n-1)-u(n)) * (1/(u(n)-20));
