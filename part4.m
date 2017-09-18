rIn=1;      % The inner radius of the cylinder in centimeters
rOut=2;     % Outer radius of the cylinder in centimeters
uLiq=450;   % Temperature of the fluid
uAmb=20;    % Ambient temperature
n=25;       % N=25 (according to project description)
margin=0.1;   % Assumption: the difference between the temperature of two runs
rCount=0;   % Count loops
uInitLoop=0;
k=1;        % Heat transfer between metal and air
diff=100;   % Dummy 

while diff >= margin
    rCount=rCount+1;
    r=linspace(rIn, rOut, n);
    h=(rOut-rIn)/n;  % Discretization according to the finite differential method

    % Inhomogeneous metal
    d=1+(3*r.^2-9.*r+6)./(r.^3-4.5*r.^2+6.*r-3);
    
    % Boundary conditions
    col=zeros(n,1);
    col(1,1)=450;
    col(n,1)=20*k*((2*r(n)/h)+d(n));

    % Matrix
    mat=zeros(n,n);
    mat(1,1)=1;
    for i=2:n-1
        mat(i,i+1)=(r(i)/h^2)+d(i)/(2*h); 
        mat(i,i)=-2*r(i)/h^2; 
        mat(i,i-1)=r(i)/h^2-d(i)/(2*h);
    end
    mat(n,n-1)=-2*r(n)/h^2;
    mat(n,n)=k*d(n)+2*r(n)/h^2+2*k*r(n)/h;

    % Solve matrix
    u=mat\col;
    nLoop(rCount)=n;
    uLoop(rCount)=u(n);
    
    % Update conditions
    n=n*2;
    diff=abs(uLoop(rCount)-uInitLoop);
    uInitLoop=uLoop(rCount);
end

% Visual plot to show temperature 
plot(r,u);
grid on;
axis([1 2 20 uLiq]);
xlabel('Radius of the cylinder (cm)');
ylabel('Temperature of the metall');

