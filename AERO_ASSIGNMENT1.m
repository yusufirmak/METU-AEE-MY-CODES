numP= 4
numBP= numP+1 %number of boundary points;
tO=(360/(numP))/2;
flagPlot= [1;
           1;
           1;
           1];

theta=linspace(0,360,numBP)';
theta= theta+tO;
theta=theta*(pi/180);

xb=cos(theta)
yb=sin(theta)

edge= zeros(numP,1);
for i=1:1:numP
    edge(i)= (xb(i+1)-xb(i))*(yb(i+1)+yb(i));
end
sumEdge= sum(edge);

if (sumEdge < 0)                                        
    xb = flipud(xb);                                                       
    yb = flipud(xb); 
end

xc= zeros(numP,1);
yc= zeros(numP,1);
s= zeros(numP,1);
phiD= zeros(numP,1);


for i = 1:1:numP
        xc(i)   = 0.5*(xb(i)+xb(i+1));
        yc(i)   = 0.5*(yb(i)+yb(i+1));  
        dx      = xb(i+1)-xb(i); 
        dy      = yb(i+1)-yb(i);  
        S(i)    = (dx^2 + dy^2)^0.5;
        	phiD(i) = atan2d(dy,dx);
        if (phiD(i) < 0)                                                       
        phiD(i) = phiD(i) + 360;
        end    
end

deltaD             = phiD + 90; 
betaD              = deltaD;
betaD(betaD > 360) = betaD(betaD > 360) - 360; 

phi  = phiD.*(pi/180);
beta = betaD.*(pi/180);

[I,J] = COMPUTE_IJ_SPM(xc,yc,xb,yb,phi,S);                                  % Compute geometric integrals
A = zeros(numP,numP);                                                   % Initialize the A matrix
for i = 1:1:numP                                                          % Loop over all i panels
    for j = 1:1:numP                                                     % Loop over all j panels
        if (i == j)                                                         % If the panels are the same
            A(i,j) = pi;                                                    % Set A equal to pi
        else                                                                % If panels are not the same
            A(i,j) = I(i,j);                                                % Set A equal to geometric integral
        end
    end
end

b = zeros(numP,1);                                                        % Initialize the b array
for i = 1:1:numP                                                          % Loop over all panels
    b(i) = -Vinf*2*pi*cos(beta(i));                                         % Compute RHS array
end

lambda  = A\b;                                                              % Compute all source strength values
fprintf('Sum of L: %g\n',sum(lambda.*S)); 

Vt = zeros(numP,1);                                                       % Initialize tangential velocity array
Cp = zeros(numP,1);                                                       % Initialize pressure coefficient array
for i = 1:1:numP                                                          % Loop over all i panels
    addVal  = 0;                                                            % Reset the summation value to zero
    for j = 1:1:numP                                                      % Loop over all j panels
        addVal = addVal + (lambda(j)/(2*pi))*(J(i,j));                      % Sum all tangential source panel terms
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal;                                     % Compute tangential velocity by adding uniform flow term
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient
end

analyticTheta = linspace(0,2*pi,200)'                                      % Analytical theta angles [rad]
analyticCP    = 1-4*sin(analyticTheta).^2 

if (flagPlot(1) == 1)
    figure(1);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    fill(XB,YB,'k');                                                        % Plot polygon
    for i = 1:1:numPan                                                      % Loop over all panels
        X(1) = XC(i);                                                       % Set X start of panel orientation vector
        X(2) = XC(i) + S(i)*cosd(betaD(i)+AoA);                             % Set X end of panel orientation vector
        Y(1) = yc(i);                                                       % Set Y start of panel orientation vector
        Y(2) = yc(i) + S(i)*sind(betaD(i)+AoA);                             % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',3);                                       % Plot panel normal vector
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

if (flagPlot(2) == 1)
    figure(2);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    plot(XB,YB,'k-','LineWidth',3);                                         % Plot panels
    p1 = plot([XB(1) XB(2)],[YB(1) YB(2)],'g-','LineWidth',3);              % Plot first panel
    p2 = plot([XB(2) XB(3)],[YB(2) YB(3)],'m-','LineWidth',3);              % Plot second panel
    pB = plot(XB,YB,'ko','MarkerFaceColor','k','MarkerSize',10);            % Plot boundary points
    pC = plot(XC,yc,'ko','MarkerFaceColor','r','MarkerSize',10);            % Plot control points
    legend([pB,pC,p1,p2],...                                                % Show legend
           {'Boundary','Control','First Panel','Second Panel'});
    ylim([-1 1]);                                                           % Set Y-limits
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    zoom reset;                                                             % Reset zoom
end

if (flagPlot(3) == 1)
    figure(3);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    pA = plot(analyticTheta,analyticCP,'k-','LineWidth',3);                 % Plot analytical pressure coefficient
    pC = plot(beta,Cp,'ks','MarkerFaceColor','r','MarkerSize',10);          % Plot compute pressure coefficient
    xlabel('Angle [rad]');                                                  % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 2*pi]);                                                         % Set X-limits
    ylim([-3.5 1.5]);                                                       % Set Y-limits
    legend([pA,pC],{'Analytical','SPM'},'Location','S');                    % Add legend
end

if (flagPlot(4) == 1)
    figure(5);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    axis equal;                                                             % Set axes equal
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    contourf(XX,YY,CpXY,100,'EdgeColor','none');                            % Plot contour
    fill(XB,YB,'k');                                                        % Plot polygon
    xlim(xVals);                                                            % Set X-limits
    ylim(yVals);                                                            % Set Y-limits
    zoom reset;                                                             % Reset zoom
end
