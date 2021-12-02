function [xDot]=snakeEquationsHeadTracking(t,x,param)

% points following + broken joints 
% robot parameters
ct = param.ct; 
cn = param.cn;
N  = param.N;
m  = param.m;
l  = param.l;
J  = m*l*l/3;

% controler parameters
Kp = param.Kp;
Kd = param.Kd;

% state vector
fi    = x(1:N);
p     = x(N+1:N+2);
fiDot = x(N+3:2*N+2);
pDot  = x(2*N+3:2*N+4);

% equation 2.37
H = ones(N,N);
H = triu(H);
theta = H*fi;
thetaDot = H*fiDot;

% auxiliary matrices A, D, e, Stheta, CTheta 
% from definition on page 41
A      = [diag(ones(N-1,1)) zeros(N-1,1)]+[zeros(N-1,1) diag(ones(N-1,1))];
D      = [diag(ones(N-1,1)) zeros(N-1,1)]-[zeros(N-1,1) diag(ones(N-1,1))];
e      = ones(N,1);
%a      = [1; zeros(N-1,1)]; % added 
a      = [zeros(N-1,1); 1]; % added 
Stheta = diag(sin(theta));
Ctheta = diag(cos(theta));

% equation 2.34c and 2.34d
K      = A'*((D*D')^(-1))*D;
V      = A'*((D*D')^(-1))*A;
Gp     = a'-a'*D'*((D*D')^(-1))*D; % added
Vp     = a'*D'*((D*D')^(-1))*A; % added
%C      = tril(-1*ones(N,N-1),-1); % added
C      = triu(ones(N,N-1));

% equation 2.34a and 2.34b
W   = m*l^2*Stheta*V*Ctheta-m*l^2*Ctheta*V*Stheta;
M   = J*eye(N)+m*l^2*Stheta*V*Stheta+m*l^2*Ctheta*V*Ctheta;

Wp  =  m*l*[Vp*Ctheta; Vp*Stheta];  % added
Mp  =  m*l*[Vp*Stheta; -Vp*Ctheta];  % added

% control algorithm
u = zeros(N-1,1);
%fi_zad=fi_zad_agg{1};
%fi_zad_dot=fi_zad_agg{2};
%fi_zad_dotdot=fi_zad_agg{3};

T1=[D; a'];

X=T1^(-1)*[-l*A*cos(theta); p(1)];
Y=T1^(-1)*[-l*A*sin(theta); p(2)];

%X=l*C*A*cos(theta)+e*p(1);
%Y=l*C*A*sin(theta)+e*p(2);

global fi0_last;
global flag;

for i=1:size(flag,2)
    %robot has not reached the point
    if ~flag(i)
        
        rx = param.Px(i)-p(1);
        ry = param.Py(i)-p(2);
         %rx = param.Px(i)-sum(X)/param.N;
         %ry = param.Py(i)-sum(Y)/param.N;
        
        r  = sqrt(rx^2+ry^2);
        if r>0.05
            theta_bar = (sum(theta)/N);
            %theta_bar=theta(1);
            theta_ref = -atan2(ry,rx);
            fi0  = param.st1*(theta_bar+theta_ref);
             if fi0 > deg2rad(40)
                 fi0=deg2rad(40);
             elseif fi0 < -deg2rad(40)
                 fi0=-deg2rad(40);
             end
            fi0_last=fi0;
        else
            flag(i)=t;
            fi0 = fi0_last;
        end
        break
    else
        fi0 = 0;
    end
end

for i = 1:N-1
    fi_z = param.ampli_u*sin(param.omega_u*t+(i-1)*param.delta_u); 
    fi_z = fi_z+fi0;
    u(i) = Kp*(fi_z-fi(i))-Kd*fiDot(i); %control law
end

if param.coupledControl
    for i = 1:N-1
        if i==1
            um(1)=u(1);%-u(i+1);
        %elseif i==(N-1)
        %    um(N-1)=u(N-1);%+u(i-1);
        else
            um(i)=u(i)+u(i-1)/2;%-u(i+1));
        end
    end

    u=um';
end

% Robot disturbances
% if param.randDist
%     u(param.distLink)=u(param.distLink).*(1+(param.disturb-0.5)*param.errLevel);
% else
    u(param.distLink)=0;
%end

% equation 2.12a and 2.12b
%XDot =  l*K'*Stheta*thetaDot+e*pDot(1);
%YDot = -l*K'*Ctheta*thetaDot+e*pDot(2);

XDot =  l*C*A*Stheta*thetaDot+e*pDot(1); % changed
YDot = -l*C*A*Ctheta*thetaDot+e*pDot(2); % changed

% equation 2.25
fR  = -[ct*(Ctheta.^2)+cn*(Stheta.^2) (ct-cn)*Stheta*Ctheta; ...
    (ct-cn)*Stheta*Ctheta ct*(Stheta.^2)+cn*(Ctheta.^2)]*[XDot;YDot];

%hx=inv(D*D')*(m*l*A*())
%disp(D'*((D*D')^(-1))*D);

% equation 2.39 in two parts - 
% 2.41a i 2.41b - actuated i unactuated joints
%MFi = [H'*M*H, zeros(N,2); ...
%      zeros(2,N), N*m*eye(2,2)];  
MFi = [H'*M*H, zeros(N,2); ...
      -Mp*H, N*m*eye(2,2)]; % changed
  

% MFi_11 = MFi(1:N-1, 1:N-1);
% MFi_12 = MFi(N:N+2, 1:3);
MFi_21 = MFi(N:N+2, 1:N-1);
MFi_22 = MFi(N:N+2, N:N+2);

%WFi = [H'*W*diag(H*fiDot)*H*fiDot; zeros(2,1)];
WFi = [H'*W*diag(H*fiDot)*H*fiDot;...
    -Wp*diag(H*fiDot)*H*fiDot]; % changed

% WFi_1 = WFi(1:N-1);
WFi_2 = WFi(N:N+2);

% equation 2.40d
%GFi = [-l*H'*Stheta*K, l*H'*Ctheta*K; ...
%      -e', zeros(1,N); ...
%      zeros(1,N), -e'];
GFi = [-l*H'*Stheta*K, l*H'*Ctheta*K; ...
      -Gp, zeros(1,N); ...
      zeros(1,N), -Gp]; % changed


% GFi_1 = GFi(1:N-1, 1:2*N);
GFi_2 = GFi(N:N+2, 1:2*N);

% equation 2.46a i 2.46b
AqFi = -MFi_22^(-1)*(WFi_2 +GFi_2*fR);
BqA  = -MFi_22^(-1)*MFi_21;

% State vector from equation 2.47
xDot=[fiDot; pDot; u; AqFi+BqA*u];