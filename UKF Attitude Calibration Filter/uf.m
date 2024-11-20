function [qe,we,be,ge,s1e,pe,P_cov]=uf(dts,tf,y)

% function [qe,we,be,ge,s1e,P_cov]=uf(dts,tf,y)
% 
% Unscented Filter for Spacecraft System Alignment Calibration, Lai, Crassidis, Harman
%
% INPUT: -
% dts       sampling interval
% tf        final sampling time
% y         measurements, 3X1 gyros, 4X1 star sensor, 3X1 payload
%
% OUTPUT: -
% qe        estimated quaternion
% we        estimated spacecraft rotational rate
% be        estimated gyro bias
% ge        estimated gyro misalignments, scale factor error
% s1e       estimated star tracker 1 misalignment
% P_cov     output covariance

% constants
T=dts;                  % sampling interval
parameters;             % loading the sensors misalignments, noise properties ...
n=18; nn=2*n; % n1=n+1; n2=n+2; n3=n+3;
kap=3-n;    % doesn't work when it is negative, 
%kap=3;
par=1;
fac=2*(par+1);  % Multiplication factor for attitude error representation

% pre-allocate space
qe_his  = zeros(4,tf/dts);      % to store all the estimated quaternions, qe
we_his  = zeros(3,tf/dts);      % estimated rotational rate
be_his  = zeros(3,tf/dts);      % estimated gyro bias
ge_his  = zeros(6,tf/dts);      % estimated gyro misalignments, scale factor error
s1e_his = zeros(3,tf/dts);      % estimated star tracker 1 misalignment
pe_his  = zeros(3,tf/dts);      % estimated payload misalignment
P_cov   = zeros(n,tf/dts);     % covariance
qq      = zeros(4,nn);
dq      = zeros(4,nn);
yez_p   = zeros(3,nn);
yez_st1 = zeros(3,nn);

% initialization
qe=e2q(pi/4,pi/2,pi/3,7)'; be=[0;0;0]*pi/180/3600; %be=bg;    %%%%% !!!!!
qe=[0.05148900000000  -0.54903900000000  -0.83061600000000   0.07734550000000]';
ge=[0;0;0;0;0;0]; s1e=[0;0;0]; pe=[0;0;0]; %ge=g; s1e=s1; %%%%%!!!!!
xest=[zeros(3,1);be;ge;s1e;pe];                % initial error vector
we=T_b_g0*(y(1:3,1)-be);
Pa=(.5*pi/180)^2; Pb=(.5*pi/180/3600)^2; Pk=(500/1e6)^2; Pm=(500/60^2*pi/180)^2; Ps=(50/60^2*pi/180)^2; Pp=(60/60^2*pi/180)^2;
Pa=(5*pi/180)^2; Pb=(.5*pi/180/3600)^2; Pk=(500/1e6)^2; Pm=(500/60^2*pi/180)^2; Ps=(50/60^2*pi/180)^2; Pp=(60/60^2*pi/180)^2;
P = diag([Pa Pa Pa Pb Pb Pb Pm Pm Pm Pk Pk Pk Ps Ps Ps Pp Pp Pp]);
R = [SIG_s1pm zeros(3); zeros(3) SIG_pm];     % measurement covariance
P_cov(:,1)=diag(P);

wait = waitbar(0,'UF running ...');
for i=1:tf/dts-1
    wgm = y(1:3,i);
    Og  = omega_g(wgm-be);
    Me  = m_fun(ge,we);
    
    % Compute error quaternion and propagate it forward
    qq0 = om(we,dts)*qe;
    %we  = T_b_g0*[wgm + (eye(3)+Me)*be + Og*ge];
    %we  = T_b_g0*[(wgm-be)+Og*ge];
    we   = T_b_g0*[eye(3)+Me]*[wgm-be];
    
    % covariance propagation
    RR  = [ T_b_g0     zeros(3,3) zeros(3,6) zeros(3,3) zeros(3,3)
            zeros(3,3) eye(3,3)   zeros(3,6) zeros(3,3) zeros(3,3)
            zeros(6,3) zeros(6,3) eye(6,6)   zeros(6,3) zeros(6,3)
            zeros(3,3) zeros(3,3) zeros(3,6) eye(3,3)   zeros(3,3) 
            zeros(3,3) zeros(3,3) zeros(3,6) zeros(3,3) ones(3,3) ];
    %QQ  = [ T/4*SIG_a+T^3/12*SIG_r+T^3/12*Og*SIG_g*Og' T^2/4*SIG_r T^2/4*Og*SIG_g zeros(3,3) zeros(3,3)
    %        T^2/4*SIG_r                                T*SIG_r     zeros(3,6)     zeros(3,3) zeros(3,3)
    %        T^2/4*SIG_g*Og'                            zeros(6,3)  T*SIG_g        zeros(6,3) zeros(6,3)
    %        zeros(3,3)                                 zeros(3,3)  zeros(3,6)     T*SIG_s1   zeros(3,3)
    %        zeros(3,3)                                 zeros(3,3)  zeros(3,6)     zeros(3,3) T*SIG_p    ];
    QQ  = [  T*SIG_a+T^3/3*SIG_r+T^3/3*Og*SIG_g*Og'  -T^2/2*SIG_r  -T^2/2*Og*SIG_g zeros(3,3)  zeros(3,3)
            -T^2*SIG_r                                T*SIG_r       zeros(3,6)     zeros(3,3)  zeros(3,3)
            -T^2*SIG_g*Og'                            zeros(6,3)    T*SIG_g        zeros(6,3)  zeros(6,3)
             zeros(3,3)                               zeros(3,3)    zeros(3,6)     T*SIG_s1    zeros(3,3)
             zeros(3,3)                               zeros(3,3)    zeros(3,6)     zeros(3,3)  T*SIG_p   ];            
    Q   = RR*QQ*RR';
    psquare = chol((n+kap)*(P+Q))';    % Matrix square root
    sigv    = real([psquare -psquare]);
    xx0     = xest;
    xx      = sigv+kron(xest,ones(1,nn));

% Compute other quaternions and propagate them forward
for j = 1:nn,
    if j==1|j==2|j==3|j==n+1|j==n+2|j==n+3
        qerr  = dp2dq(xx(1:3,j),fac,par);
        qinit = q_mult(qerr,qe);
    else
        qinit = qe;
    end
    bex=xx(4:6,j); gex=xx(7:12,j); 
    Mex=m_fun(gex,we); %s1ex=xx(16:18,j);
    Ogx       = omega_g(wgm-bex);
    %wex       = T_b_g0*[(wgm-bex)+Ogx*gex]; % try to comment this
    wex       = T_b_g0*[eye(3)+Mex]*[wgm-bex]; % try to comment this
    qq(:,j)   = om(wex,dts)*qinit;
    dq(:,j)   = q_mult(qq(:,j),[-qq0(1:3);qq0(4)]); % find err quaternion, dq(:,j)=delta_q(qq0,qq(:,j));
    xx(1:3,j) = fac*dq(1:3,j)/(par+dq(4,j));
end
% Mean estimate -- note that xx0(1:3) is still zero from the reset.
xest = 1/(n+kap)*(kap*xx0+.5*sum(xx,2));
% Mean covariance
pp0  = kap*(xx0-xest)*(xx0-xest)';
pmat = xx-kron(xest,ones(1,nn));
P    = 1/(n+kap)*(.5*pmat*pmat');

% Get mean observation quantities
% vector star tracker 1
for j=1:nn
    T_s1e_s01    = eye(3)-cross(xx(13:15,j));
    T_s1e_b      = T_s1e_s01*T_s01_b;
    yez_st1(:,j) = mea_star2(qq(:,j),T_s1e_b,pr_st1,zeros(3,1));
end
T_s1e_s01 = eye(3)-cross(xx0(13:15,:));
T_s1e_b   = T_s1e_s01*T_s01_b;
ye0_st1   = mea_star2(qq0,T_s1e_b,pr_st1,zeros(3,1));
% vector payload
for j=1:nn
    T_pe_p0    = eye(3)-cross(xx(16:18,j));
    T_pe_b     = T_pe_p0*T_p0_b;    
    yez_p(:,j) = mea_payload(qq(:,j),T_pe_b,pr,zeros(3,1));
end
T_pe_p0 = eye(3)-cross(xx0(16:18,:));
T_pe_b  = T_pe_p0*T_p0_b;
ye0_p   = mea_payload(qq0,T_pe_b,pr,zeros(3,1));
% Mean observation
yez = [yez_st1;yez_p];
ye0 = [ye0_st1;ye0_p];
ye  = 1/(n+kap)*(kap*ye0+.5*sum(yez,2));

% Compute covariances
pyy0   = kap*(ye0-ye)*(ye0-ye)';
pyymat = yez-kron(ye,ones(1,nn));
pyy    = 1/(n+kap)*(.5*pyymat*pyymat');
pxy0   = kap*(xx0-xest)*(ye0-ye)';
pxy    = 1/(n+kap)*(.5*pmat*pyymat');
pvv    = pyy+R;

% Update state vector and covariance
gain = real(pxy*inv(pvv));
P    = P-gain*pvv*gain';
xest = xest+gain*(y(4:9,i+1)-ye);

% Update quaternion and bias and reset atitude error to zero
qerr      = dp2dq(xest(1:3),fac,par);
qe        = q_mult(qerr,qq0); %qe=qe/norm(qe);
xest(1:3) = zeros(3,1);
be=xest(4:6); ge=xest(7:12); s1e=xest(13:15); pe=xest(16:18); 
Me= m_fun(ge,we);
Og        = omega_g(wgm-be);
%we        = T_b_g0*[(wgm-be)+Og*ge];
we        = T_b_g0*[eye(3)+Me]*[wgm-be];
qe_his(:,i+1)=qe; we_his(:,i+1)=we; be_his(:,i+1)=be; ge_his(:,i+1)=ge; 
s1e_his(:,i+1)=s1e; pe_his(:,i+1)=pe; P_cov(:,i+1)=diag(P);

waitbar(i/(tf/dts),wait)
end
close(wait)

% these are the parameters to be sent back
qe=qe_his; we=we_his; be=be_his; ge=ge_his; s1e=s1e_his; pe=pe_his;
