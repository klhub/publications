function [qe,be,p_cov]=usque(im,bm,wgm,r,sigu,sigv,poa,pog,av,dt,x0,par,kap);

%[qe,be,p_cov]=usque(im,bm,wgm,r,sigu,sigv,poa,pog,av,dt,x0,par,kap);
%
% This program determines the attitude of a satellite using multiple sensors
% in an UnScented QUaternion Estimator with multiplicative quaternion errors
% (see Lefferts, Markley, and Shuster, JGCD Sept.-Oct. '82).
% This algorithm uses the discrete propagated error covariance and the
% Cholesky square root of the covariance matrix.
%
% The inputs are:
%     im = inertial measurements [mx(3*s)], s = number of sensors
%     bm = body measurements  [mx(3*s)]
%    wgm = gyro measurements (rad/sec) [mx3]
%      r = measurement covariance [(3*s)x(3*s)]
%   sigu = gyro bias noise standard deviation (rad/sec^1.5) [3x3]
%   sigv = gyro noise standard deviation (rad/sec^0.5) [3x3]
%    poa = initial error covariance of attitude (rad^2) [3x3]
%    pog = initial error covariance of gyro bias (rad^2/sec^2) [3x3]
%     av = 1 for sensor available, 0 for not  [mxs]
%     dt = sampling interval (sec)
%     x0 = initial estimates of quaternions and biases ([q0 b0]), [1x7]
%    par = 0 for 2 x Gibbs vector attitude error representation
%        = 1 for 4 x MRP attitude error representation
%          Try other values, too
%    kap = weight for central point (kap = -3 recommended)
%
% The outputs are:
%     qe = estimated quaternions [mx4]
%     be = estimated gyro biases (rad/sec) [mx3]
%  p_cov = diagonal covariances [mx6]

% John L Crassidis & F L Markley 3/1/02

% constants and conversions
x0=x0(:)';
[m,sen]=size(bm); sen=sen/3;
i500=0;
n=3;           % number of non-attitude components of state vector
ntot=n+3;      % total number of components of state vector
npts=2*ntot;   % number of sigma points
n1=ntot+1; n2=ntot+2; n3=ntot+3;

% Multiplication factor for attitude error representation
fac=2*(par+1);

% pre-allocate space
qe=zeros(m,4);
be=zeros(m,3);
we=zeros(m,3);
p_cov=zeros(m,ntot);
qq=zeros(4,npts);
dq=zeros(4,npts);

% initial bias, quaternion estimate, error vector
be(1,:)=x0(1,5:7);
qe(1,:)=x0(1,1:4);
xest=[zeros(3,1);be(1,:)'];

% initial state and discrete process noise covariances
p=[poa zeros(3);zeros(3) pog];
p_cov(1,:)=diag(p)';
qb=sigu.^2;
qc=0.5*dt*[sigv.^2-qb*dt^2/6 zeros(3,3); zeros(3,3) qb];

% main loop
wait = waitbar(0,'USQUE running ...');
for i=1:m-1,

% display when every 500th point is reached
if (i500==500),
    disp(sprintf('      Filter has reached point %5i',i-1))
    i500=0;
end
i500=i500+1;

% Matrix square root
psquare=chol((ntot+kap)*(p+qc))';
sigv=real([psquare -psquare]);
xx0=xest;
xx=sigv+kron(xest,ones(1,npts));

% Compute error quaternion and propagate it forward
we=wgm(i,:)'-xest(4:6);
w=norm(we);
co=cos(0.5*w*dt);
nsi=we*sin(0.5*w*dt)/w;
om=[co*eye(3)-crossm(nsi) nsi;-nsi' co];
qq0=(om*qe(i,:)');

% Compute other quaternions and propagate them forward
for j = 1:npts,
    if j==1|j==2|j==3|j==n1|j==n2|j==n3
       xxsq=xx(1:3,j)'*xx(1:3,j);
       qerr4=(fac*sqrt(fac^2+(1-par^2)*xxsq)-par*xxsq)/(fac^2+xxsq);
       qerr=[(par+qerr4)/fac*xx(1:3,j); qerr4];
       qmat=[qerr(4)*eye(3)-crossm(qerr(1:3)) qerr(1:3);-qerr(1:3)' qerr(4)];
       qinit=qmat*qe(i,:)';
    else
       qinit=qe(i,:)';
    end
    we=wgm(i,:)'-xx(4:6,j);
    w=norm(we);
    co=cos(0.5*w*dt);
    nsi=we*sin(0.5*w*dt)/w;
    om=[co*eye(3)-crossm(nsi) nsi;-nsi' co];
    qq(:,j)=om*qinit;

% Get attitude error vector from new quaternions and qq0
    qmatr=[qq(4,j)*eye(3)-crossm(qq(1:3,j)) qq(1:3,j);-qq(1:3,j)' qq(4,j)];
    dq(:,j)=qmatr*[-qq0(1:3);qq0(4)];
    xx(1:3,j)=fac*dq(1:3,j)/(par+dq(4,j));
end

% Mean estimate -- note that xx0(1:3) is still zero from the reset.
xest=1/(ntot+kap)*(kap*xx0'+0.5*sum(xx,2)')';

% Mean covariance
pp0=kap*(xx0-xest)*(xx0-xest)';
pmat=xx-kron(xest,ones(1,npts));
p=1/(ntot+kap)*(pp0+0.5*pmat*pmat')+qc;

% Get mean observation quantities
yez=[];pbee=[];ye0=[];r1=[];ym=[];
% See which sensors are available
[i1,j1]=find(av(i+1,:)==1);
for nn=1:length(j1),
    for j = 1:npts,
       a=attm(qq(:,j));
       pbe=a*im(i+1,j1(nn)*3-2:j1(nn)*3)';
       pbee = [pbee pbe];
    end
    yez=[yez;pbee];
    pbee=[];
    a=attm(qq0);
    pbe0=a*im(i+1,j1(nn)*3-2:j1(nn)*3)';
    ye0=[ye0;pbe0];
    r1(nn*3-2:nn*3,nn*3-2:nn*3)=r(j1(nn)*3-2:j1(nn)*3,j1(nn)*3-2:j1(nn)*3);
    ym(nn*3-2:nn*3,1)=bm(i+1,j1(nn)*3-2:j1(nn)*3)';
end

% Mean observation
ye=1/(ntot+kap)*(kap*ye0+0.5*sum(yez,2));

% Compute covariances
pyy0=kap*(ye0-ye)*(ye0-ye)';
pyymat=yez-kron(ye,ones(1,npts));
pyy=1/(ntot+kap)*(pyy0+0.5*pyymat*pyymat');

pxy0=kap*(xx0-xest)*(ye0-ye)';
pxy=1/(ntot+kap)*(pxy0+0.5*pmat*pyymat');

pvv=pyy+r1;

% Update state vector and covariance
gain=real(pxy*inv(pvv));
p=p-gain*pvv*gain';
xest=xest+gain*(ym-ye);

% Update quaternion and bias and reset atitude error to zero
xestsq=xest(1:3)'*xest(1:3);
qerr4=(fac*sqrt(fac^2+(1-par^2)*xestsq)-par*xestsq)/(fac^2+xestsq);
qerr=[(par+qerr4)/fac*xest(1:3); qerr4];
qmat=[qerr(4)*eye(3)-crossm(qerr(1:3)) qerr(1:3);-qerr(1:3)' qerr(4)];
qe(i+1,:)=(qmat*qq0)';
be(i+1,:)=xest(4:6)';
p_cov(i+1,:)=diag(p)';
xest(1:3)=zeros(3,1);

waitbar(i/m,wait)
end
close(wait)