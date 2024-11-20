% gyro misalignment 
lambda_x = 500/1e6 ;                % from ppm to fraction
delta_yz = 200/60^2*pi/180 ;        % from arc-s to rad
delta_zy = 300/60^2*pi/180 ;
delta_xz = 0/60^2*pi/180 ;
lambda_y = 500/1e6 ;
delta_zx = 400/60^2*pi/180 ;
delta_xy = 0/60^2*pi/180 ;
delta_yx = 0/60^2*pi/180 ;
lambda_z = 500/1e6 ;
mu_x     = .100/1e6 ;                % asymmetric scale factor error
mu_y     = .100/1e6 ;
mu_z     = .100/1e6 ;

%lambda_x_bar =  lambda_x ;         % these only for misalignment wrt to payload
%delta_yz_bar = -delta_yz ; 
%delta_zy_bar =  delta_zy ; 
%delta_xz_bar =  delta_xz ; 
%lambda_y_bar =  lambda_y ;
%delta_zx_bar = -delta_zx ; 
%delta_xy_bar = -delta_xy ; 
%delta_yx_bar =  delta_yx ; 
%lambda_z_bar =  lambda_z ; 
%M           = [ lambda_x_bar delta_yz_bar delta_zy_bar
%                delta_xz_bar lambda_y_bar delta_zx_bar
%                delta_xy_bar delta_yx_bar lambda_z_bar ];
                
xi_x    = -delta_zx ;               % -400 arc-sec
xi_y    =  delta_xy-delta_zy ;      % -300 arc-sec
xi_z    = -delta_yz ;               % -200 arc-sec

T_b_g0      = eye(3) ;                      % from body coordinate to gyro coordinate
%t1=45*pi/180; t2=45*pi/180; t3=45*pi/180;     %%%%%% !!!!!!! doen't work yet
%T_b_g0      = [  cos(t3)*cos(t2)  cos(t3)*sin(t2)*sin(t1)+sin(t3)*cos(t1) cos(t3)*sin(t2)*cos(t1)+sin(t3)*sin(t1)
%                -sin(t3)*cos(t2) -sin(t3)*sin(t2)*sin(t1)+cos(t3)*cos(t1) sin(t3)*sin(t2)*cos(t1)+cos(t3)*sin(t1)
%                 sin(t2)         -cos(t2)*sin(t1)                         cos(t2)*cos(t1)                         ];
T_g0_b      = T_b_g0';
%bg          = [ .2 .3 .2 ]'*pi/180/3600 ;   % initial gyro biases, rad
bg          = [ 0 0 0 ]'*pi/180/3600 ;   % since the w is already the measurement

% star tracker - quaternion
delta_s1    = [-20;-20; 20]/60^2 ;          % convert to degree only, not rad
delta_s2    = [ 20; 20; 20]/60^2 ;
delta_s1    = [0;0;0];              % !!! since w is already misaligned
delta_s2    = [0;0;0];              % !!! since w is already misaligned
q_s01_s1    = [(1/2)*delta_s1;1] ;  q_s01_s1=q_s01_s1/norm(q_s01_s1) ;  % need to normalize ?
q_s02_s2    = [(1/2)*delta_s1;1] ;  q_s02_s2=q_s02_s2/norm(q_s02_s2) ;
q_b_s01     = e2q(pi/1,pi/2,pi/3,7)' ;              % arbitrary sensor mounting, from 123 Euler sequences to quaternion
q_b_s02     = e2q(pi/4,pi/5,pi/6,7)' ;
q_b_s1      = q_mult(q_s01_s1,q_b_s01) ;
q_b_s2      = q_mult(q_s02_s2,q_b_s02) ;

% star tracker - vector
del_st1_x   = -20/60^2*pi/180 ;     % from arc-sec to radian
del_st1_y   = -20/60^2*pi/180 ;     % from arc-sec to radian
del_st1_z   =  20/60^2*pi/180 ;     % from arc-sec to radian
del_st1     = [del_st1_x;del_st1_y;del_st1_z];
del_st1     = [0;0;0];              % !!! since w is already misaligned
T_s1_s01    = eye(3)-cross(del_st1);  % from nominal sensor mounting to true sensor mounting  %%%% normalize?
T_s01_b     = q2att_mat(e2q(20*pi/180,45*pi/180,145*pi/180,7)');  % nominal sensor mounting
T_s1_b      = T_s1_s01*T_s01_b;  % from assumed (nominal) sensor mounting to actual moutning
pr_st1      = [0.64032;0.7071;0.3] ;                % arbitrary reference in unit vector

% payload
del_p_x     =  25/60^2*pi/180 ;     % from arc-sec to radian
del_p_y     =  25/60^2*pi/180 ;     % from arc-sec to radian
del_p_z     = -20/60^2*pi/180 ;     % from arc-sec to radian
del_p       = [del_p_x;del_p_y;del_p_z];
del_p       = [0;0;0];              % !!! since w is already misaligned
T_p_p0      = eye(3)-cross(del_p);  % from nominal sensor mounting to true sensor mounting  %%%% normalize?
T_p0_b      = q2att_mat(e2q(pi/4,pi/4,pi/4,7)) ;    % arbitrary sensor mounting, also from 123 Euler sequences
T_p_b       = T_p_p0*T_p0_b;
%T_p_b       = T_p0_b;
pr          = [0.7071;0.3;0.64032] ;                % arbitrary reference in unit vector

% Misalignment
g           = [ xi_x
                xi_y 
                xi_z
                lambda_x
                lambda_y
                lambda_z
                mu_x
                mu_y 
                mu_z     ] ;
s1          = delta_s1 ;
s2          = delta_s2 ;        
g           = zeros(6,1);           % !!! since w is already misaligned

% noise characteristics
sig_q       = 1e-8 ;                % quaternion noise, process noise, sigma        % from Dr Crassidis
sig_a       = 3.49065850398866E-06;             % gyro noise sigma, (rad^2/sec)^.5
sig_r       = 6.28462179216066E-09;             % gyro drift rate sigma, (rad^2/sec^3)^.5
sig_g       = [ .001/60^2*pi/180    % gyro misalignment
                .001/60^2*pi/180 
                .001/60^2*pi/180
                .001/1e6
                .001/1e6
                .001/1e6
                .001/1e6
                .001/1e6
                .001/1e6         ];
sig_s1      = [ .001/60^2*pi/180      % star tracker 1 alignment process noise
                .001/60^2*pi/180 
                .001/60^2*pi/180 
                .001/60^2*pi/180 ] ;
sig_s2      = [ 1/60^2*pi/180      % star tracker 2 alignment process noise
                1/60^2*pi/180 
                1/60^2*pi/180 
                1/60^2*pi/180 ] ;            
sig_p       = [ .001/60^2*pi/180      % Payload alignment process noise
                .001/60^2*pi/180 
                .001/60^2*pi/180 
                .001/60^2*pi/180 ] ;
sig_s1m     = [ 6/60^2*pi/180       % star tracker 1 measurement noise, quaternion measurement
                6/60^2*pi/180 
                6/60^2*pi/180 
                6/60^2*pi/180 ] ;
sig_s1pm    = [ 5/60^2*pi/180       % star tracker 1 measurement noise, vector measurement
                5/60^2*pi/180 
                5/60^2*pi/180 ] ;            
sig_s2m     = [ 1/60^2*pi/180       % star tracker 2 measurement noise
                1/60^2*pi/180 
                1/60^2*pi/180 
                1/60^2*pi/180 ] ;            
sig_pm      = [ .5/60^2*pi/180      % payload measurement noise
                .5/60^2*pi/180 
                .5/60^2*pi/180 ] ;
SIG_q       = sig_q^2*eye(3) ; 
SIG_a       = sig_a^2*eye(3) ;
SIG_r       = sig_r^2*eye(3) ;
SIG_g       = [ sig_g(1)^2 0       0       0       0       0       0       0       0
                   0    sig_g(2)^2 0       0       0       0       0       0       0
                   0       0    sig_g(3)^2 0       0       0       0       0       0
                   0       0       0    sig_g(4)^2 0       0       0       0       0
                   0       0       0       0    sig_g(5)^2 0       0       0       0
                   0       0       0       0       0    sig_g(6)^2 0       0       0
                   0       0       0       0       0       0    sig_g(7)^2 0       0
                   0       0       0       0       0       0       0    sig_g(8)^2 0
                   0       0       0       0       0       0       0       0    sig_g(9)^2 ] ;
SIG_s1      = [ sig_s1(1)^2  0        0
                    0    sig_s1(2)^2  0
                    0        0    sig_s1(3)^2 ] ;
SIG_s2      = [ sig_s2(1)^2  0        0
                    0    sig_s2(2)^2  0
                    0        0    sig_s2(3)^2 ] ;                
SIG_p       = [ sig_p(1)^2  0       0
                    0   sig_p(2)^2  0
                    0       0   sig_p(3)^2 ] ;                                
SIG_s1m     = [ sig_s1m(1)^2  0         0
                    0     sig_s1m(2)^2  0
                    0         0     sig_s1m(3)^2 ] ;
SIG_s1pm    = [ sig_s1pm(1)^2   0         0
                    0      sig_s1pm(2)^2  0
                    0           0     sig_s1pm(3)^2 ] ;                
SIG_s2m     = [ sig_s2m(1)^2  0         0
                    0     sig_s2m(2)^2  0
                    0         0     sig_s2m(3)^2 ] ;
SIG_pm      = [ sig_pm(1)^2   0        0
                   0     sig_pm(2)^2   0
                   0          0   sig_pm(2)^2 ] ;               
Qc          = [ SIG_a/1    zeros(3,3) zeros(3,9) zeros(3,3) zeros(3,3)
                zeros(3,3) SIG_r      zeros(3,9) zeros(3,3) zeros(3,3)
                zeros(9,3) zeros(9,3) SIG_g      zeros(9,3) zeros(9,3)
                zeros(3,3) zeros(3,3) zeros(3,9) SIG_s1     zeros(3,3)
                zeros(3,3) zeros(3,3) zeros(3,9) zeros(3,3) SIG_s2     ] ;
sig_g       = sig_g(1:6,1);
SIG_g       = SIG_g(1:6,1:6);       
