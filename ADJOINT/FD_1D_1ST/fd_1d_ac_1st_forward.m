
%==========================================================================
% 1d fd acoustic modeling
% V. Etienne - Dec 2019
%
% This code compute forward and adjoint modelling for the
% 1st order wave equation, pressure-veloicty
%
% There are 4 distinct computations:
% #1 forward modelling with standard explicit FD scheme
% #2 forward modelling with matrix-vector products
% #3 adjoint modelling with standard explicit FD scheme
% #4 adjoint modelling with matrix-vector products
%
% the forward euqation:
% dP /dt = coef1 dvz/dz + S
% dvz/dt = coef2 dP/dz
%
% with
%
% coef1 = kappa = (rho * vp * vp) at pr location
% coef2 = 1/rho
% where rho is CONSTANT
%
% reference grid is P (z=0, t=0)
%
% discrete scheme
%
%  P^(n+1)_i - P^n_i           V^(n+0.5)_(i+0.5) - V^(n+0.5)_(i-0.5)
%  ----------------- = coef1_i ------------------------------------- + S^(n+0.5)_i
%         DT                                    DZ
%
%  V^(n+0.5)_(i+0.5) - V^(n-0.5)_(i+0.5)                 P^n_(i+1) - P^n_i
%  ------------------------------------- = coef2_(i+0.5) ------------------
%                   DT                                          DZ
%
%==========================================================================
clear all;
close all ;

% ***************************** input parameters **************************
NZ_MED    = 31 ;      % number of grid points in medium
global DZ ;
DZ        = 10 ;      % spatial smapling (m)
NPML_ZBEG = 10 ;      % number of grid point cpml z-

RCOEF     = 1.e-10 ;  % cpml reflection coefficient
NT        = 401 ;     % number of time steps
DT        = 0.001 ;   % time step (s)
RHO       = 1000 ;    % rho is constant
VP        = 4000 ;    % vp is constant
FREQ      = 20 ;      % frequency of source Ricker function
PERC      = 0.2 ;     % saturation coef. for figures

% please do not change the following parameters
NPML_ZEND = 0 ;       % number of grid point cpml z+ (NOT YET SUPPORTED)
LSTENCIL  = 2 ;       % half FD stencil length (ONLY FD O4 SUPPORTED)
% ***************************** end parameters *****************************

% compute usefull indexes
nz = NZ_MED + NPML_ZBEG + NPML_ZEND + 2*LSTENCIL ;
izBeg  = 1 ;
izBeg1 = izBeg  + LSTENCIL ;
izBeg2 = izBeg1 + NPML_ZBEG ;
izEnd2 = izBeg2 + NZ_MED - 1  ;
izEnd1 = izEnd2 + NPML_ZEND ;
izEnd  = izEnd1 + LSTENCIL ;

% allocate unknown arrays
NZ          = izEnd ;
pr          = zeros(NT,NZ) ;
vz          = zeros(NT,NZ) ;
npml_zBeg   = NPML_ZBEG + 2*LSTENCIL ;
mem_pr_zBeg = zeros(npml_zBeg,1) ;
mem_vz_zBeg = zeros(npml_zBeg,1) ;

% allocate static variables
coef1 = zeros(NZ,1) ;
coef2 = zeros(NZ,1) ;
apml_zBeg      = zeros(npml_zBeg,1) ;
bpml_zBeg      = zeros(npml_zBeg,1) ;
apml_half_zBeg = zeros(npml_zBeg,1) ;
bpml_half_zBeg = zeros(npml_zBeg,1) ;

% initialize model coef
for iz=1:NZ
    coef1(iz) = RHO * VP * VP * DT ;
    coef2(iz) = 1. / RHO * DT ;
end

% initialize PML
cpml_alpha_max = FREQ * pi ;
cpml_vmax      = VP ;
cpml_rcoef     = RCOEF ;

d0  = -(2+1.)*cpml_vmax*log10(cpml_rcoef) / (2.*npml_zBeg*DZ) ;
iz1 = izBeg1 ;
iz2 = izBeg2 ;

for iz=izBeg1:izBeg2-1
    ipml=iz;
    
    coef1pml = DT * VP * VP * RHO ;
    
    dnorm = (izBeg2 - iz)/(npml_zBeg) ;
    if (dnorm < 0)
        dnorm = 0. ;
    end
    alpha = cpml_alpha_max * (1. - dnorm) ;
    dd = d0 * dnorm * dnorm ;
    
    bpml_zBeg(ipml) = exp(-(dd+alpha)*DT) ;
    apml_zBeg(ipml) = dd * (bpml_zBeg(ipml)-1.) / (dd+alpha) ;
    
    dnorm = (izBeg2 - iz -0.5)/(npml_zBeg) ;
    if (dnorm < 0)
        dnorm = 0. ;
    end
    alpha = cpml_alpha_max * (1. - dnorm) ;
    dd = d0 * dnorm * dnorm ;
    
    bpml_half_zBeg(ipml) = exp(-(dd+alpha)*DT) ;
    apml_half_zBeg(ipml) = dd * (bpml_half_zBeg(ipml)-1.) / (dd+alpha) ;
end

% initialize source
% Ricker
source = zeros(NT,NZ) ;
da  = pi* FREQ ;
t00 = 1.5 * sqrt(6.) / da  ;
xsrc = floor(NZ/3) ;
for it=1:NT
    tt = it * DT ;
    aa = pi * FREQ * (tt - t00) ;
    a2 = aa * aa ;
    source (it,xsrc) = (1. - 2. * a2) * exp(-a2) ;
end

%==================================================================================
%
%            Phase  1: forward modelling with standard FD scheme
%
%==================================================================================

% loop on time steps
for it=2:NT
    
    % compute velocity at time = (it-1)*dt -1/2dt
    for iz=izBeg1:izEnd1-1
        vz(it,iz) = vz(it-1,iz) + coef2(iz) * D_Z(pr(it-1,:), iz+1) ;
    end
    
    % update cpml velocity (z-) at time = (it-1)*dt -1/2dt
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        d_pr_z = D_Z(pr(it-1,:), iz+1) ;
        mem_pr_zBeg(ipml) = bpml_half_zBeg(ipml) * mem_pr_zBeg(ipml) + apml_half_zBeg(ipml) * d_pr_z ;
        vz(it,iz) = vz(it,iz) + coef2(iz) * mem_pr_zBeg(ipml) ;
    end
    
    % compute pressure at time = (it-1)*dt
    for iz=izBeg1:izEnd1
        pr(it,iz) = pr(it-1,iz) + coef1(iz) * D_Z(vz(it,:), iz) ;
    end
    
    % update cpml pressure (z-) at time = (it-1)*dt
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        d_vz_z = D_Z(vz(it,:), iz) ;
        mem_vz_zBeg(ipml) = bpml_zBeg(ipml) * mem_vz_zBeg(ipml) + apml_zBeg(ipml) * d_vz_z ;
        pr(it,iz) = pr(it,iz) + coef1(iz) * mem_vz_zBeg(ipml) ;
    end
    
    % add source
    for iz=1:NZ
        pr(it,iz) = pr(it,iz) + source(it,iz) ;
    end
        
end

pr_forward_standard = pr ;
vz_forward_standard = vz ;

%==================================================================================
%
%            Phase  2: forward modelling with matrix vector products
%
%==================================================================================

% first build matrices
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% associate each unknown with an unique id
var_id = 1 ;
for iz=1:NZ
    pr(:,iz) = var_id ;
    var_id = var_id + 1 ;
end
for iz=1:NZ
    vz(:,iz) = var_id ;
    var_id = var_id + 1 ;
end
for iz=1:npml_zBeg
    mem_pr_zBeg(iz) = var_id ;
    var_id = var_id + 1 ;
end
for iz=1:npml_zBeg
    mem_vz_zBeg(iz) = var_id ;
    var_id = var_id + 1 ;
end
global nvar
nvar = var_id - 1 ;
fprintf("Total variables %d\n", nvar) ;

MAT_P     = eye(nvar) ;
MAT_V     = eye(nvar) ;
MAT_MEM_P = eye(nvar) ;
MAT_MEM_V = eye(nvar) ;

% add entries in matrices following computations done standard FD scheme
for it=2:NT
    
    % compute velocity at time = (it-1)*dt -1/2dt
    for iz=izBeg1:izEnd1-1
        MAT_V = add_der_op(MAT_V, vz(it,iz), pr(it,:), iz+1, coef2(iz)) ;
        %vz(it,iz) = vz(it,iz) + coef2(iz) * D_Z(pr(it-1,:), iz+1) ;
    end
    
    % update cpml velocity (z-) at time = (it-1)*dt -1/2dt
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        MAT_V = add_1term(MAT_V, vz(it,iz), mem_pr_zBeg(ipml), coef2(iz)) ;
        MAT_MEM_P = replace_1term(MAT_MEM_P, mem_pr_zBeg(ipml), mem_pr_zBeg(ipml), bpml_half_zBeg(ipml)) ;
        MAT_MEM_P = add_der_op(MAT_MEM_P, mem_pr_zBeg(ipml), pr(it,:), iz+1, apml_half_zBeg(ipml)) ;
        %d_pr_z = D_Z(pr(it-1,:), iz+1) ;
        %mem_pr_zBeg(ipml) = bpml_half_zBeg(ipml) * mem_pr_zBeg(ipml) + apml_half_zBeg(ipml) * d_pr_z ;
        %vz(it,iz) = vz(it,iz) + coef2(iz) * mem_pr_zBeg(ipml) ;
    end

    % compute pressure at time = (it-1)*dt
    for iz=izBeg1:izEnd1
        MAT_P = add_der_op(MAT_P, pr(it,iz), vz(it-1,:), iz, coef1(iz)) ;
        %pr(it,iz) = pr(it-1,iz) + coef1(iz) * D_Z(vz(it,:), iz) ;
    end
    
    % update cpml pressure (z-) at time = (it-1)*dt
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        MAT_P = add_1term(MAT_P, pr(it,iz), mem_vz_zBeg(ipml), coef1(iz)) ;
        MAT_MEM_V = replace_1term(MAT_MEM_V, mem_vz_zBeg(ipml), mem_vz_zBeg(ipml), bpml_zBeg(ipml)) ;
        MAT_MEM_V = add_der_op(MAT_MEM_V, mem_vz_zBeg(ipml), vz(it-1,:), iz, apml_zBeg(ipml)) ;
        %d_vz_z = D_Z(vz(it,:), iz) ;
        %mem_vz_zBeg(ipml) = bpml_zBeg(ipml) * mem_vz_zBeg(ipml) + apml_zBeg(ipml) * d_vz_z ;
        %pr(it,iz) = pr(it,iz) + coef1(iz) * mem_vz_zBeg(ipml) ;
    end
    
    % add source
    % nothing to do
        
    % one time step is enough
    break
end

% plot matrices

figure('Position',[100 100 900 700])
subplot(2,2,1)
%plot_matrix(MAT_P, "Forward operator P")
MAT_P_BIN=binary_matrix(MAT_P) ;
plot_matrix(MAT_P_BIN, "Forward operator P (binary)", 1)
subplot(2,2,2)
%plot_matrix(MAT_V, "Forward operator V")
MAT_V_BIN=binary_matrix(MAT_V) ;
plot_matrix(MAT_V_BIN, "Forward operator V (binary)", 1)
subplot(2,2,3)
%plot_matrix(MAT_MEM_P, "Forward operator MEM P")
MAT_MEM_P_BIN=binary_matrix(MAT_MEM_P) ;
plot_matrix(MAT_MEM_P_BIN, "Forward operator MEM P (binary)", 1)
subplot(2,2,4)
%plot_matrix(MAT_MEM_V, "Forward operator MEM V")
MAT_MEM_V_BIN=binary_matrix(MAT_MEM_V) ;
plot_matrix(MAT_MEM_V_BIN, "Forward operator MEM V (binary)", 1)

%plot_matrix(MAT_P_BIN+MAT_V_BIN+MAT_MEM_P_BIN+MAT_MEM_V_BIN, ...
%    "Forward operator Sum (binary)", 1)

% second, perform modeling with the matrices
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% allocate new variables
pr2          = zeros(NT,NZ) ;
vz2          = zeros(NT,NZ) ;
mem_pr_zBeg2 = zeros(npml_zBeg,1) ;
mem_vz_zBeg2 = zeros(npml_zBeg,1) ;

u_next = zeros(nvar,1) ;
u_prev = zeros(nvar,1) ;

% loop on time steps
for it=2:NT
    
    % update cpml velocity (z+) at time = (it-1)*dt -1/2dt
    u_next = MAT_MEM_P * u_prev ;
    u_prev = u_next ;
    
    % compute velocity at time = (it-1)*dt -1/2dt
    u_next = MAT_V * u_prev ;
    u_prev = u_next ;
    
    % store wavefield
    for iz=1:NZ
        vz2(it,iz) = u_next(vz(it,iz)) ;
    end 
    
    % update cpml pressure (z-) at time = (it-1)*dt
    u_next = MAT_MEM_V * u_prev ;
    u_prev = u_next ;
    
    % compute pressure at time = (it-1)*dt
    u_next = MAT_P * u_prev ;
    u_prev = u_next ;
    
    % add source
    for iz=1:NZ
        u_next(pr(it,iz)) = u_prev(pr(it,iz)) + source(it,iz) ;
    end
    u_prev = u_next ;
    
    % store wavefield
    for iz=1:NZ
        pr2(it,iz) = u_next(pr(it,iz)) ;
    end
           
end

% display components
pr_forward_matrix = pr2;
vz_forward_matrix = vz2;

% compute NRMS between phase 1 and phase 2
pr_nrms = compute_nrms(pr_forward_standard, pr_forward_matrix) 
vz_nrms = compute_nrms(vz_forward_standard, vz_forward_matrix) 

figure('Position',[100 100 900 700])
subplot(2,2,1)
plot_matrix(pr_forward_standard, "pr computed with standard FD scheme", 0.5)
subplot(2,2,2)
plot_matrix(vz_forward_standard, "vz computed with standard FD scheme", 0.5)
subplot(2,2,3)
plot_matrix(pr_forward_matrix, "pr with matrix-vector products", 0.5)
subplot(2,2,4)
plot_matrix(vz_forward_matrix, "vz with matrix-vector products", 0.5)

% test matrix combinations
A = MAT_V * MAT_MEM_P ;
B = MAT_P * MAT_MEM_V ;

figure('Position',[100 100 1500 500])
subplot(1,2,1)
plot_matrix_lines(binary_matrix(A), "A (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)
subplot(1,2,2)
plot_matrix_lines(binary_matrix(A'), "A' (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)

figure('Position',[100 100 1500 500])
subplot(1,2,1)
plot_matrix_lines(binary_matrix(B), "B (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)
subplot(1,2,2)
plot_matrix_lines(binary_matrix(B'), "B' (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)

BA = B*A;
figure('Position',[100 100 1500 500])
subplot(1,2,1)
plot_matrix_lines(binary_matrix(BA), "BA (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)
subplot(1,2,2)
plot_matrix_lines(binary_matrix((BA)'), "(BA)' (binary)", 1, pr, vz, mem_pr_zBeg, mem_vz_zBeg)

%==================================================================================
%
%            Phase  3: ADJOINT modelling with standard FD scheme
%
%==================================================================================

%==================================================================================
%
%            Phase  4: ADJOINT modelling with matrix vector products
%
%==================================================================================


%**********************************************************************************
%
%                        FUNCTIONS DEFINITIONS
%
%**********************************************************************************

% compute (spatial) FD derivative
function [der] = D_Z(U, iz)
global DZ ;
% FD O2
%der = (u(iz) - u(iz-1)) / DZ ;
% FD O4
der = (1.125 * (U(iz) - U(iz-1)) -1/24 * (U(iz+1) - U(iz-2))) / DZ ;
end

% add entries in matrix related to partial derivaties
function M = add_der_op(M, destid, source, iz, coef)
global DZ ;

% FD point #1
sourceid = source(iz - 2) ;
M(destid, sourceid) = M(destid, sourceid) + 1/24 / DZ * coef ;

% FD point #2
sourceid = source(iz - 1) ;
M(destid, sourceid) = M(destid, sourceid) - 1.125/ DZ * coef ;

% FD point #3
sourceid = source(iz) ;
M(destid, sourceid) = M(destid, sourceid) + 1.125 / DZ * coef ;

% FD point #4
sourceid = source(iz + 1) ;
M(destid, sourceid) = M(destid, sourceid) - 1/24 / DZ * coef ;
end

% add single entry in matrix
function M = add_1term(M, destid, sourceid, coef)
M(destid, sourceid) = M(destid, sourceid) + coef ;
end

% replace single entry in matrix
function M = replace_1term(M, destid, sourceid, coef)
M(destid, sourceid) = coef ;
end

% plot the values within a matrix
% beware that small values night not be visible
function [] = plot_matrix(M, name, perc)
%figure
hold on
title(name)
imagesc(M)
max_val=max(max(abs(M))) ;
caxis([-perc*max_val perc*max_val])
colorbar; colormap(gray) ;
axis tight ; axis ij ;
end

% plot the values within a matrix
% beware that small values night not be visible
function [] = plot_matrix_lines(M, name, perc, pr, vz, mem_pr_zBeg, mem_vz_zBeg)
global nvar 
hold on
title(name)
imagesc(M)
max_val=max(max(abs(M))) ;
caxis([-perc*max_val perc*max_val])
colorbar; colormap(gray) ;
axis tight ; axis ij ;

xline=[1 nvar] ;
zline=[max(max(pr)) max(max(pr))] ;
plot(xline, zline, '--y')
plot(zline, xline, '--y')
zline=[max(max(vz)) max(max(vz))] ;
plot(xline, zline, '--r')
plot(zline, xline, '--r')
zline=[max(max(mem_pr_zBeg)) max(max(mem_pr_zBeg))] ;
plot(xline, zline, '--b')
plot(zline, xline, '--b')
%zline=[max(max(mem_vz_zBeg)) max(max(mem_vz_zBeg))] ;
%plot(xline, zline, '--g')
%plot(zline, xline, '--g')
end

% convert a matrix M to sign(M)
% usefull to visualize the non-zero terms of a matrix
function M2 = binary_matrix(M)
size_mat=size(M) ;
M2=M;
for ix=1:size_mat(1)
    for iy=1:size_mat(2)
        M2(ix,iy) = sign(M(ix,iy)) ;
    end
end
end

% compute NRMS between 2 wavefield
function nrms = compute_nrms(u, v)
size_mat=size(u) ;
diff = 0 ;
ref  = 0 ;
for ix=1:size_mat(1)
    for iy=1:size_mat(2)
        diff = diff + (u(ix, iy) - v(ix, iy)) * (u(ix, iy) - v(ix, iy)) ;
        ref  = ref  + u(ix, iy) * u(ix, iy) ;
    end
end
nrms = sqrt(diff) / sqrt(ref) ;
end