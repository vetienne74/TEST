
%==========================================================================
% 1d fd acoustic modeling

% 1st order wave equation, pressure-veloicty

% dP /dt = coef1 dvz/dz + S
% dvz/dt = coef2 dP/dz

% with

% coef1 = kappa = (rho * vp * vp) at pr location
% coef2 = 1/rho
% where rho is CONSTANT

% reference grid is P (z=0, t=0)

% discrete scheme

%  P^(n+1)_i - P^n_i           V^(n+0.5)_(i+0.5) - V^(n+0.5)_(i-0.5)
%  ----------------- = coef1_i ------------------------------------- + S^(n+0.5)_i
%         DT                                    DZ

%  V^(n+0.5)_(i+0.5) - V^(n-0.5)_(i+0.5)                 P^n_(i+1) - P^n_i
%  ------------------------------------- = coef2_(i+0.5) ------------------
%                   DT                                          DZ

%==========================================================================
clear all; close all ;

% ***************************** input parameters **************************
BUILD_MATRIX = 1 ;

NZ_MED    = 31 ;
DZ        = 10 ;
NPML_ZBEG = 10 ;
NPML_ZEND = 0 ;
RCOEF     = 1.e-10 ;
NT        = 401 ;
DT        = 0.001 ;
RHO       = 1.0 ;
VP        = 4000 ;
FREQ      = 20 ;
LSTENCIL  = 2 ;
PERC      = 0.2 ;
% ***************************** end parameters *****************************

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

if (BUILD_MATRIX)
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
    
    MAT_P     = zeros(nvar, nvar) ;
    MAT_V     = zeros(nvar, nvar) ;
    MAT_MEM_P = zeros(nvar, nvar) ;
    MAT_MEM_V = zeros(nvar, nvar) ;
end


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

% loop on time steps
for it=2:NT
    
    % compute pressure
    for iz=izBeg1:izEnd1
        if (BUILD_MATRIX)
            MAT_P = fill_mat(MAT_P, pr(it,iz), vz(it-1,:), iz, coef1(iz)) ;
        else
            pr(it,iz) = pr(it-1,iz) + coef1(iz) * D_Z(vz(it-1,:), iz, DZ) ;
        end
    end
    
    % update cpml pressure (z-)
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        if (BUILD_MATRIX)
            MAT_P = fill_mat2(MAT_P, pr(it,iz), mem_vz_zBeg(ipml), coef1(iz)) ;
            MAT_MEM_V = fill_mat2(MAT_MEM_V, mem_vz_zBeg(ipml), mem_vz_zBeg(ipml), bpml_zBeg(ipml)) ;
            MAT_MEM_V = fill_mat(MAT_MEM_V, mem_vz_zBeg(ipml), vz(it-1,:), iz, apml_zBeg(ipml)) ;
        else
            d_vz_z = D_Z(vz(it-1,:), iz, DZ) ;
            mem_vz_zBeg(ipml) = bpml_zBeg(ipml) * mem_vz_zBeg(ipml) + apml_zBeg(ipml) * d_vz_z ;
            pr(it,iz) = pr(it,iz) + coef1(iz) * mem_vz_zBeg(ipml) ;
        end
    end
    
    % add source
    if ~(BUILD_MATRIX)
        for iz=1:NZ
            pr(it,iz) = pr(it,iz) + source(it-1,iz) ;
        end
    end
    
    % compute velocity
    for iz=izBeg1:izEnd1-1
        if (BUILD_MATRIX)
            MAT_V = fill_mat(MAT_V, vz(it,iz), pr(it,:), iz+1, coef2(iz)) ;
        else
            vz(it,iz) = vz(it-1,iz) + coef2(iz) * D_Z(pr(it,:), iz+1, DZ) ;
        end
    end
    
    % update cpml velocity (z+)
    for iz=izBeg1:izBeg2-1
        ipml = iz ;
        if (BUILD_MATRIX)
            MAT_V = fill_mat2(MAT_V, vz(it,iz), mem_pr_zBeg(ipml), coef2(iz)) ;
            MAT_MEM_P = fill_mat2(MAT_MEM_P, mem_pr_zBeg(ipml), mem_pr_zBeg(ipml), bpml_half_zBeg(ipml)) ;
            MAT_MEM_P = fill_mat(MAT_MEM_P, mem_pr_zBeg(ipml), pr(it,:), iz+1, apml_half_zBeg(ipml)) ;
        else
            d_pr_z = D_Z(pr(it,:), iz+1, DZ) ;
            mem_pr_zBeg(ipml) = bpml_half_zBeg(ipml) * mem_pr_zBeg(ipml) + apml_half_zBeg(ipml) * d_pr_z ;
            vz(it,iz) = vz(it,iz) + coef2(iz) * mem_pr_zBeg(ipml) ;
        end
    end
    
    if (BUILD_MATRIX)
        break
    end
end

if (BUILD_MATRIX)
    
    % plot matrices
    
    %plot_matrix(MAT_P, "Forward operator P")
    MAT_P_BIN=binary_matrix(MAT_P) ;
    plot_matrix(MAT_P_BIN, "Forward operator P (binary)", 1)
    %plot_matrix(MAT_V, "Forward operator V")
    MAT_V_BIN=binary_matrix(MAT_V) ;
    plot_matrix(MAT_V_BIN, "Forward operator V (binary)", 1)
    
    %plot_matrix(MAT_MEM_P, "Forward operator MEM P")
    MAT_MEM_P_BIN=binary_matrix(MAT_MEM_P) ;
    plot_matrix(MAT_MEM_P_BIN, "Forward operator MEM P (binary)", 1)
    %plot_matrix(MAT_MEM_V, "Forward operator MEM V")
    MAT_MEM_V_BIN=binary_matrix(MAT_MEM_V) ;
    plot_matrix(MAT_MEM_V_BIN, "Forward operator MEM V (binary)", 1)
    
    plot_matrix(MAT_P_BIN+MAT_V_BIN+MAT_MEM_P_BIN+MAT_MEM_V_BIN, ...
        "Forward operator Sum (binary)", 1)
    
    % perform modeling with the matrices
    %************************************
    
    % allocate new variables
    pr2          = zeros(NT,NZ) ;
    vz2          = zeros(NT,NZ) ;
    mem_pr_zBeg2 = zeros(npml_zBeg,1) ;
    mem_vz_zBeg2 = zeros(npml_zBeg,1) ;
    
    u_next = zeros(nvar,1) ;
    u_prev = zeros(nvar,1) ;
    
    % loop on time steps
    for it=2:NT
        
        % update cpml pressure (z-)
        %u_next = MAT_MEM_V * u_prev ;
        %u_prev = u_next ;
        
        % compute pressure
        u_prev(xsrc) = 1.0 ;
        u_next = MAT_P * u_prev ;                
        u_prev = u_next ;
        
        % add source
        for iz=1:NZ
            %u_prev(pr(it,iz)) = u_prev(pr(it,iz)) + source(it-1,iz) ;
        end          
        
        % store wavefield
        for iz=1:NZ
            pr2(it,iz) = u_next(pr(it,iz)) ;
        end
        
        % update cpml velocity (z+)
        %u_next = MAT_MEM_P * u_prev ;
        %u_prev = u_next ;
        
        % compute velocity
        u_next = MAT_V * u_prev ;
        u_prev = u_next ;
        
        % store wavefield
        for iz=1:NZ
            vz2(it,iz) = u_next(vz(it,iz)) ;
        end
        
    end
end

% display components
plot_matrix(pr, "pr id", 1)
plot_matrix(pr2, "pr computed with matrices", 0.1)
plot_matrix(vz2, "vz computed with matrices", 0.1)

%==========================================================================
%                        FUNCTIONS DEFINITIONS
%==========================================================================

% compute (spatial) FD derivative
function [der] = D_Z(U, iz, DZ)
% FD O2
%der = (u(iz) - u(iz-1)) / DZ ;
% FD O4
der = (1.125 * (U(iz) - U(iz-1)) -1/24 * (U(iz+1) - U(iz-2))) / DZ ;
end

% add entries in matrix related to partial derivaties
function M = fill_mat(M, destid, source, iz, coef)

% FD point #1
sourceid = source(iz - 2) ;
M(destid, sourceid) = M(destid, sourceid) + 1/24 * coef ;

% FD point #2
sourceid = source(iz - 1) ;
M(destid, sourceid) = M(destid, sourceid) - 1.125 * coef ;

% FD point #3
sourceid = source(iz) ;
M(destid, sourceid) = M(destid, sourceid) + 1.125 * coef ;

% FD point #4
sourceid = source(iz + 1) ;
M(destid, sourceid) = M(destid, sourceid) - 1/24 * coef ;
end

% add single entry in matrix
function M = fill_mat2(M, destid, sourceid, coef)
M(destid, sourceid) = M(destid, sourceid) + coef ;
end

% plot the values within a matrix
% beware that small values night not be visible
function [] = plot_matrix(M, name, perc)
figure
hold on
title(name)
imagesc(M)
max_val=max(max(abs(M))) ;
caxis([-perc*max_val perc*max_val])
colorbar; colormap(gray) ;
axis tight ; axis ij ;
end

% convert a matrix M to sign(M)
% usefull to visualize the non-zero terms of a matrix
function M2 = binary_matrix(M)
global nvar
M2=M;
for ix=1:nvar
    for iy=1:nvar
        M2(ix,iy) = sign(M(ix,iy)) ;
    end
end
end