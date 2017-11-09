type PMLCoef
     vz::SparseMatrixCSC{Float64, Int64}
     vx::SparseMatrixCSC{Float64, Int64}
     pz::SparseMatrixCSC{Float64, Int64}
     px::SparseMatrixCSC{Float64, Int64}
end

function InitCpml(nz::Int64, nx::Int64, ext::Int64, dz::Float64, dx::Float64, vmax::Float64, iflag::Int64)
      if ext == 10
         R = 0.01
      elseif ext == 20
         R = 0.001
      elseif ext == 30
         R = 0.0001
      else
         error("unsupported damping layers")
      end
      ll = 2;
      if iflag == 1
         Nz = nz + 2*ext
         Nx = nx + 2*ext
         vz = spzeros(Nz, Nx)
         vx = spzeros(Nz, Nx)
         pz = spzeros(Nz, Nx)
         px = spzeros(Nz, Nx)
         #  upper and left
         for i = 1 : ext
             vz[i,:] = -1.5 * vmax/(ext*dz) * log(R) *  ((ext+1-i)/ext)^ll
             vx[:,i] = -1.5 * vmax/(ext*dx) * log(R) *  ((ext+1-i)/ext)^ll
             pz[i,:] = -1.5 * vmax/(ext*dz) * log(R) *  ((ext+1-i)/ext)^ll
             px[:,i] = -1.5 * vmax/(ext*dx) * log(R) *  ((ext+1-i)/ext)^ll
         end
         #  lower
         for iz = nz+ext+1 : Nz
             vz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz-ext)/ext)^ll
             pz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz-ext)/ext)^ll
         end
         #  right
         for ix = nx+ext+1 : Nx
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
         end
      elseif iflag == 2
         Nz = nz +   ext
         Nx = nx + 2*ext
         vz = spzeros(Nz, Nx)
         vx = spzeros(Nz, Nx)
         pz = spzeros(Nz, Nx)
         px = spzeros(Nz, Nx)
         #  lower
         for iz = nz+1 : Nz
             vz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz)/ext)^ll
             pz[iz,:] = -1.5 * vmax/(ext*dz) * log(R) * ((iz-nz)/ext)^ll
         end
         #  left
         for ix = 1 : ext
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ext+1-ix)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ext+1-ix)/ext)^ll
         end
         #  right
         for ix = nx+ext+1 : Nx
             vx[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
             px[:,ix] = -1.5 * vmax/(ext*dx) * log(R) * ((ix-nx-ext)/ext)^ll
         end
      end
      Cpml = PMLCoef(vz, vx, pz, px)
      return Cpml
end

# =============================partial derivative matrix========================
type FidMtx
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dz :: Float64
     dx :: Float64
     dt :: Float64
     MvzBvz :: SparseMatrixCSC{Float64,Int64}
     MvzBp  :: SparseMatrixCSC{Float64,Int64}
     MvxBvx :: SparseMatrixCSC{Float64,Int64}
     MvxBp  :: SparseMatrixCSC{Float64,Int64}
     MpzBpz :: SparseMatrixCSC{Float64,Int64}
     MpzBvz :: SparseMatrixCSC{Float64,Int64}
     MpxBpx :: SparseMatrixCSC{Float64,Int64}
     MpxBvx :: SparseMatrixCSC{Float64,Int64}
end

function modExpand(par::Array{Float64,2}, ext::Int64, iflag::Int64)
    temp1 = repmat(par[:,  1], 1, ext)
    temp2 = repmat(par[:,end], 1, ext)
    par   = hcat(temp1, par, temp2)
    if iflag == 1
       temp1 = repmat(par[1,  :]', ext, 1)
       temp2 = repmat(par[end,:]', ext, 1)
       par   = vcat(temp1, par, temp2)
    elseif iflag == 2
       temp2 = repmat(par[end,:]', ext, 1)
       par   = vcat(par, temp2)
    end
    return par
end

function Mvx(vx_pml::SparseMatrixCSC{Float64,Int64}, dx::Float64, dt::Float64, ext::Int64)
    (m,n) = size(vx_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vx_pml[iz,ix]/2
            b = 1/dt - vx_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for ix = 2: n-2
        tmp[ix,ix  ] = -c1; tmp[ix, ix+1] = c1;
        tmp[ix,ix-1] = -c2; tmp[ix, ix+2] = c2;
    end
    tmp[n-1,n-1] = -c1; tmp[n-1, n  ] = c1; tmp[n-1,n-2] = -c2;
    tmp[n  ,n  ] = -c1; tmp[n  , n-1] =-c2;
    MvxBvx = spdiagm(C3)
    MvxBp  = spdiagm(denum) * kron(tmp, speye(m))
    return MvxBvx, MvxBp
end

function Mvz(vz_pml::SparseMatrixCSC{Float64,Int64}, dz::Float64, dt::Float64, ext::Int64)
    (m,n) = size(vz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + vz_pml[iz,ix]/2
            b = 1/dt - vz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = 1 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] = -c1; tmp[1,2] = c1; tmp[1,3] = c2;
    for iz = 2: m-2
        tmp[iz,iz  ] = -c1; tmp[iz, iz+1] = c1;
        tmp[iz,iz-1] = -c2; tmp[iz, iz+2] = c2;
    end
    tmp[m-1,m-1] = -c1; tmp[m-1, m  ] = c1; tmp[m-1, m-2] = -c2;
    tmp[m  ,m  ] = -c1; tmp[m  , m-1] =-c2;
    MvzBvz = spdiagm(C3)
    MvzBp  = spdiagm(denum) * kron(speye(n), tmp)
    return MvzBvz, MvzBp
end

function Mpx(px_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dx::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(px_pml)
    a1 = 9/8  ; a2 = -1/24;
    c1 = a1/dx; c2 = a2/dx;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + px_pml[iz,ix]/2
            b = 1/dt - px_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(n, n)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for ix = 3: n-1
        tmp[ix,ix-1] = -c1; tmp[ix, ix  ] = c1;
        tmp[ix,ix-2] = -c2; tmp[ix, ix+1] = c2;
    end
    tmp[n,n-1] = -c1; tmp[n,n] = c1; tmp[n,n-2] = -c2;
    MpxBpx = spdiagm(C3)
    MpxBvx = spdiagm(denum) * kron(tmp, speye(m))
    return MpxBpx, MpxBvx
end

function Mpz(pz_pml::SparseMatrixCSC{Float64,Int64}, v::Array{Float64,2}, dz::Float64, dt::Float64, ext::Int64, iflag::Int64)
    v = modExpand(v, ext, iflag)
    (m, n) = size(pz_pml)
    a1 = 9/8;   a2 = -1/24;
    c1 = a1/dz; c2 = a2/dz;
    C3 = zeros(m*n)
    denum = zeros(m*n)
    for ix = 1 : n
        for iz = 1 : m
            a = 1/dt + pz_pml[iz,ix]/2
            b = 1/dt - pz_pml[iz,ix]/2
            denum[(ix-1)*m+iz] = (v[iz,ix])^2 / a
            C3[(ix-1)*m+iz]= b / a
        end
    end
    tmp = spzeros(m, m)
    tmp[1,1] =  c1; tmp[1,2] = c2;
    tmp[2,1] = -c1; tmp[2,2] = c1; tmp[2,3] = c2;
    for iz = 3: m-1
        tmp[iz,iz-1] = -c1; tmp[iz, iz  ] = c1;
        tmp[iz,iz-2] = -c2; tmp[iz, iz+1] = c2;
    end
    tmp[m,m-1] = -c1; tmp[m,m] = c1; tmp[m, m-2] = -c2
    MpzBpz = spdiagm(C3)
    MpzBvz = spdiagm(denum) * kron(speye(n), tmp)
    return MpzBpz, MpzBvz
end

function DispStable!(vmax::Float64, vmin::Float64, f0::Float64, dz::Float64, dx::Float64, dt::Float64)
    h = vmin / 5 / f0
    if dz >= h || dx >=h
       warn("large grid size may produce frequency dispersion, h should less than $h")
    end
    h = dx
    if dz < dx
       h = dz
    end
    tt = 6*h / (7*sqrt(3)*vmax)
    if dt >= tt
       warn("large time step may cause unstable, dt should less than $tt")
    end
    return nothing
end

function InitFidMtx(nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dz::Float64, dx::Float64, dt::Float64, vmax::Float64, vmin::Float64, f0::Float64, v::Array{Float64,2})
    DispStable!(vmax, vmin, f0, dz, dx, dt)
    Cpml = InitCpml(nz, nx, ext, dz, dx, vmax, iflag)
    (MvzBvz, MvzBp)  = Mvz(Cpml.vz,    dz, dt, ext       )
    (MvxBvx, MvxBp)  = Mvx(Cpml.vx,    dx, dt, ext       )
    (MpzBpz, MpzBvz) = Mpz(Cpml.pz, v, dz, dt, ext, iflag)
    (MpxBpx, MpxBvx) = Mpx(Cpml.px, v, dx, dt, ext, iflag)
    fidMtx = FidMtx(nz, nx, ext, iflag, dz, dx, dt, MvzBvz, MvzBp, MvxBvx, MvxBp, MpzBpz, MpzBvz, MpxBpx, MpxBvx)
    return fidMtx
end

# =========================forward modeling ====================================
type SnapShot
     isz:: Int64
     isx:: Int64
     nz :: Int64
     nx :: Int64
     ext:: Int64
     iflag :: Int64
     dt :: Float64
     it :: Int64
     vz :: Array{Float64, 1}
     vx :: Array{Float64, 1}
     pz :: Array{Float64, 1}
     px :: Array{Float64, 1}
end

function InitSnapShot(isz::Int64, isx::Int64, nz::Int64, nx::Int64, ext::Int64, iflag::Int64, dt::Float64, it::Int64)
    if iflag == 1
       Nz = nz + 2*ext
       Nx = nx + 2*ext
    elseif iflag == 2
       Nz = nz +   ext
       Nx = nx + 2*ext
    end
    spt = SnapShot(isz, isx, nz, nx, ext, iflag, dt, it, zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx), zeros(Nz*Nx))
    return spt
end

function Base.copy(snapShot1::SnapShot, snapShot2::SnapShot)
    if snapShot2.nz != snapShot1.nz || snapShot2.nx != snapShot1.nx || snapShot2.ext != snapShot1.ext || snapShot2.dt != snapShot1.dt || snapShot2.iflag != snapShot1.iflag
       error("the two snapShot are different")
    end
    snapShot1.isz = snapShot2.isz
    snapShot1.isx = snapShot2.isx
    snapShot1.it  = snapShot2.it
    snapShot1.vz[:] = snapShot2.vz[:]
    snapShot1.vx[:] = snapShot2.vx[:]
    snapShot1.pz[:] = snapShot2.pz[:]
    snapShot1.px[:] = snapShot2.px[:]
    return nothing
end

# =========================forward modeling ====================================
type Source
     iz :: Int64
     ix :: Int64
     dt :: Float64
     nt :: Int64
     p  :: Array{Float64,1}
end

function Ricker(f0::Float64, dt::Float64)
	  nw = 2.2/f0/dt
	  nw = 2*floor(Int,nw/2)+1
	  nc = floor(Int,nw/2)
	  w  = zeros(nw)
	  k  = collect(1:nw)
    k  = vec(k)
	  alpha = (nc-k+1)*f0*dt*pi
	  beta = alpha.^2
	  w = (1.-beta.*2).*exp(-beta)
	  return w
end

function InitSource(isz::Int64, isx::Int64, f0::Float64, dt::Float64)
    p = Ricker(f0, dt)
    nt= length(p)
    src = Source(isz, isx, dt, nt, p)
    return src
end

function AddSource!(spt::SnapShot, src::Source)
    dt = spt.dt
    if spt.dt != src.dt
       error("time sampling dismatch")
    end
    nz = spt.nz
    nx = spt.nx
    ext= spt.ext
    iflag = spt.iflag
    if iflag == 1
       zupper = ext
       Nz  = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz  = nz +   ext
    end
    Nx  = nx + 2*ext
    it  = spt.it
    if 0.0 <= (it-1)*dt <= (src.nt-1)*dt
       indt = it
       indz = src.iz + zupper
       indx = src.ix + ext
       pos = (indx-1) * Nz + indz
      #  spt.vz[pos] = spt.vz[pos] + src.p[indt]
       spt.pz[pos] = spt.pz[pos] + (src.p[indt])/2
       spt.px[pos] = spt.px[pos] + (src.p[indt])/2
    end
    return nothing
end

type ShotGather
     isz :: Int64
     isx :: Int64
     irz :: Array{Int64, 1}
     irx :: Array{Int64, 1}
     ot  :: Float64
     nt  :: Int64
     dt  :: Float64
     d   :: Array{Float64, 2}
end

function InitShotGather(isz::Int64, isx::Int64, pos::Array{Int64,2}, ot::Float64, nt::Int64, dt::Float64)
    nrec = size(pos, 1)
    irecz= pos[:, 1]
    irecx= pos[:, 2]
    d    = zeros(nt, nrec)
    shot = ShotGather(isz, isx, irecz, irecx, ot, nt, dt, d)
    return shot
end

function spt2shot!(shot::ShotGather, spt::SnapShot, Nz::Int64, Nx::Int64, ext::Int64, zupper::Int64)
    nrec = length(shot.irx)
    it = spt.it
    for irec = 1 : nrec
        iz = shot.irz[irec] + zupper
        ix = shot.irx[irec] + ext
        ind= (ix-1)*Nz + iz
        shot.d[it, irec] = spt.pz[ind] + spt.px[ind]
    end
    return nothing
end

function OneStepForward!(spt2::SnapShot, spt1::SnapShot, fidMtx::FidMtx)
    spt2.it  = spt1.it + 1
    spt2.vz = fidMtx.MvzBvz * spt1.vz + fidMtx.MvzBp  * (spt1.pz+spt1.px)
    spt2.vx = fidMtx.MvxBvx * spt1.vx + fidMtx.MvxBp  * (spt1.pz+spt1.px)
    spt2.pz = fidMtx.MpzBpz * spt1.pz + fidMtx.MpzBvz *  spt2.vz
    spt2.px = fidMtx.MpxBpx * spt1.px + fidMtx.MpxBvx *  spt2.vx
    return nothing
end

function MultiStepForward(pos::Array{Int64,2}, src::Source, fidMtx::FidMtx; tmax=2.0, interval=500, print_flag=false)
    nz = fidMtx.nz ; nx = fidMtx.nx;
    ext= fidMtx.ext; iflag = fidMtx.iflag;
    if iflag == 1
       zupper = ext
       Nz = nz + 2*ext
    elseif iflag == 2
       zupper = 0
       Nz = nz +   ext
    end
    Nx = nx + 2*ext
    dt = fidMtx.dt
    nt = round(Int64, tmax/dt)+1
    isz= src.iz; isx= src.ix;
    shot = InitShotGather(isz, isx, pos, 0.0, nt, dt)
    tl = 0.0; tu = tl + (src.nt-1)*dt;
    spt1 = InitSnapShot(src.iz, src.ix, nz, nx, ext, iflag, dt, 1)
    spt2 = InitSnapShot(src.iz, src.ix, nz, nx, ext, iflag, dt, 2)
    AddSource!(spt1, src)
    spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
    for it = 2: nt
        OneStepForward!(spt2, spt1, fidMtx)
        if tl <= (it-1)*dt <= tu   #Add sources
           AddSource!(spt2, src)
        end
        copy(spt1, spt2)
        spt2shot!(shot, spt1, Nz, Nx, ext, zupper)
        if print_flag
           if mod(it, interval) == 0
              println("the current step: $it")
           end
        end
    end
    return shot
end

"""
    fidMtx = AcousticSetup(v, dz, dx, dt, f0, ext=20, iflag=1)

Setup sparse matrix for computing spatial derivation, return a composite type FidMtx

# Arguments
* `v   :: Array{Float64,2}`: P-wave velocity model
* `dz  :: Float64`: vertical grid size
* `dx  :: Float64`: horizontal grid size
* `dt  :: Float64`: size of time step
* `f0  :: Float64`: dominant frequency of Ricker wavelet

# keywords Arguments
* `ext=20` : number of PML absorbing boundry layers
* `iflag=1`: four sides(iflag=1) or three sides(iflag=2) absorbing boundries

# Output
* `fidMtx :: FidMtx`: composite type of sparse matrix
"""
function AcousticSetup(v::Array{Float64,2}, dz::Float64, dx::Float64, dt::Float64, f0::Float64; ext=30, iflag=1)
    (nz, nx) = size(v)
    vmax = maximum(v); vmin = minimum(v)
    fidMtx = InitFidMtx(nz, nx, ext, iflag, dz, dx, dt, vmax, vmin, f0, v)
    return fidMtx
end

"""
    shot = SeisAcousticWave(fidMtx, pos, isz, isx, f0, dt, tmax=2.0)

finite difference modeling of acoustic wave field, generate a common shot gather

# Arguments
* `fidMtx :: FidMtx`        : composite type of sparse matrix
* `pos    :: Array{Int64,2}`: index of receivers, first column is the vertical index, second column contains horizontal index.
* `isz    :: Int64`         : vertical index of source
* `isz    :: Int64`         : horizontal index of source
* `f0     :: Float64`       : dominant frequency of Ricker wavelet
* `dt     :: Float64`       : size of time step

# keywords Arguments
* `tmax=1.0` : length of simulation

# Output
* `shot :: ShotGather`: composite type for common shot gather
"""
function SeisAcousticWave(fidMtx::FidMtx, pos::Array{Int64,2}, isz::Int64, isx::Int64, f0::Float64, dt::Float64; tmax=1.0, print_flag=false, interval=100)
    src  = InitSource(isz, isx, f0, dt)
    shot = MultiStepForward(pos, src, fidMtx, tmax=tmax, print_flag=print_flag, interval=interval)
    return shot
end
