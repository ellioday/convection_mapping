;----------------------------------------------------------------------
; FUNCTION  NORM_THETA
;
; PURPOSE:
;   This function returns the adjusted values of the 
;   angle theta, such that the full range of theta values run
;   from 0 to Theta_Limit, where Theta_Limit is either pi or pi/2
;
; $Log: norm_theta.pro,v $
; Revision 1.1  1997/08/05 14:50:53  baker
; Initial revision
;
;
function norm_theta,theta,thetamax,theta_limit=theta_limit,alpha=alpha
;
if (NOT keyword_set(theta_limit)) then theta_limit = !pi
theta_limit=(NOT keyword_set(theta_limit))?!pi :theta_limit
;
alpha = theta_limit/thetamax
;
theta_prime = alpha * theta
return,theta_prime
end

;----------------------------------------------------------------------

; This function creates a Lat/Long grid and then flattens it
; to return a list of lat/long position.
;
; The grid spacing is 2 degrees in latitude and 5 degrees in longitude
;
; $Log: create_grid.pro,v $
; Revision 1.2  1997/11/21 16:10:59  baker
; added provision for creating a complete 0-360 degree grid.
;
;
function create_grid,latmin,dlat,dlon,complete=complete
if (n_params() LE 1) then begin
  d_lat = 2
  d_lon = 5
endif else begin
  d_lat = dlat
  d_lon = dlon
endelse

nlats = fix((90 - latmin)/d_lat)
nlongs = fix(360/d_lon)
nlongs = (keyword_set(complete))?nlongs + 1:nlongs
pos = fltarr(2, nlats, nlongs)

for i = 0,nlats-1 do begin
  lat = i*d_lat + latmin
  for j = 0,nlongs-1 do begin
    lon = j*d_lon
    pos(0,i,j)=float(lat)
    pos(1,i,j)=float(lon)
  endfor
endfor
return,reform(pos,2,nlats*nlongs)
end

;----------------------------------------------------------------------

;+
;  FUNCTION EVAL_LEGENDRE
;
;  Purpose:  evalaute all the Associated Legendre Polynomials
;            from L=0 to L=Lmax, at a set of points
;
;  Calling Sequence:
;    PLM = eval_legendre(Lmax, pos)
;
;          where Lmax give the maximum order of the Legendre
;          polynomials, and pos gives the positions where the
;          the polynomials are to be evaluated.
; 
;          pos = pos(2,N) where the first index indicates
;                the theta,phi position and the second index
;                lists the points where we have data.
;
;
;-
; $Log: eval_legendre.pro,v $
; Revision 1.1  1997/08/05 14:42:55  baker
; Initial revision
;
;
function dbang,Mmax
result = replicate(1.0d0,Mmax+1)
for i=1,MMax do $
  for j = i, Mmax do result(j) = result(j)*(2*i - 1.0)
return, result
end

;----------------------------------------------------------------------

;+
; PURPOSE
;  This routine converts a Legendre polynoiomial index pair (l,m)
;  into a single index (k).
;
; FUNCTION:  INDEX_LEGENDRE
;
; Calling Sequence:
;   k = index_legendre(l,m,[/SH],[/SNGL],[/DBL])
;
; The keywords SH, SNGL, and DBL have the following
; meanings:
;   /SH:  We are doing Spherical harmonics where m runs from -l to + l
;
;   /SINGLE:  We are doing Associated Legendre Polynomials with m=0,l
;
;   /DOUBLE:  We are doing Associated Legendre Polynomials
;             but for each polynomial we have two coefficients 
;             one for cos(phi) and one for sin(phi), as before, m
;             runs from 0 to l.  Basically, /DOUBLE means we
;             are doing spherical harmonics for a real valued 
;             function using sin(phi) and cos(phi) rather than
;             exp(i*phi).
;-
; $Log: index_legendre.pro,v $
; Revision 1.1  1997/08/05 14:47:30  baker
; Initial revision
;
;
function index_legendre,l,m,sh=sh,sngl=sngl,dbl=dbl
if (m GT l) then return,-1
;
if (keyword_set(SH)) then begin
  return, fix(m+l*(l+1))
endif $
else if(keyword_set(SNGL)) then begin
  return, fix(m + l*(l+1)/2)
endif $
else if(keyword_set(DBL)) then begin
  if (l eq 0) then return,0 $
  else if (m eq 0) then return, L*L $
  else return, L*L + 2*m - 1
endif $
else begin
  print,'INDEX_LEGENDRE:  you must specify one and only one of'
  print,'the keywords /SH (spherical harmonics),'
  print,'/SNGL (single Legendre polynomial),'
  print,'/DBL (cos(phi),sin(phi) pairs with Legendre Polynomials'
endelse
end

;------------------------------------------------------------------

;+
; FUNCTION EVAL_POTENTIAL
;
; PURPOSE:  evaluate the electric potential on a set of
;           points, given the coefficients of the spherical
;           harmonic expansion.
;
; Calling Sequence:
;
;   pot = eval_potential,a,plm,phi
;
;     where 'a' is the set of coefficients given as a vector
;               indexed by k = index_legendre(l,m,/dbl)
;
;     plm is an array (N,Lmax,Lmax) where N = number of points
;         where the potential is to be evaluated, and the
;         other two indices give the associated Legendre polynomial
;         to be evaluated at each position.
;
;     phi is the azimuthal coordinate at each evaluation point.
;
;-
function eval_potential,a,plm,phi
ss=size(plm)
lmax=ss(2)-1
v = replicate(0.,n_elements(phi))
for m=0,lmax do begin
  for L=m,Lmax do begin
  k = index_legendre(L,m,/dbl)
  v = (m eq 0) ? v + a(k)*plm(*,L,0) : $
    v + a(k)*cos(m*phi)*plm(*,l,m) + a(k+1)*sin(m*phi)*plm(*,l,m)
  endfor
endfor
;
return,v
end

;----------------------------------------------------------------------

function eval_legendre,Lmax,x_in,pm=pm,pp=pp
x = double(reform(x_in,n_elements(x_in)))
N = n_elements(x)
xx = (1-x^2)
xx_L = xx # replicate(1.0d0,Lmax+1)
mover2 = dindgen(Lmax+1)/2.0
mover2 = replicate(1.0d0,N) # mover2
xx_Mover2 = xx_L^mover2
; xx_Mover2 is the matrix of all the (1-x^2) values raised to
; all the powers of m/2 with m running from 0 to Lmax.
;
two_m_dbang = dbang(Lmax)
two_m_dbang = replicate(1.0d0,N) # two_m_dbang
pwrm = replicate(1.0d0,N) # ((-1)^indgen(Lmax+1))
pmm = xx_Mover2*pwrm*two_m_dbang
; 
;  we have now computed the value of Pmm at each point for
;  each value of m from 0 to Lmax
;
;
; We have now computed P(m+1,m) at each point and for each
; value of m from 0 to Lmax.
;
;
;  p(m+1,m) = x(2m+1)P(m,m)
;

pmmp1 = (x # replicate(1.0d0,Lmax+1))*(replicate(1.0d0,N) # $
                                     (dindgen(Lmax+1)*2.0+1.0)) * pmm
;
; OK, now we have pmm and pmmp1 at every point.
; Next we have to compute the rest of the plm values from
; the recursion relation:
; (l-m)P(l,m) = x(2l-1)p(l-1,m) - (l+m-1)P(l-2,m)
;
plm = replicate(0.0d0,N,lmax+1,lmax+1)
for l = 0,lmax do plm(*,l,l) = pmm(*,l)
for l = 0,lmax-1 do plm(*,l+1,l) = pmmp1(*,l)
for l = 0, lmax -2 do $
   for k=l+2,Lmax do plm(*,k,l) = $
               (1.0d0/(k-l))*((x*plm(*,k-1,l)*(2*k -1))-(k+l-1)*plm(*,k-2,l))


if (keyword_set(pm)) then pm = pmm
if (keyword_set(pp)) then pp = pmmp1
return,plm
end

;----------------------------------------------------------------------

;+
; PURPOSE
;  This set of routines is used to evaluate the Electric field
;  given a set of coefficients defining the potential.
;
;-
function eval_etheta_coeffs,Pcoeffs,theta,Lmax=Lmax,Latmin=Latmin
;
theta_max = (90.0d0-Latmin)*!dtor
alpha = 1.0d0
theta_prime = norm_theta(theta,theta_max,alpha=alpha)
;
Re = 6357.0d0*1000.
n = n_elements(theta)
;
kmax = index_legendre(Lmax,Lmax,/dbl)
ecoeffs = dblarr(kmax+2,n)
q = where(theta_prime ne 0.0)
;
for m=0,Lmax do begin
  for L=m,Lmax do begin
    k3 = index_legendre(L,m,/dbl)
    k4 = index_legendre(L,m,/dbl)
    if (k3 ge 0) then begin
      ecoeffs(k4,q) = ecoeffs(k4,q) - $
        (Pcoeffs(k3)*alpha*L*cos(theta_prime(q))/sin(theta_prime(q)))/Re
    endif
;
    k1 = (L LT Lmax) ? index_legendre(L+1,m,/dbl) : -1
    k2 = index_legendre(L,m,/dbl)
    if (k1 ge 0) then begin
      ecoeffs(k2,q) = ecoeffs(k2,q) $
        + (Pcoeffs(k1)*alpha*(L+1+m)/sin(theta_prime(q)))/Re
    endif
;
    if (m gt 0) then begin
      if (k3 ge 0) then k3 = k3 + 1
      k4 = k4 + 1
      if (k1 ge 0) then k1 = k1 + 1
      k2 = k2 + 1
      if (k3 ge 0) then begin
        ecoeffs(k4,q) = ecoeffs(k4,q) - $
          (Pcoeffs(k3)*alpha*L*cos(theta_prime(q))/sin(theta_prime(q)))/Re
      endif
      if (k1 ge 0) then begin
        ecoeffs(k2,q) = ecoeffs(k2,q) $
          + (Pcoeffs(k1)*alpha*(L+1+m)/sin(theta_prime(q)))/Re
      endif
    endif
  endfor
endfor
return,ecoeffs
end
;
;--------------------------------------------------------------
function eval_ephi_coeffs,Pcoeffs,theta,Lmax=Lmax,Latmin=Latmin
;
Re = 6357.0d0*1000.
n = n_elements(theta)
;
kmax = index_legendre(Lmax,Lmax,/dbl)
ecoeffs = dblarr(kmax+2,n)
q = where(theta ne 0)
for m=1,Lmax do begin
  for L=m,Lmax do begin
    k3 = index_legendre(L,m,/dbl)
    k4 = index_legendre(L,m,/dbl)
    if (k3 ge 0) then begin
      ecoeffs(k4,q) = ecoeffs(k4,q) - Pcoeffs(k3+1)*m/sin(theta(q))/Re
      ecoeffs(k4+1,q) = ecoeffs(k4+1,q) + Pcoeffs(k3)*m/sin(theta(q))/Re
    endif
  endfor
endfor
return,ecoeffs
end
;
;--------------------------------------------------------------------
function eval_component,ecoeffs,plm,pos,Lmax=Lmax
theta = (90.0d0 - pos(0,*))*!dtor
phi = pos(1,*)*!dtor
n = n_elements(theta)
ecomp = dblarr(n)
;
for m = 0,Lmax do begin
  for L=m,Lmax do begin
    k = index_legendre(L,m,/dbl)
    if (m eq 0) then ecomp = ecomp + ecoeffs(k,*)*plm(*,l,m) $
    else ecomp = ecomp + ecoeffs(k,*)*plm(*,l,m)*cos(m*phi) + $
                           ecoeffs(k+1,*)*plm(*,l,m)*sin(m*phi)
  endfor
endfor
return,ecomp
end
;
;-----------------------------------------------------------------------
function eval_efield,pcoeffs,plm,pos,Lmax=Lmax,Latmin=Latmin
theta = pos(0,*)
theta = (90.0d0 - theta)*!dtor
etc = eval_etheta_coeffs(pcoeffs,theta,Lmax=Lmax,Latmin=Latmin)
epc = eval_ephi_coeffs(pcoeffs,theta,Lmax=Lmax,latmin=latmin)
;
etheta = eval_component(etc,plm,pos,Lmax=Lmax)
ephi = eval_component(epc,plm,pos,Lmax=Lmax)
n = n_elements(theta)
e_field = fltarr(2,n)
e_field(0,*)=etheta
e_field(1,*)=ephi
return,e_field
end
;
;-----------------------------------------------------------------------
function calc_efield,pos,solution,latmin,lon_shft,lat_shft,order

theta = (90.0 - pos(0,*))*!dtor
thetamax = (90.0-latmin)*!dtor
theta_prime = norm_theta(theta,thetamax)
x = cos(theta_prime)
plm = eval_legendre(order, x)
Lmax=order
pcoeffs=solution(2,*)

etc = eval_etheta_coeffs(pcoeffs,theta,Lmax=Lmax,Latmin=Latmin)
epc = eval_ephi_coeffs(pcoeffs,theta,Lmax=Lmax,latmin=latmin)
;
etheta = eval_component(etc,plm,pos,Lmax=Lmax)
ephi = eval_component(epc,plm,pos,Lmax=Lmax)
n = n_elements(theta)
e_field = fltarr(2,n)
e_field(0,*)=etheta
e_field(1,*)=ephi
return,e_field
end
;
;---------------------------------------------------------------------
function eval_vel,pcoeffs,plm,pos,Lmax=Lmax,Latmin=Latmin
;
e_field = eval_efield(pcoeffs,plm,pos,Lmax=Lmax,Latmin=Latmin)

Re = 6375.0d*1000.
Altitude = 300.0*1000.0
bpolar = -.62e-4
phi=90-pos(0,*)
bmag = bpolar*(1.0-3.0*Altitude/Re)*sqrt(3.0*cos(phi*!dtor )^2+1.0)/2.0
vel = e_field
vel(0,*)= e_field(1,*)/bmag
vel(1,*) = - e_field(0,*)/bmag
return,vel
end

;-------------------------------------------------------------------

;+
;  This routine calculates the 2-d velocities at
;  a set of points.
;
function calc_vels,pos,solution,latmin,order
;
theta = (90.0 - pos(0,*))*!dtor
thetamax = (90.0-latmin)*!dtor
theta_prime = norm_theta(theta,thetamax)
x = cos(theta_prime)
plm = eval_legendre(order, x)
vvec = eval_vel(solution(2,*),plm,pos,Lmax=order,latmin=latmin)
return,vvec
end

;-------------------------------------------------------------

function find_gradV,pos,solution,latmin,lon_shft,lat_shft,order

; First shift coordinates into 'model' reference (pole shifted 4 deg nightwards)
posx=pos ;lat,lon
if (lat_shft ne 0) then begin
  npnts=N_ELEMENTS(posx)/2
  kaz=FLTARR(npnts)
  crd_shft,lon_shft,lat_shft,npnts,posx,kaz
endif

; Calculate vectors
vvec = calc_vels(posx,solution,latmin,order)
vmag = reform(sqrt(vvec(0,*)^2 + vvec(1,*)^2))

q = where (vmag ne 0, qc)
if (qc eq 0) then begin
  print,'%%plot_gradV - all vectors have 0 length'
  return,0
endif
vaz    =  fltarr(n_elements(vmag))
vaz(q) =  atan(vvec(1,q),-vvec(0,q))*!radeg

; Now shift back into 'real word'
if (lat_shft ne 0) then begin
  xat_shft = -lat_shft
  npnts    =  n_elements(vmag)
  crd_shft,lon_shft,xat_shft,npnts,posx,vaz
endif

return,{mag:vmag,az:vaz}
end
