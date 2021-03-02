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

;get values
latmin = 49.0
order = 8
lmax = order
lon_shft = 0
lat_shft = 0

mlats = [50.5, 50.5, 50.5, 51.5, 51.5, 51.5, 52.5, 52.5, 52.5, 52.5, 53.5, 51.5, 52.5, 53.5, 52.5, 52.5, 53.5, 52.5, 53.5, 53.5, 54.5, 52.5, 54.5, 53.5, 56.5, 54.5, 56.5, 55.5, 57.5, 54.5, 56.5, 55.5, 57.5, 55.5, 57.5, 58.5, 55.5, 54.5, 55.5, 77.5, 77.5, 53.5, 54.5, 54.5, 55.5, 54.5, 55.5, 56.5, 54.5, 56.5, 57.5, 57.5, 57.5, 55.5, 56.5, 57.5, 58.5, 58.5, 58.5, 55.5, 57.5, 57.5, 58.5, 59.5, 59.5, 57.5, 58.5, 59.5, 60.5, 58.5, 58.5, 59.5, 59.5, 60.5, 60.5, 61.5, 61.5, 56.5, 57.5, 59.5, 60.5, 61.5, 56.5, 57.5, 58.5, 59.5, 60.5, 57.5, 58.5, 60.5, 61.5, 56.5, 58.5, 59.5, 60.5, 61.5, 57.5, 59.5, 61.5, 58.5, 60.5, 61.5, 57.5, 59.5, 60.5, 62.5, 58.5, 59.5, 61.5, 62.5, 57.5, 60.5, 61.5, 62.5, 58.5, 60.5, 57.5, 61.5, 62.5, 58.5, 59.5, 60.5, 61.5, 62.5, 57.5, 60.5, 62.5, 61.5, 62.5, 62.5, 63.5, 64.5, 65.5, 64.5, 63.5, 65.5, 62.5, 62.5, 64.5, 61.5, 63.5, 63.5, 62.5, 63.5, 62.5, 72.5, 60.5, 61.5, 60.5, 61.5, 61.5, 59.5, 60.5, 61.5, 81.5, 82.5, 85.5, 86.5, 59.5, 60.5, 60.5, 61.5, 63.5, 61.5, 63.5, 62.5, 60.5, 61.5, 62.5, 63.5, 65.5, 66.5, 66.5, 66.5, 82.5, 82.5, 79.5, 81.5, 78.5, 80.5, 81.5, 82.5, 79.5, 80.5, 78.5, 79.5, 80.5, 81.5, 77.5, 80.5, 78.5, 79.5, 77.5, 78.5, 79.5, 61.5, 62.5, 63.5, 64.5, 65.5, 65.5, 63.5, 64.5, 61.5, 65.5, 62.5, 63.5, 64.5, 63.5, 62.5, 61.5, 54.5, 53.5, 54.5, 56.5, 57.5, 54.5, 53.5, 56.5, 56.5, 53.5, 55.5, 55.5, 52.5, 55.5, 52.5, 51.5, 52.5, 51.5, 50.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 63.5, 63.5, 63.5]
mlons = [348.2095947265625, 346.6375427246094, 345.06549072265625, 344.7321472167969, 347.9464416503906, 346.33929443359375, 342.7397155761719, 341.09588623046875, 339.4520568847656, 344.3835754394531, 340.6542053222656, 349.5535583496094, 346.02740478515625, 342.3364562988281, 347.6712341308594, 349.3150634765625, 347.3831787109375, 350.9588928222656, 349.0654296875, 350.7476501464844, 348.8038330078125, 352.6027526855469, 350.52630615234375, 352.4299011230469, 348.2412109375, 352.2488098144531, 350.0502624511719, 352.058837890625, 349.7409362792969, 353.9712829589844, 351.8592834472656, 353.8235168457031, 351.6062316894531, 355.5882263183594, 353.47149658203125, 353.2978820800781, 357.3529357910156, 359.1387634277344, 359.1176452636719, 320.76922607421875, 316.1538391113281, 291.8691711425781, 286.79425048828125, 285.07177734375, 285.0, 291.96173095703125, 286.76470947265625, 284.92462158203125, 293.6842041015625, 286.7336730957031, 284.4559631347656, 282.5906677246094, 280.72540283203125, 292.058837890625, 288.542724609375, 286.3212585449219, 284.3616943359375, 282.4468078613281, 280.53192138671875, 293.8235168457031, 290.05181884765625, 288.1865234375, 286.2765808105469, 284.2622985839844, 282.2950744628906, 291.9170837402344, 288.1914978027344, 286.2295227050781, 283.72882080078125, 292.0212707519531, 290.10638427734375, 290.1639404296875, 288.19671630859375, 287.7966003417969, 285.7627258300781, 285.6976623535156, 283.6046447753906, 295.7789001464844, 293.7823791503906, 292.1311340332031, 289.83050537109375, 287.79071044921875, 297.58795166015625, 295.6476745605469, 293.9361572265625, 294.0983581542969, 291.8644104003906, 297.512939453125, 295.85107421875, 293.8983154296875, 291.97674560546875, 299.39697265625, 297.7659606933594, 296.0655822753906, 295.93218994140625, 294.06976318359375, 299.37823486328125, 298.03277587890625, 296.16278076171875, 299.68084716796875, 297.9660949707031, 298.2558288574219, 301.2435302734375, 300.0, 300.0, 298.1927795410156, 301.5957336425781, 301.96722412109375, 300.3488464355469, 300.3614501953125, 303.1087951660156, 302.0339050292969, 302.4418640136719, 302.5301208496094, 303.5106506347656, 304.06781005859375, 304.9740905761719, 304.5348815917969, 304.69879150390625, 305.425537109375, 305.9016418457031, 306.1016845703125, 306.6278991699219, 306.8674621582031, 306.8393859863281, 308.1355895996094, 309.0361328125, 21.976743698120117, 22.77108383178711, 20.60240936279297, 21.24223518371582, 19.74193572998047, 20.53691291809082, 22.064516067504883, 23.478260040283203, 22.953020095825195, 24.93975830078125, 27.108434677124023, 29.032258987426758, 26.162790298461914, 27.95030975341797, 30.18633460998535, 29.277109146118164, 103.97515869140625, 105.18072509765625, 281.6666564941406, 350.8474426269531, 350.5813903808594, 352.88134765625, 352.6744079589844, 354.7674560546875, 266.557373046875, 267.4576416015625, 266.8604736328125, 275.0943298339844, 271.9148864746094, 276.4285583496094, 286.3636474609375, 296.0655822753906, 295.93218994140625, 293.8983154296875, 294.06976318359375, 291.8012390136719, 296.16278076171875, 294.0372619628906, 296.02410888671875, 297.9660949707031, 298.2558288574219, 298.1927795410156, 298.5093078613281, 68.85906219482422, 68.75, 71.25, 88.75, 340.85107421875, 348.5106506347656, 346.3636474609375, 349.8113098144531, 347.5, 350.8474426269531, 356.603759765625, 356.17022705078125, 351.81817626953125, 356.94915771484375, 352.5, 357.2727355957031, 3.0508475303649902, 3.396226406097412, 353.0769348144531, 9.152542114257812, 357.5, 2.7272727489471436, 357.69232177734375, 2.5, 8.181818008422852, 317.093017578125, 317.7108459472656, 316.39752197265625, 317.0322570800781, 315.302001953125, 317.7181091308594, 318.633544921875, 319.3548278808594, 319.18603515625, 320.13421630859375, 319.8795166015625, 320.86956787109375, 321.67742919921875, 323.1055908203125, 322.0481872558594, 321.2790832519531, 4.306220054626465, 5.887850284576416, 6.028707981109619, 8.140703201293945, 8.393782615661621, 7.751196384429932, 7.570093631744385, 9.949748992919922, 11.758793830871582, 9.252336502075195, 11.470588684082031, 13.235294342041016, 9.041095733642578, 15.0, 10.684931755065918, 10.446428298950195, 12.328766822814941, 12.053571701049805, 11.790392875671387, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 69.5, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 73.5999984741211, 65.4000015258789, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 77.69999694824219, 61.29999923706055, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 81.80000305175781, 57.20000076293945, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 85.9000015258789, 53.10000228881836, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 49.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 50.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 51.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 52.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 53.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 54.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 55.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 56.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 57.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 58.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 59.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 60.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 61.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 62.5, 63.5, 63.5, 63.5]

N = [0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]
N1 = [0.0, 0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 2.0, -2.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 5.0, -5.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 5.0, -5.0, 6.0, -6.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 5.0, -5.0, 6.0, -6.0, 7.0, -7.0, 0.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 5.0, -5.0, 6.0, -6.0, 7.0, -7.0, 8.0, -8.0]
N2 = [-1274.4440251679125, 1715.63115254453, -2037.6788086475722, -6896.116197325692, 3788.365984834491, 484.83133918825115, 176.0034319218624, 202.5913070835169, 58.31526676858237, 232.2041194717932, 1169.470920534148, 1391.5295741768746, -264.1147229979086, -347.07505136997565, -2.8034124380879626, 26.63348785318426, -38.924848605016905, -484.86331630360814, -269.1186668186956, -70.393402598624, -29.96580245143626, 43.33675540827561, 37.765956548798975, -6.388646644060415, -6.172807803822487, -308.5831500093806, -108.61747074039233, 261.56608594829527, 119.21159128484429, -61.6605935538417, 0.4045736526113573, -0.9794070177077709, -0.44157488405039436, 0.1919230635251737, 0.14333335908002282, 0.08983804181461968, -972.4714908993577, 112.94909069029556, -54.47904917267027, -46.72737628336666, 6.241382007141547, 2.1633232124309045, 11.62637390059969, -0.1530818863720634, -1.177785752299541, -0.15439049871570243, 0.2103704792996827, -0.013568263818995858, 0.019502184885273468, 261.41504546977166, 182.8459797596068, -55.611054393258456, 19.287066651156408, 10.14013732029748, -1.7024316511792306, -2.377143165381214, -0.6885879982065207, -0.19625369769465262, 0.09024159152221713, 0.05886011745167435, 0.007886448698615869, -0.0408868668890502, -0.00409078587600272, 0.0005133230575983134, 398.1415473145109, 163.4643497215226, 121.38365207464389, -0.9617128519446103, 4.58989685507772, -1.0832601569368516, 2.7429494461450235, 0.43792570637469025, -0.08176833752252667, -0.029680642790771074, -0.029512214936027903, -0.006175299643301253, 0.00012398701227707003, 0.0016958170536139992, 0.0008357687111500042, -0.00032659031303755066, 9.67693167443687e-05]
N3 = [0.0, 1686.8295253198437, 1649.2589988610277, 1397.3171382123937, 1143.4305694131504, 812.5546678560895, 864.9453663519236, 195.94836699453967, 246.3324231500826, 195.4237538875359, 199.62063463516822, 630.7063206943805, 689.7753103727093, 152.18281511480689, 155.704274147934, 87.09943891970357, 126.06190229079245, 518.4520633069928, 519.4363333621102, 199.39080469906244, 205.4277946712654, 148.23528241675868, 179.21326367640827, 221.5232442311592, 100.57892618452863, 59.870839170893596, 76.94617862038693, 278.87143701209243, 380.4272351799503, 149.39438143918676, 173.63984865090995, 30.31345012353916, 65.33982332690516, 81.81170046762972, 40.694049145764964, 259.32597548185686, 312.3682805381788, 61.20078472506768, 33.04092395836498, 62.35133321838244, 13.21008680113961, 174.20590784260006, 183.7690805509789, 54.16059143766868, 47.0655832951608, 42.892163469322156, 32.24927201716272, 180.69625615210302, 93.65681480487446, 138.99974109408205, 111.05190295633071, 90.93195098802204, 44.50299498260312, 47.271420973002364, 19.05647906099678, 113.39484263856787, 162.18598960030144, 111.70105664858129, 76.83381253829765, 113.78516321824412, 33.8304040763824, 67.62692647327526, 47.596126947586264, 62.92431365301553, 30.83485384683912, 65.35716427843252, 74.80282219952842, 20.87028892961725, 5.713676478836523, 21.71361823381488, 86.26747230826743, 118.58303514832917, 4.703984737615943, 22.473959212266866, 30.424931081756764, 36.16572000656307, 6.404803957868929, 1.8587258293744513, 1.6696424228100053, 3.028252607192842, 1.187370939166972]

pos = fltarr(2, n_elements(mlats))
for i=0,n_elements(mlats)-1 do begin
	PRINT, i
	pos[0, i] = mlats(i)
	pos[1, i] = mlons(i)
endfor

solution = fltarr(4, n_elements(N))
for i=0, n_elements(N)-1 do begin
	solution[0, i] = N(i)
	solution[1, i] = N1(i)
	solution[2, i] = N2(i)
	solution[3, i] = N3(i)
endfor

gradV = find_gradV(pos, solution, latmin, lon_shft, lat_shft, order)

;find_gradV
posx = pos
vvec1 = calc_vels(posx, solution, latmin, order)

vmag = reform(sqrt(vvec1(0,*)^2 + vvec1(1,*)^2))

q = where (vmag ne 0, qc)
if (qc eq 0) then begin
  print,'%%plot_gradV - all vectors have 0 length'
endif
vaz    =  fltarr(n_elements(vmag))
vaz(q) =  atan(vvec1(1,q),-vvec1(0,q))*!radeg

; Now shift back into 'real word'
if (lat_shft ne 0) then begin
  xat_shft = -lat_shft
  npnts    =  n_elements(vmag)
  crd_shft,lon_shft,xat_shft,npnts,posx,vaz
endif

;calc_vel
theta = (90.0 - pos(0,*))*!dtor
thetamax = (90.0-latmin)*!dtor
theta_prime = norm_theta(theta,thetamax)
x = cos(theta_prime)
plm = eval_legendre(order, x)
vvec = eval_vel(solution(2,*),plm,pos,Lmax=order,latmin=latmin)

;eval_vel
pcoeffs = solution(2,*)
Lmax = order
e_field = eval_efield(pcoeffs,plm,pos,Lmax=Lmax,Latmin=Latmin)
Re = 6375.0d*1000.
Altitude = 300.0*1000.0
bpolar = -.62e-4
phi=90-pos(0,*)
bmag = bpolar*(1.0-3.0*Altitude/Re)*sqrt(3.0*cos(phi*!dtor )^2+1.0)/2.0
vel = e_field
vel(0,*)= e_field(1,*)/bmag
vel(1,*) = - e_field(0,*)/bmag

;eval_efield
theta = pos(0,*)
theta = (90.0d0 - theta)*!dtor
etc = eval_etheta_coeffs(pcoeffs,theta,Lmax=Lmax,Latmin=Latmin)
epc = eval_ephi_coeffs(pcoeffs,theta,Lmax=Lmax,latmin=latmin)

;eval_etheta_coeffs
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

END
