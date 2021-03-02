;Below are 'internal' functions and procedures which are 'invisible' to the user
;----------------------------------------------------------------------

function mappot_time_inv_f,time_index

COMMON map_pot

arr_time=[map_data(time_index).shr*100+map_data(time_index).smt, $
    	    map_data(time_index).ehr*100+map_data(time_index).emt]

return,arr_time

END

;--------------------------------------------------------------------

function mappot_time_inv_f3,time_index

COMMON map_pot

arr_time=[fill_digit(map_data(time_index).shr,2)+':'+fill_digit(map_data(time_index).smt,2)+':'+fill_digit(map_data(time_index).ssc,2), $
    	    fill_digit(map_data(time_index).ehr,2)+':'+fill_digit(map_data(time_index).emt,2)+':'+fill_digit(map_data(time_index).esc,2)]

return,arr_time

END

;--------------------------------------------------------------------
function mappot_time_f,arr_time

COMMON map_pot

IF (size(arr_time))(1) NE 7 THEN BEGIN
  time_hrs=arr_time/100
  time_mins=arr_time-time_hrs*100
  time_index=where(map_data.shr EQ time_hrs AND map_data.smt EQ time_mins)
ENDIF ELSE BEGIN
  time_hrs=fix(strmid(arr_time,0,2))
  time_mins=fix(strmid(arr_time,2,2))
  time_secs=fix(strmid(arr_time,4,2))
  time_index=where(map_data.shr EQ time_hrs AND map_data.smt EQ time_mins AND map_data.ssc EQ time_secs)
ENDELSE

return,time_index(0)

END

;--------------------------------------------------------------------
function gen_tim_tit,shr,smin,ehr,emin,ssec,esec

IF n_params() EQ 4 THEN $
tim_tit=fill_digit(shr,2)+':'+fill_digit(smin,2)+':00 - '+ $
fill_digit(ehr,2)+':'+fill_digit(emin,2)+':00 UT' $
ELSE $
tim_tit=fill_digit(shr,2)+':'+fill_digit(smin,2)+':'+fill_digit(ssec,2)+' - '+ $
fill_digit(ehr,2)+':'+fill_digit(emin,2)+':'+fill_digit(esec,2)+' UT'

return,tim_tit
END

;----------------------------------------------------------------------------------------
function gen_dat_tit,dy,mo,yr

month2=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
dat_tit=fill_digit(dy,2)+' '+month2(mo-1)+' '+fill_digit(yr,4)

return,dat_tit
END

;----------------------------------------------------------------------
function tag_pos,xl,yl

theta=atan(xl,yl)
r=sqrt(xl^2+yl^2)+3
nx=r*sin(theta)
ny=r*cos(theta)

tags=[nx,ny]
return,tags

end
;--------------------------------------------------------------
function get_model_vpc,bx,by,bz,Kp

vpc=intarr(4,2,8)
vpc(0,*,*)=[[-18,5],[-21,11],[-29,16],[-39,20],[-43,23],[-33,21],[-28,17],[-24,11]]
vpc(1,*,*)=[[-16,5],[-22,7],[-28,14],[-37,18],[-39,17],[-33,16],[-26,12],[-21,6]]
vpc(2,*,*)=[[-15,3],[-19,8],[-31,14],[-39,23],[-42,26],[-33,24],[-28,15],[-22,7]]
vpc(3,*,*)=[[-12,3],[-13,9],[-30,16],[-50,27],[-53,19],[-36,29],[-30,21],[-23,10]]

IF N_PARAMS() EQ 0 THEN return,vpc

IF N_PARAMS() EQ 1 THEN return,vpc(0,*,*)

  angle=atan(by,bz)/!dtor
  angle=angle+22.5
  angle=(angle LT 0)?angle+360:angle
  j=fix(angle/45)

  bt=sqrt(by^2+bz^2)
  IF Bt GE 0 AND Bt LT 4 THEN return,vpc(1,*,j)
  IF Bt GE 4 AND Bt LT 6 THEN return,vpc(2,*,j)
  IF Bt GE 6 THEN return,vpc(3,*,j)

END
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
        Pcoeffs(k3)*alpha*L*cos(theta_prime(q))/sin(theta_prime(q))/Re
    endif
;
    k1 = (L LT Lmax) ? index_legendre(L+1,m,/dbl) : -1
    k2 = index_legendre(L,m,/dbl)
    if (k1 ge 0) then begin
      ecoeffs(k2,q) = ecoeffs(k2,q) $
        + Pcoeffs(k1)*alpha*(L+1+m)/sin(theta_prime(q))/Re
    endif
;
    if (m gt 0) then begin
      if (k3 ge 0) then k3 = k3 + 1
      k4 = k4 + 1
      if (k1 ge 0) then k1 = k1 + 1
      k2 = k2 + 1
      if (k3 ge 0) then begin
        ecoeffs(k4,q) = ecoeffs(k4,q) - $
          Pcoeffs(k3)*alpha*L*cos(theta_prime(q))/sin(theta_prime(q))/Re
      endif
      if (k1 ge 0) then begin
        ecoeffs(k2,q) = ecoeffs(k2,q) $
          + Pcoeffs(k1)*alpha*(L+1+m)/sin(theta_prime(q))/Re
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

;-------------------------------------------------------------
function perp_projection,vector,n
;
vn = reform(total(vector*n,1))
np = n_elements(vn)
vv = fltarr(2,np)
for i=0,np-1 do vv(*,i) = vn(i)*n(*,i)
return,vector - vv
end
;------------------------------------------------------------
;  This function calculates a "true" vector by combining the
;  real line-of-sight component from the doppler data with
;  the component perpendicular to the line-of-sight derived from
;  the potential map.
;
function calc_truevecs,vlos,kvect,pos,solution,latmin,lon_shft,lat_shft,order
;
; First shift coordinates into 'model' reference (pole shifted 4 deg nightwards)

posx=pos
if (lat_shft ne 0) then begin
  npnts=N_ELEMENTS(posx)/2
  kaz=FLTARR(npnts)
  crd_shft,lon_shft,lat_shft,npnts,posx,kaz
endif

tkvect=fltarr(2,n_elements(kvect))
FOR j=0,n_elements(kvect)-1 DO BEGIN
    tkvect(*,j) = [-cos(kvect(0,j)*!dtor ),sin(kvect(0,j)*!dtor )]
ENDFOR

; Calculate vectors
vvect = calc_vels(pos,solution,latmin,order)
tvect = perp_projection(vvect,tkvect)

for i=0,n_elements(vlos)-1 do begin
  tvect(*,i) = tvect(*,i) + vlos(i)*tkvect(*,i)
endfor
vmag = sqrt(tvect(0,*)^2 + tvect(1,*)^2)
vaz = atan(tvect(1,*),-tvect(0,*)) * !radeg

;
; if we have done a pole shift we now have to shift 
; back to the standard magnetic coordinate system
;
if (lat_shft NE 0) then begin
  npnts    =  n_elements(vmag)
  crd_shft,lon_shft,-lat_shft,npnts,pos,vaz
endif

return,{mag:vmag,az:vaz}
end
;-------------------------------------------------------------
;  This function calculates a "merge" vector by combining the
;  real line-of-sight components from the doppler data from 2 or more
;  overlapping radars.
;
function calc_mergevecs,vlos,kvect,pos,latmin,lon_shft,lat_shft
;
; First shift coordinates into 'model' reference (pole shifted 4 deg nightwards)

posx=pos
if (lat_shft ne 0) then begin
  npnts=N_ELEMENTS(posx)/2
  kaz=FLTARR(npnts)
  crd_shft,lon_shft,lat_shft,npnts,posx,kaz
endif

;Determine where we have overlapping vectors
got_data=where(abs(vlos) LT 9999.0, count)
merged=fltarr(2, 2000)
merged(*,*)=-9999.90
new_pos=fltarr(2, 2000)
new_pos(*,*)=-9999.90
vmag=fltarr(2000)
vmag(*)=-9999.90
vaz=fltarr(2000)
vaz(*)=-9999.90
m1=0
if count gt 0 then begin
 for i=0,n_elements(got_data)-1 do begin
  diff=pos(*,got_data)-rebin(pos(*,i),2,n_elements(got_data),/sample)
  comm=where(diff(0,*) EQ 0.0 AND diff(1,*) EQ 0.0, count)
  if count ge 2 then begin
   ;we have at least 2 overlapping vectors so calculate merged vectors:
   ;what permutations of vector pairs do we have?:
   vec_pairs=intarr(2,count^2)
   for k=0,count-1 do for l=0,count-1 do vec_pairs(*,(k*count)+l)=[k,l]
   vec_pairs(*,where(vec_pairs(0,*) ge vec_pairs(1,*)))=-1
   vec_pairs=vec_pairs(*,where(vec_pairs(0,*) ne -1))
   
   for j=0,n_elements(vec_pairs(0,*)) -1 do begin
    merged(*,m1+j) = [vlos(0,comm(vec_pairs(0,j)))*sin(kvect(0,comm(vec_pairs(0,j)))*!dtor ) + $
    	    	vlos(0,comm(vec_pairs(1,j)))*sin(kvect(0,comm(vec_pairs(1,j)))*!dtor ) , $
		vlos(0,comm(vec_pairs(0,j)))*cos(kvect(0,comm(vec_pairs(0,j)))*!dtor ) + $
    	    	vlos(0,comm(vec_pairs(1,j)))*cos(kvect(0,comm(vec_pairs(1,j)))*!dtor ) ]
    new_pos(*,m1+j) = pos(*,comm(vec_pairs(1,j)))
   endfor
   m1=m1+n_elements(vec_pairs(0,*))
   ;merged(*,j) = [-cos(kvect(0,comm(j))*!dtor ),sin(kvect(0,comm(j))*!dtor )]
    

; Calculate vectors
;vvect = calc_vels(pos,solution,latmin,order)
;tvect = perp_projection(vvect,tkvect)

;for i=0,n_elements(vlos)-1 do begin
;  tvect(*,i) = tvect(*,i) + vlos(i)*tkvect(*,i)
;endfor

;vmag = sqrt(tvect(0,*)^2 + tvect(1,*)^2)
;vaz = atan(tvect(1,*),-tvect(0,*)) * !radeg


  endif
 endfor
endif else begin
  print,'No data at current time'
  return,0
endelse

IF m1 Gt 0 THEN BEGIN
  vmag(0:m1-1) = sqrt(merged(0,0:m1-1)^2 + merged(1,0:m1-1)^2)
  vaz(0:m1-1) = atan(merged(0,0:m1-1),merged(1,0:m1-1)) * !radeg
ENDIF

;
; if we have done a pole shift we now have to shift 
; back to the standard magnetic coordinate system
;
if (lat_shft NE 0) then begin
  npnts    =  n_elements(vmag)
  crd_shft,lon_shft,-lat_shft,npnts,pos,vaz
endif

return,{mag:vmag,az:vaz,pos:new_pos}
end
;-------------------------------------------------------------
function calc_losvecs,model=model

COMMON new_to_old

IF keyword_set(model) THEN return,{mag:mdata(3,*),az:mdata(2,*)} ELSE $
return,{mag:vdata(6,*),az:vdata(2,*)}
end
;-------------------------------------------------------------
function get_electric_potentials,in_data,index=index,ephi=ephi,etheta=etheta,etot=etot

COMMON new_to_old

IF ~KEYWORD_SET(map_open_flag) AND N_ELEMENTS(in_data) EQ 0 THEN BEGIN
  print,'You have to read in the MAP file first!'
  print,'OPEN_MAP_POT,fname   (where fname includes the directory string if not in current directory)'
  return,-1
ENDIF

IF KEYWORD_SET(index) THEN set_up_arrays,index
IF ~KEYWORD_SET(time_spec_flag) AND N_ELEMENTS(in_data) EQ 0 THEN set_up_arrays,0

IF N_ELEMENTS(in_data) GT 0 THEN BEGIN
  latzer_ref=in_data[0].latmin
;  solution=transpose(rebin(in_data.coeffs,49,3))
; think this must be from using some other form of in_data besides map_data...
  solution=in_data[0].coeff
  lat_shft = 0
  order = 6
ENDIF

plot_lat_min=latzer_ref
latmin=latzer_ref


latstep   =  1.0				;points in plotting coords
longstep  =  2.0
nlats     =  fix((90.-plot_lat_min)/latstep)
nlongs    =  fix(360./longstep)
lats      =  findgen(nlats)*latstep + plot_lat_min
longs     =  findgen(nlongs+1)*longstep
grid      =  create_grid(plot_lat_min, latstep, longstep)

xxx_arr   =  reform(grid(1,*))
zon_arr   = [xxx_arr(uniq(xxx_arr)),360.]
zat_arr   =  reform(grid(0,0:nlats-1))

pot_arr1   =  fltarr(nlongs+1,nlats)

if (lat_shft eq 0) then iflg_coord = 1
if (lat_shft ne 0) then iflg_coord = 0

if (iflg_coord eq 1) then begin			;plot in primed coords

  convert_pos,grid,th,ph

  tmax      = (90.0 - latmin)*!dtor
  tprime    =  norm_theta(th,tmax)
  x         =  cos(tprime)
  plm       =  eval_legendre(order, x)
  v         =  eval_potential(solution(2,*),plm,ph)
  IF KEYWORD_SET(ephi) THEN BEGIN
     epc    =  eval_ephi_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     v		=  eval_component(epc,plm,grid,Lmax=order)*1e6
  ENDIF
  IF KEYWORD_SET(etheta) THEN BEGIN
     etc    =  eval_etheta_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     v		=  eval_component(etc,plm,grid,Lmax=order)*1e6
  ENDIF
  IF KEYWORD_SET(etot) THEN BEGIN
     etc    =  eval_etheta_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     epc    =  eval_ephi_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     vt		=  eval_component(etc,plm,grid,Lmax=order)*1e6
     vp		=  eval_component(epc,plm,grid,Lmax=order)*1e6
     v		=  sqrt(vt^2+vp^2)
  ENDIF

  for i = 0,nlongs-1 do begin
    jl  =   i*nlats
    ju  =   i*nlats + (nlats-1)
    pot_arr1(i,*) = v(jl:ju)/1000.
  endfor

  pot_arr1(nlongs,*) =  pot_arr1(0,*)

  q = where(zat_arr le latmin, qc)		;set to zero below latmin
  if (qc ne 0) then pot_arr1(*,q(0):q(qc-1)) = 0.

endif

if (iflg_coord eq 0) then begin			;plot in unprimed coords

  npnts      =  nlongs * nlats			
  kaz        =  fltarr(npnts)
  crd_shft,lon_shft,lat_shft,npnts,grid,kaz

  xon_arr_p  =  fltarr(nlongs+1,nlats)
  xat_arr_p  =  fltarr(nlongs+1,nlats)
  for i = 0,nlongs-1 do begin			;find prime coords (latmin chk)
    jl  = i*nlats
    ju  = jl + nlats-1
    xon_arr_p(i,*) = grid(1,jl:ju)
    xat_arr_p(i,*) = grid(0,jl:ju)
  endfor

  convert_pos,grid,th,ph

  tmax      = (90.0 - latmin)*!dtor
  tprime    =  norm_theta(th,tmax)
  x         =  cos(tprime)
  plm       =  eval_legendre(order, x)
  v         =  eval_potential(solution(2,*),plm,ph)
  IF KEYWORD_SET(ephi) THEN BEGIN
     epc    =  eval_ephi_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     v		=  eval_component(epc,plm,grid,Lmax=order)*1e6
  ENDIF
  IF KEYWORD_SET(etheta) THEN BEGIN
     etc    =  eval_etheta_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     v		=  eval_component(etc,plm,grid,Lmax=order)*1e6
  ENDIF
  IF KEYWORD_SET(etot) THEN BEGIN
     etc    =  eval_etheta_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     epc    =  eval_ephi_coeffs(solution(2,*),th,Lmax=order,latmin=latmin)
     vt		=  eval_component(etc,plm,grid,Lmax=order)*1e6
     vp		=  eval_component(epc,plm,grid,Lmax=order)*1e6
     v		=  sqrt(vt^2+vp^2)
  ENDIF
  
  for i = 0,nlongs-1 do begin
    jl  =   i*nlats
    ju  =   i*nlats + (nlats-1)
    pot_arr1(i,*) = v(jl:ju)/1000.
  endfor

  q = where(xat_arr_p le latmin, qc)		;zero pot below latmin
  if (qc ne 0) then pot_arr1(q) = 0.

  pot_arr1(nlongs,*)      =  pot_arr1(0,*)

endif

return, {potarr:pot_arr1, zatarr:zat_arr, zonarr:zon_arr}

end

;-------------------------------------------------------------

function convert_electric_potentials,pot_data,polar_x,polar_y

COMMON new_to_old

pot_arr=pot_data.potarr
zat_arr=pot_data.zatarr
zon_arr=pot_data.zonarr

; Convert to polar grid
polar_x=FLTARR(n_elements(zon_arr),n_elements(zat_arr))
polar_y=FLTARR(n_elements(zon_arr),n_elements(zat_arr))

FOR j=0,n_elements(zat_arr)-1 DO BEGIN
	FOR i=0,n_elements(zon_arr)-1 DO BEGIN
		polar_x(i,j)=-(90-zat_arr(j))*SIN((zon_arr(i)+lon_shft)*!dtor)
		polar_y(i,j)= (90-zat_arr(j))*COS((zon_arr(i)+lon_shft)*!dtor)
	ENDFOR
ENDFOR

RETURN,1

end

;-------------------------------------------------------------

function convert_electric_potentials_rect,pot_data,polar_x,polar_y,polar_p2

COMMON new_to_old

pot_arr=pot_data.potarr
zat_arr=pot_data.zatarr
zon_arr=pot_data.zonarr

; Convert to polar grid
polar_x=2*(90-latzer_ref)*findgen(51)/50.-(90-latzer_ref)
polar_y=2*(90-latzer_ref)*findgen(51)/50.-(90-latzer_ref)
polar_p=FLTARR(51,51,2)

FOR j=0,n_elements(zat_arr)-1 DO BEGIN
	FOR i=0,n_elements(zon_arr)-1 DO BEGIN
		px =-(90-zat_arr(j))*SIN((zon_arr(i)+lon_shft)*!dtor)
		py = (90-zat_arr(j))*COS((zon_arr(i)+lon_shft)*!dtor)
		res = min(abs(polar_x-px),wpx)
		res = min(abs(polar_y-py),wpy)
		polar_p(wpx,wpy,1)=polar_p(wpx,wpy,1)+1
		polar_p(wpx,wpy,0)=((polar_p(wpx,wpy,1)-1)*polar_p(wpx,wpy,0)+pot_arr(i,j))/polar_p(wpx,wpy,1)
	ENDFOR
ENDFOR
polar_p2=polar_p(*,*,0)

RETURN,1

end

;--------------------------------------------------------------

function read_map,fname,map_data,size=size

COMMON new_to_old

ON_IOERROR, BAD

openr,unit,fname,/get_lun   	;open the MAP file, include 
    	    	    	    	;directory in fname if 
    	    	    	    	;not in correct directory

blah=' ' & source=' ' & mod_dir=' ' & mod_mag=' '     ;define strings

;read upto the point we get npnt
    readf,unit,format='(12i0)',syr,smo,sdy,shr,smt,ssc,eyr,emo,edy,ehr,emt,esc
    readf,unit,nblocks
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,format='(a10,1x,i10,1x,i10)',source, major_rev, minor_rev
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,bx,by,bz
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,format='(i10,1x,i10,1x,i10,1x,a10,1x,a10)',doping_level,imf_flag,imf_delay,mod_dir,mod_mag
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,lon_shft,lat_shft,hemisphere,order,latmin,error_wt,model_wt
    
    ;lon_shft doesn't seem to contain the longitude shift due to UT, 
    ;so calculate it manually:
    lon_shft=MLT(syr,cnvtime(syr,smo,sdy,shr,smt,ssc),0)*15+180
    
    readf,unit,num_coeff
    IF num_coeff GT 300 THEN BEGIN
      print,'Number of Coefficients is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    coeff=fltarr(4,num_coeff > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF num_coeff GT 0 THEN BEGIN
      FOR i=0,num_coeff-1 DO BEGIN
        readf,unit,n,n1,n2,n3
        coeff(*,i)=[n,n1,n2,n3]
      ENDFOR
    ENDIF
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,chi_sqr,chi_sqr_dat,rms_err
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,pot_drop,pot_drop_err
    readf,unit,pot_max,pot_max_err
    readf,unit,pot_min,pot_min_err
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,mlt_start,mlt_end,mlt_avg
    
    readf,unit,num_model
    IF num_model GT 2000 THEN BEGIN
      print,'Number of model points is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    mdata=fltarr(6,num_model)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    FOR i=0,num_model-1 DO BEGIN
      readf,unit,gmlomg,gmlat,kvect,vlos
      mdata(0:3,i)=[gmlomg,gmlat,kvect,vlos]
    ENDFOR
    
    readf,unit,num_bnd
    IF num_bnd GT 73 THEN BEGIN
      print,'Number of boundary points is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    bnd=fltarr(2,num_bnd > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF num_bnd EQ 0 THEN bnd(*,0)=-9999.9 ELSE BEGIN
     FOR i=0,num_bnd-1 DO BEGIN
      readf,unit,bnd_lon,bnd_lat
      bnd(*,i)=[bnd_lon,bnd_lat]
     ENDFOR
    ENDELSE
    
    readf,unit,snpnt,snprm
    sdata=fltarr(snprm,snpnt > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF snpnt EQ 0 THEN sdata(*,0)=-9999.9 ELSE BEGIN
     FOR i=0,snpnt-1 DO BEGIN
      readf,unit,st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max 
      sdata(*,i)=[st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max]
     ENDFOR
    ENDELSE
    
    readf,unit,npnt,nprm
    data_elements=npnt > 2000
close,unit
free_lun,unit

IF keyword_set(size) THEN data_elements=size
IF data_elements GT 32000 THEN npnt1=0L ELSE npnt1=0

openr,unit,fname,/get_lun   	;open the MAP file, include 

;define a structure which will contain all map_parameters

  r = { syr:0, smo:0, sdy:0, shr:0, smt:0, ssc:0, $
    	eyr:0, emo:0, edy:0, ehr:0, emt:0, esc:0, $
	nblocks:0, source:' ', major_rev:0, minor_rev:0, $
	bx:0.0, by:0.0, bz:0.0, doping_level:0, imf_flag:0, $
	imf_delay:0, mod_dir:' ', mod_mag:' ', lon_shft:0.0, $
	lat_shft:0.0, hemisphere:0, order:0, latmin:0.0, $
	error_wt:0, model_wt:0, num_coeff:0, coeff:fltarr(4,300), $
	chi_sqr:0.0, chi_sqr_dat:0.0, rms_err:0.0, pot_drop:0.0, $
	pot_drop_err:0.0, pot_max:0.0, pot_max_err:0.0, pot_min:0.0, $
	pot_min_err:0.0, mlt_start:0.0, mlt_end:0.0, mlt_avg:0.0, $
	num_model:0, mdata:fltarr(6,2000), num_bnd:0, bnd:fltarr(2,73), $
	snpnt:0, snprm:0, sdata:fltarr(18,50), npnt:npnt1, nprm:0, $
	data:fltarr(12,data_elements) }

;map_data = replicate(r,1441)

j=-1
WHILE NOT EOF(unit) DO BEGIN

    r.coeff(*,*)=9999.9
    r.mdata(*,*)=9999.9
    r.sdata(*,*)=9999.9
    r.data(*,*)=9999.9

    j=j+1

    readf,unit,format='(12i0)',syr,smo,sdy,shr,smt,ssc,eyr,emo,edy,ehr,emt,esc
    readf,unit,nblocks
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,format='(a10,1x,i10,1x,i10)',source, major_rev, minor_rev
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,bx,by,bz
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,format='(i10,1x,i10,1x,i10,1x,a10,1x,a10)',doping_level,imf_flag,imf_delay,mod_dir,mod_mag
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,lon_shft,lat_shft,hemisphere,order,latmin,error_wt,model_wt
    
    ;lon_shft doesn't seem to contain the longitude shift due to UT, 
    ;so calculate it manually:
    lon_shft=MLT(syr,cnvtime(syr,smo,sdy,shr,smt,ssc),0)*15+180
    
    readf,unit,num_coeff
    IF num_coeff GT 300 THEN BEGIN
      print,'Number of Coefficients is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    coeff=fltarr(4,num_coeff > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF num_coeff GT 0 THEN BEGIN
      FOR i=0,num_coeff-1 DO BEGIN
        readf,unit,n,n1,n2,n3
        coeff(*,i)=[n,n1,n2,n3]
      ENDFOR
    ENDIF
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,chi_sqr,chi_sqr_dat,rms_err
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,pot_drop,pot_drop_err
    readf,unit,pot_max,pot_max_err
    readf,unit,pot_min,pot_min_err
    
    FOR i=0,3 DO readf,unit,format='(a1)',blah
    readf,unit,mlt_start,mlt_end,mlt_avg
    
    readf,unit,num_model
    IF num_model GT 2000 THEN BEGIN
      print,'Number of model points is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    mdata=fltarr(6,num_model)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    FOR i=0,num_model-1 DO BEGIN
      readf,unit,gmlomg,gmlat,kvect,vlos
      mdata(0:3,i)=[gmlomg,gmlat,kvect,vlos]
    ENDFOR
    
    readf,unit,num_bnd
    IF num_bnd GT 73 THEN BEGIN
      print,'Number of boundary points is too high'
      print,'Going to bugger up structure complilation'
      return,0
    ENDIF
    bnd=fltarr(2,num_bnd > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF num_bnd EQ 0 THEN bnd(*,0)=-9999.9 ELSE BEGIN
     FOR i=0,num_bnd-1 DO BEGIN
      readf,unit,bnd_lon,bnd_lat
      bnd(*,i)=[bnd_lon,bnd_lat]
     ENDFOR
    ENDELSE
    
    readf,unit,snpnt,snprm
    sdata=fltarr(snprm,snpnt > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF snpnt EQ 0 THEN sdata(*,0)=-9999.9 ELSE BEGIN
     FOR i=0,snpnt-1 DO BEGIN
      readf,unit,st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max 
		
      sdata(*,i)=[st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max]
     ENDFOR
    ENDELSE
    
    readf,unit,npnt,nprm
    data=fltarr(nprm,npnt >1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF npnt EQ 0 THEN data(*,0)=0.0 ELSE BEGIN
     IF nprm LE 8 THEN BEGIN
      FOR i=0,npnt-1 DO BEGIN
        readf,unit,gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd 
        data(*,i)=[gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd]
      ENDFOR
     ENDIF ELSE BEGIN
      FOR i=0,npnt-1 DO BEGIN
      	readf,unit,gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd,pwr,pwr_sd,wdt,wdth_sd
        data(*,i)=[gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd,pwr,pwr_sd,wdt,wdth_sd]
      ENDFOR    
     ENDELSE
    ENDELSE
    
    
    ;Fill up structure "r"
    
    r.syr = syr & r.smo = smo & r.sdy = sdy & r.shr = shr & r.smt = smt & r.ssc = ssc
    r.eyr = eyr & r.emo = emo & r.edy = edy & r.ehr = ehr & r.emt = emt & r.esc = esc
    r.nblocks = nblocks & r.source = source & r.major_rev = major_rev
    r.minor_rev = minor_rev & r.bx = bx & r.by = by & r.bz = bz
    r.doping_level = doping_level & r.imf_flag = imf_flag & r.imf_delay = imf_delay
    r.mod_dir = mod_dir & r.mod_mag = mod_mag & r.lon_shft = lon_shft
    r.lat_shft = lat_shft & r.hemisphere = hemisphere & r.order = order 
    r.latmin = latmin & r.error_wt = error_wt & r.model_wt = model_wt
    r.num_coeff = num_coeff & r.chi_sqr = chi_sqr 
    r.chi_sqr_dat = chi_sqr_dat & r.rms_err = rms_err & r.pot_drop = pot_drop
    r.pot_drop_err = pot_drop_err & r.pot_max = pot_max & r.pot_max_err = pot_max_err
    r.pot_min = pot_min & r.pot_min_err = pot_min_err & r.mlt_start = mlt_start
    r.mlt_end = mlt_end & r.mlt_avg = mlt_avg & r.num_model = num_model
    r.num_bnd = num_bnd & r.snpnt = snpnt
    r.snprm = snprm & r.npnt = npnt & r.nprm = nprm
    
    FOR k=0,3 DO BEGIN
      r.coeff(k,0:n_elements(coeff(0,*))-1) = coeff(k,*)
    ENDFOR
    
    FOR k=0,3 DO BEGIN
      r.mdata(k,0:n_elements(mdata(0,0:num_model-1))-1) = mdata(k,0:num_model-1)
    ENDFOR
    r.num_model = num_model
    
    FOR k=0,1 DO BEGIN
      r.bnd(k,0:n_elements(bnd(0,*))-1) = bnd(k,*)
    ENDFOR

    FOR k=0,17 DO BEGIN
      r.sdata(k,0:n_elements(sdata(0,*))-1) = sdata(k,*)
    ENDFOR
    
    FOR k=0,nprm-1 DO BEGIN
      r.data(k,0:n_elements(data(0,*))-1) = data(k,*)
    ENDFOR

;    map_data(j)=r     ;add "r" to day-structure "map_data"
	IF j EQ 0 THEN map_data=r ELSE map_data=[map_data,r]    ;build map_data as you go...

ENDWHILE
close,unit
free_lun,unit

map_data=map_data(0:j)      ;resize "map_data"
map_open_flag=1
vector_flag=1

return,0

BAD:  return,-1

END

;-------------------------------------------------------------
function read_grid,gname,size=size

openr,unit,gname,/get_lun   	;open the GRID file, include 
    	    	    	    	;directory in fname if 
    	    	    	    	;not in correct directory

; Declare error label.  
ON_IOERROR, BAD  

blah=' ' & source=' ' & mod_dir=' ' & mod_mag=' '     ;define strings

;read upto the point we get npnt
    readf,unit,format='(12i0)',syr,smo,sdy,shr,smt,ssc,eyr,emo,edy,ehr,emt,esc
    readf,unit,nblocks
    readf,unit,snpnt,snprm
    sdata=fltarr(snprm,snpnt > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF snpnt EQ 0 THEN sdata(*,0)=-9999.9 ELSE BEGIN
     FOR i=0,snpnt-1 DO BEGIN
      readf,unit,st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max 
      sdata(*,i)=[st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max]
     ENDFOR
    ENDELSE
    
    readf,unit,npnt,nprm
    data_elements=npnt > 2000
    IF KEYWORD_SET(size) THEN data_elements=size
    IF data_elements GT 32000 THEN npnt1=0L ELSE npnt1=0

close,unit
free_lun,unit
openr,unit,gname,/get_lun

;define a structure which will contain all map_parameters

  r = { syr:0, smo:0, sdy:0, shr:0, smt:0, ssc:0, $
    	eyr:0, emo:0, edy:0, ehr:0, emt:0, esc:0, $
	nblocks:0, snpnt:0, snprm:0, sdata:fltarr(18,50), $
	npnt:npnt1, nprm:0, data:fltarr(12,data_elements) }

grid_data = replicate(r,1441)

j=-1
WHILE NOT EOF(unit) DO BEGIN

    r.sdata(*,*)=9999.9
    r.data(*,*)=9999.9

    j=j+1

    readf,unit,format='(12i0)',syr,smo,sdy,shr,smt,ssc,eyr,emo,edy,ehr,emt,esc
    readf,unit,nblocks
    readf,unit,snpnt,snprm
    sdata=fltarr(snprm,snpnt > 1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF snpnt EQ 0 THEN sdata(*,0)=-9999.9 ELSE BEGIN
    
     point_lun,-1*unit,lunpos
     FOR i=0,snpnt-1 DO BEGIN
      valid=0
      readf,unit,st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max
      valid=1
      sdata(*,i)=[st_id,chn,nvec,freq0,major_rev,minor_rev,prog_id, $
      	    	    noise_mean,noise_sd,gsct,v_min,v_max,p_min,p_max, $
		    w_min,w_max,ve_min,ve_max]
     ENDFOR
     BAD: IF ~ valid THEN BEGIN
	   point_lun,unit,lunpos
       temp_string=''
       FOR i=0,snpnt-1 DO BEGIN
         readf,unit,temp_string
	     ast=strpos(temp_string,'*')
	     eplus=strpos(temp_string,'e+')
	     poss=strsplit(temp_string)
	     IF ast(0) NE -1 OR eplus(0) NE -1 THEN BEGIN
	       FOR z=0,6 DO sdata(z,i)=float(strmid(temp_string,poss(z),poss(z+1)-poss(z)))
	       sdata(7:17,i)=[10000.,10000., 1,35.0000,2000.00,3.00000,50.0000,10.0000,1000.00,0.00000,200.000]
	     ENDIF ELSE BEGIN
	       FOR z=0,16 DO sdata(z,i)=float(strmid(temp_string,poss(z),poss(z+1)-poss(z)))
	       sdata(17,i)=float(strmid(temp_string,poss(17),strlen(temp_string)-1))
	     ENDELSE
       ENDFOR
     END

    ENDELSE
    
    readf,unit,npnt,nprm
    data=fltarr(nprm,npnt >1)
    FOR i=0,2 DO readf,unit,format='(a1)',blah
    IF npnt EQ 0 THEN data(*,0)=0.0 ELSE BEGIN
     IF nprm LE 8 THEN BEGIN
      FOR i=0,npnt-1 DO BEGIN
        readf,unit,gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd 
        data(*,i)=[gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd]
      ENDFOR
     ENDIF ELSE BEGIN
      FOR i=0,npnt-1 DO BEGIN
      	readf,unit,gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd,pwr,pwr_sd,wdt,wdth_sd
        data(*,i)=[gmlong,gmlat,kvect,st_id,chn,grid_index,vlos,vlos_sd,pwr,pwr_sd,wdt,wdth_sd]
      ENDFOR    
     ENDELSE
    ENDELSE
    
    
    ;Fill up structure "r"
    
    r.syr = syr & r.smo = smo & r.sdy = sdy & r.shr = shr & r.smt = smt & r.ssc = ssc
    r.eyr = eyr & r.emo = emo & r.edy = edy & r.ehr = ehr & r.emt = emt & r.esc = esc
    r.nblocks = nblocks & r.snpnt = snpnt & r.snprm = snprm
    r.npnt = npnt & r.nprm = nprm
    
    FOR k=0,17 DO BEGIN
      r.sdata(k,0:n_elements(sdata(0,*))-1) = sdata(k,*)
    ENDFOR
    
    FOR k=0,nprm-1 DO BEGIN
      r.data(k,0:n_elements(data(0,*))-1) = data(k,*)
    ENDFOR

    grid_data(j)=r     ;add "r" to day-structure "map_data"

ENDWHILE
close,unit
free_lun,unit

grid_data=grid_data(0:j)      ;resize "grid_data"
return,grid_data

END

;-------------------------------------------------------------------------------

FUNCTION vel_time_series,lats,lons,mlt=mlt,polar=polar,max=max

;returns a 2-D array (2,n) containing velocity time series for current file (n records)
;determined from the average in a box constrained by lats and lons (2 2-elements arrays)
;if mlt keyword is set then specify mlts in place of lons
;components returned are Vnorth,Veast unless polar keyword is set, then returns Vmag,Vaz

COMMON map_pot
COMMON new_to_old

m=n_elements(map_data)
vel_series1=fltarr(2,m)

FOR j=0,m-1 DO BEGIN
  
  IF keyword_set(mlt) THEN lons1=[mlt_to_maglon(lons(0),map_data(j).shr+map_data(j).smt/60.),mlt_to_maglon(lons(1),map_data(j).shr+map_data(j).smt/60.)] ELSE lons1=lons
  box=where(map_data(j).data(0,*) GE lons1(0) AND map_data(j).data(0,*) LT lons1(1) AND $
    map_data(j).data(1,*) GE lats(0) AND map_data(j).data(1,*) LT lats(1))
  IF lons1(1) LT lons1(0) THEN $
  box=where((map_data(j).data(0,*) GE lons1(0) AND map_data(j).data(0,*) LT 999.9) OR map_data(j).data(0,*) LT lons1(1) AND $
    map_data(j).data(1,*) GE lats(0) AND map_data(j).data(1,*) LT lats(1))

  IF n_elements(box) GE 2 THEN BEGIN
    IF vector_flag EQ 1 THEN velt=find_gradV([map_data(j).data(1,*),map_data(j).data(0,*)],map_data(j).coeff,map_data(j).latmin,map_data(j).lon_shft,map_data(j).lat_shft,map_data(j).order)
    IF vector_flag EQ 2 THEN velt=calc_truevecs(map_data(j).data(6,*),map_data(j).data(2,*),[map_data(j).data(1,*),map_data(j).data(0,*)],map_data(j).coeff,map_data(j).latmin,map_data(j).lon_shft,map_data(j).lat_shft,map_data(j).order)
    IF vector_flag EQ 3 THEN velt=calc_losvecs()
    IF keyword_set(max) THEN vel_series1(*,j)=[max(velt.mag(box)*cos(velt.az(box)*!dtor)), max(velt.mag(box)*sin(velt.az(box)*!dtor))] ELSE $
    vel_series1(*,j)=[mean(velt.mag(box)*cos(velt.az(box)*!dtor)), mean(velt.mag(box)*sin(velt.az(box)*!dtor))] ;NORTH,EAST
    IF keyword_set(polar) THEN $
    vel_series1(*,j)=[sqrt(vel_series1(0,j)^2 + vel_series1(1,j)^2),atan(vel_series1(1,j),vel_series1(0,j))*!radeg]
  ENDIF

ENDFOR

RETURN,vel_series1
END

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;F
;F FUNCTION B_TRACE
;F
;F WRITTEN:	08/11/00		
;F MODIFIED:	22/01/01	
;F AUTHOR: 	Jim Wild	
;F COMMENTS:	Function that uses a c wrapper to call the T89c or T98 Tsyganenko 
;F		model Fortran code to calculate the footprint position of a
;F		spacecraft
;F
;F REFERENCES:	
;F
;F INPUTS:	xgse,ygse,zgse:	s/c GSE coords in km
;F		year: 	eg, 2000
;F		dayno:  the daynumber of the year
;F		hour:   ut hour of the day
;F		min:   	ut min of the day	
;F		sc:   	ut sc of the day
;F		input_param:	if this is a single value then it is taken as the
;F				desired kp value and the T89 model is used. If it
;F				is a 4 element array, the T96 model is used and
;F				these are taken as pdyn,dst,By and Bz in that order. 
;F				In this case:
;F				pdyn - solar wind dynamic pressure in nPa (ideally 
;F				       between 0.5 and 10)
;F				dst - Dst index value (ideadlly between -100 and +20)
;F                               BYIMF and BZIMF - IMF components in nT (ideally between 
;F							-10 and +10)
;F				If input_param is a 6 element array then the T01 model 
;F				is used. The first 4 parameters are as the T96 case above.
;F				the 5th and 6th elements of the array represent the g1 and
;F				g2 indices for the model. See the Tsyganenlo 2001 papers for 
;F				details
;F 
;F		alt:	altitude above the surface of the Earth (km) at which to 
;F			terminate the fieldline tracing
;F		dir:	direction of trace. -1.0 for northern hemisphere, 1.0 for
;F			southern
;F		/file:  option that if specifies creates a file B_TRACE.DAT containing
;F			the GSE coordinates of the last fieldline trace
;F
;F
;F RETURNS:	A eight element array containing geog lat, geog lon, mag lat (PACE),
;F		mag lon, MLT, final trace position radial distance from the centre of 
;F		the Earth, number of points invoved in the trace and the L-value of
;F		the footprint position in that order  (** ignore this L-value **) 
;F
;F CALLING SEQUENCE:	result=B_TRACE(xgse,ygse,zgse,year,dayno,hour,min,sc,input_param,alt,dir,/file)
;------------------------------------------------------------------------------

   FUNCTION B_TRACE,xgse,ygse,zgse,year,dayno,hour,min,sc,input_param,alt,dir,file=file

; Create the output array and fill the lat and lon values with floating point values
; ready for output from the c wrapper and Fortran code  

   	pos=fltarr(8)
   	np=-9999.0
   	out_r=-9999.0
   	lat=-99999.0
   	lon=-99999.0
   	o_file=0.0 
	IF KEYWORD_SET(file) THEN BEGIN
		result=FINDFILE('B_TRACE.DAT',count=count)
		IF count gt 0 THEN SPAWN,'/bin/rm B_TRACE.DAT'		
		o_file=1.0
	ENDIF	
 
; Make sure all arguement that are expected as real by the Fortran are in fact 
; floating point  	
   	xgse=FLOAT(xgse)
   	ygse=FLOAT(ygse)
   	zgse=FLOAT(zgse)
   	alt=FLOAT(alt)
   	dir=FLOAT(dir)

; Make sure all arguement that are expected as int by the c-wrapper are in fact 
; longword integers   	
   	year=LONG(year)
   	dayno=LONG(dayno)
  	hour=LONG(hour)
  	min=LONG(min)
  	sc=LONG(sc)
   	o_file=LONG(o_file) 	

; Need to select the appropriate model and make sure that the input parameters are in the 
; correct format
	IF N_ELEMENTS(input_param) EQ 1 THEN BEGIN
		;PRINT,'***Calling Tsyganenko 89 model***'
  		kp=LONG(input_param(0))
  		status=CALL_EXTERNAL('/Users/ag27/Idl/geopack2008/btrace89.so','BTRACE89_idl',$
  			xgse,ygse,zgse,year,dayno,hour,min,sc,kp,alt,dir,o_file,lat,lon,out_r,np)	
  	ENDIF
  	IF N_ELEMENTS(input_param) EQ 4 THEN BEGIN
		;PRINT,'***Calling Tsyganenko 96 model***'
		pdyn=FLOAT(input_param(0))
		dst=FLOAT(input_param(1))
		byimf=FLOAT(input_param(2))
		bzimf=FLOAT(input_param(3))
		status=CALL_EXTERNAL('/Users/ag27/Idl/geopack2008/btrace96.so','BTRACE96_idl',$
  			xgse,ygse,zgse,year,dayno,hour,min,sc,pdyn,dst,byimf,bzimf,alt,dir,o_file,lat,lon,out_r,np)
	ENDIF
	IF N_ELEMENTS(input_param) EQ 6 THEN BEGIN
		;PRINT,'***Calling Tsyganenko 01 model***'
		pdyn=FLOAT(input_param(0))
		dst=FLOAT(input_param(1))
		byimf=FLOAT(input_param(2))
		bzimf=FLOAT(input_param(3))
		g1=FLOAT(input_param(4))
		g2=FLOAT(input_param(5))
		status=CALL_EXTERNAL('/Users/ag27/Idl/geopack2008/btrace01.so','BTRACE01_idl',$
  			xgse,ygse,zgse,year,dayno,hour,min,sc,pdyn,dst,byimf,bzimf,g1,g2,alt,dir,o_file,lat,lon,out_r,np)
	ENDIF

; Put the geog lat and lon into the output array     
   	pos(0)=lat
   	pos(1)=lon
   	
; Calculate the magnetic positions... 	
   	mag_pos=cnvcoord(lat,lon,alt)
   	pos(2)=mag_pos(0)
   	pos(3)=mag_pos(1)
   	
; ...and the MLT   	
   	ut_secs=((dayno-1.)*86400.)+(hour*3600.)+(min*60.)+sc
   	pos(4)=mlt(year,ut_secs,mag_pos(1))

; and the output alltitude of the final trace position
	pos(5)=out_r

; and the number of steps in the trace
	pos(6)=np
	
; and the L value of the footprint position
	pos(7)=1./((COS(pos(2)*!DTOR))^2.0)	
;	print,pos 	
; Now return the output
   	RETURN,pos
   	END
   
;------------------------------------------------------------------------------
;-------------------------------------------------------------------------------

PRO convert_pos,pos,theta,phi
theta=reform(pos(0,*))
phi=reform(pos(1,*))
theta = (90.-theta)*!dtor
phi = phi*!dtor
return
END

;----------------------------------------------------------------------

PRO legendre_k_to_lm_index,k,l,m,sh=sh,sngl=sngl,dble=dble

if (keyword_set(sh)) then begin
  l = fix(sqrt(k))
  m = k - (l*(l+1))
  return
endif $
else if (keyword_set(sngl)) then begin
  l = fix (-0.5 + 0.5*sqrt(1.0 + 8.0*k))
  m = k - (l*(l+1))/2
  return
endif $
else if (keyword_set(dble)) then begin
  L = fix(sqrt(float(k)))
  if (L*L eq k) then m = 0 else m = (k + 1 - L*L)/2
  if (index_legendre(L,m,/dbl) ne k) then m = -m
  return
endif $
else begin
  print,'k_to_lm_index:  You must specify one and only one of'
  print,'the keywords /SH (spherical harmonics),'
  print,'/SNGL (single Legendre polynomial),'
  print,'/DBL (cos(phi),sin(phi) pairs with Legendre Polynomials'
  return
endelse
END

;-------------------------------------------------------------------------

PRO crd_shft,lon_shft,lat_shft,npnts,pos,kaz

  lat   =  fltarr(npnts)
  lon   =  fltarr(npnts)
  vazm  =  fltarr(npnts)
  lat   =  pos(0,0:npnts-1)
  lon   =  pos(1,0:npnts-1)
  vazm  =  kaz(0:npnts-1)

  lon_min  =    0.
  lon_max  =  360.

  lon      =  lon + lon_shft
  check_lon,lon_min,lon_max,lon 

  if (lat_shft eq 0.) then begin
    print,' lat_shift is set to zero!'
    print,' this call to crd_shft is unneccesary..but continuing..'
  endif else begin 

    a_side = (90.-lat) * !pi/180.				;colat
    B_angl = (180.-lon) * !pi/180.				;lon
    d_side =  lat_shft * !pi/180.				;shift

    arg    =  cos(a_side)*cos(d_side) + sin(a_side)*sin(d_side)*cos(B_angl)
    q = where (abs(arg) gt 1., qcnt)
    if (qcnt ne 0) then arg(q) = arg(q)/abs(arg(q)) 

    b_side =  acos(arg)


    q = where(b_side eq 0., qc)			;adjust for point on pole
    if (qc ne 0) then b_side(q) = 0.1

    arg    = (cos(a_side)-cos(b_side)*cos(d_side))/(sin(b_side)*sin(d_side))

    q      =  where (abs(arg) gt 1., qcnt)
    if (qcnt ne 0) then arg(q) = arg(q)/abs(arg(q)) 

    A_angl =   acos(arg)

    q      =  where (lon gt 180., qcnt)			;acos ambiguity

    if (qcnt ne 0) then A_angl(q) = 2*!pi - A_angl(q)

    lon    =  A_angl * 180/!pi
    lat    =  (!pi/2.-b_side) * 180./!pi
  

    arg    = (cos(d_side)-cos(a_side)*cos(b_side))/(sin(a_side)*sin(b_side))
    q      =  where (abs(arg) gt 1., qcnt)
    if (qcnt ne 0) then arg(q) = arg(q)/abs(arg(q)) 
    C_angl =  acos(arg)

    sign_d =  d_side/abs(d_side)

    q      =  where (lon le 180., qcnt)			
    if (qcnt ne 0) then vazm(q) = vazm(q) - sign_d*C_angl(q)*180./!pi

    q      =  where (lon gt 180., qcnt)			
    if (qcnt ne 0) then vazm(q) = vazm(q) + sign_d*C_angl(q)*180./!pi

  endelse

  lon    =  lon - lon_shft
  check_lon,lon_min,lon_max,lon 

;  q = where(lon ge 180., qc)
;  if (qc ne 0) then lon(q) = lon(q) - 360.
;  q = where(lon lt -180., qc)
;  if (qc ne 0) then lon(q) = lon(q) + 360.

  pos(0,0:npnts-1) = lat(*)
  pos(1,0:npnts-1) = lon(*)
  kaz(*)           = vazm(*)

return
END

;-------------------------------------------------------------

PRO check_lon,xl,xu,lon

q = where(lon lt xl, qc)
if (qc ne 0) then lon(q) = lon(q) + 360.

q = where(lon gt xu, qc)
if (qc ne 0) then lon(q) = lon(q) - 360.

if (min(lon) lt xl or max(lon) gt xu) then begin
  print,' Invalid longitude = ',min(lon),max(lon)
  stop
endif

return
END

;----------------------------------------------------------------------

PRO plot_col_bar,ymaps,ymap,position=position,leg_pos=leg_pos, $
		legend=legend,tim=tim,level_format=level_format,min_int=min_int, $
		no_ticks=no_ticks,charsize=charsize,normal=normal,rconv=rconv

	COMMON holly_prefs

; Allow several colour bars to be stacked
	IF N_PARAMS() NE 2 THEN BEGIN & ymaps=1 & ymap=0 & ENDIF

; Initialize colour bar position
    maxval=1000
    minval=0
	no_colours=10
	IF NOT KEYWORD_SET(no_ticks) THEN no_ticks=no_colours
	ysize=0.83
	xorigin=0.05
	yorigin=0.0
	IF (format AND 8) EQ 8 THEN BEGIN
		ysize=0.73
		yorigin=0.1
	ENDIF
	xpos=0.85+xorigin
	xbox=0.02/ymaps
	ypos=(ymaps-ymap-0.75)*ysize/ymaps+yorigin
	ybox_cols =0.6*ysize/(ymaps*no_colours)
	ybox_ticks=0.6*ysize/(ymaps*no_ticks)
	ybox_gnd  =0.6*ysize/(ymaps*10)
	IF keyword_set(charsize) THEN char=charsize ELSE char=MAX([1.0/ymaps,min_charsize])
	
	IF KEYWORD_SET(position) THEN BEGIN
		xpos= position(0)
		ypos= position(1)
		xbox= position(2)-position(0)
		ybox_cols =1.0*(position(3)-position(1))/no_colours
		ybox_ticks=1.0*(position(3)-position(1))/no_ticks
		ybox_gnd  =1.0*(position(3)-position(1))/10
		IF keyword_set(charsize) THEN char=charsize ELSE char=0.75*!P.CHARSIZE
	ENDIF

; Colours
	update_colour_info
	IF !D.NAME EQ 'PS' AND ps_font EQ 1 THEN BEGIN
		!P.FONT=0
		degree='!9'+STRING("260B)+'!3'
	ENDIF ELSE BEGIN
		!P.FONT=-1
		degree='!9%!3'
	ENDELSE

	rcin=FIX(FINDGEN(no_colours)*(ncol-4)/(no_colours-1))+1
	IF NOT KEYWORD_SET(rconv) THEN rcin=reverse(rcin)
	lvl=minval+FINDGEN(no_colours)*(maxval-minval)/no_colours

; Switch colour map for velocity plot, to maintain the red shift
; convention as well as the radar convention
	rcin=ROTATE(rcin,2)

; Draw coloured boxes
	FOR level=0,no_colours-1 DO					$
		POLYFILL,[xpos,xpos+xbox,xpos+xbox,xpos],		$
			 [ypos+ybox_cols*level,ypos+ybox_cols*level,		$
			  ypos+ybox_cols*(level+1),ypos+ybox_cols*(level+1)],	$
			  COLOR=rcin(level),NORMAL=NORMAL
			  
; If min_int specified, clear lower portion of bar
	IF KEYWORD_SET(min_int) THEN					$
		POLYFILL,[xpos,xpos+xbox,xpos+xbox,xpos],		$
			 [0,0,1,1]*ybox_cols*min_int/((maxval-minval)/no_colours)+ypos,		$
			  COLOR=background,NORMAL=NORMAL

; Draw outline
	FOR level=0,no_ticks-1 DO 					$
		PLOTS,xpos+xbox*[0,1,1,0,0],ypos+ybox_ticks*(level+[0,0,1,1,0]), $
			COLOR=foreground,NORMAL=NORMAL

; Plot levels
	IF NOT KEYWORD_SET(level_format) THEN BEGIN
		IF FIX((maxval-minval)/no_ticks) NE FLOAT((maxval-minval))/no_ticks THEN BEGIN
			level_format='(F10.1)'
		ENDIF ELSE BEGIN
			level_format='(I)'
		ENDELSE
	ENDIF
	lvl=minval+FINDGEN(no_ticks+1)*(maxval-minval)/no_ticks
	FOR level=1,no_ticks DO BEGIN
		numb=STRTRIM(FIX(ABS(lvl(level))),2)+'.'+STRTRIM(ABS(FIX(lvl(level)*10)) MOD 10,2)
		IF lvl(level) LT 0 THEN numb='-'+numb
		numb=STRTRIM(STRING(lvl(level),FORMAT=level_format),2)
		XYOUTS,xpos+1.4*xbox,ypos+ybox_ticks*level-0.25*ybox_ticks,		$
			numb,COLOR=foreground,CHARSIZE=0.8*char,NORMAL=NORMAL
	ENDFOR

; Plot title
	title='Velocity (ms!E-1!N)'
	IF KEYWORD_SET(legend) THEN title=legend
	IF NOT KEYWORD_SET(leg_pos) THEN leg_pos=1
	XYOUTS,xpos+0.075*leg_pos,ypos+no_colours*ybox_cols*0.5,title,COLOR=foreground,	$
		ORIENTATION=270,ALIGNMENT=0.5,CHARSIZE=char,NORMAL=NORMAL

; ok, ok, so the keyword TIM is a complete fudge so that Tim's pedant of a referee
; was satisfied... but see if I care...
	IF NOT KEYWORD_SET(tim) THEN tim=0
	line1=1.0+tim*0.5
	line2=1.5

	END

;-----------------------------------------------------------------------
PRO create_mp_header_file

;PRO create_mp_header_file
;-----------------------------------------------------------------------
;creates an html document containing all of the header info 
;from the potential_routines file.
;

openr,unit,getenv('MP_GO_PATH')+'potential_routines_new.pro',/get_lun
openw,unit2,getenv('MP_GO_PATH')+'idl_auto/top.temp',/get_lun
openw,unit3,getenv('MP_GO_PATH')+'idl_auto/bottom.temp',/get_lun
ch=''
printf,unit3,'<BR><BR>'
k=-1
WHILE NOT eof(unit) DO BEGIN
  readf,unit,format='(a100)',ch
  IF strmid(ch,0,2) EQ ';;' THEN BEGIN
    IF strmid(ch,2,3) EQ 'FUN' THEN BEGIN
      k=k+1
      IF k GT 0 THEN printf,unit3,'<A HREF="http://www.ion.le.ac.uk/~ag27/map_potential/idl.html">Back to top</A> <BR><BR><BR>'
      tit_len=strpos(ch,',')
      IF tit_len EQ -1 THEN tit_len=strpos(ch,' ',6)
      printf,unit2,'<A HREF="#'+strmid(ch,11,tit_len-11)+'">'+strmid(ch,11,tit_len-11)+'</A><BR>'
      printf,unit3,'<A NAME="'+strmid(ch,11,tit_len-11)+'"></A>'+strmid(ch,2,100)+'<BR>'
    ENDIF ELSE IF strmid(ch,2,3) EQ 'PRO' THEN BEGIN
      k=k+1
      IF k GT 0 THEN printf,unit3,'<A HREF="http://www.ion.le.ac.uk/~ag27/map_potential/idl.html">Back to top</A> <BR><BR><BR>'
      tit_len=strpos(ch,',')
      IF tit_len EQ -1 THEN tit_len=strpos(ch,' ',6)
      printf,unit2,'<A HREF="#'+strmid(ch,6,tit_len-6)+'">'+strmid(ch,6,tit_len-6)+'</A><BR>'
      printf,unit3,'<A NAME="'+strmid(ch,6,tit_len-6)+'"></A>'+strmid(ch,2,100)+'<BR>'
    ENDIF ELSE BEGIN
      printf,unit3,strmid(ch,2,100)+'<BR>'
    ENDELSE
  ENDIF
ENDWHILE
printf,unit3,'<A HREF="http://www.ion.le.ac.uk/~ag27/map_potential/idl.html">Back to top</A> <BR><BR><BR>'

close,unit
close,unit2
close,unit3
free_lun,unit
free_lun,unit2
free_lun,unit3
spawn,'cat /people/ag27/POTMAP_C/IDL/idl_auto/head.temp /people/ag27/POTMAP_C/IDL/idl_auto/top.temp /people/ag27/POTMAP_C/IDL/idl_auto/bottom.temp /people/ag27/POTMAP_C/IDL/idl_auto/foot.temp > /people/ag27/POTMAP_C/IDL/idl_auto/idl.html'

END

;-------------------------
pro test_nmodpnts

for L=1,7 do begin

    num_model=2*L
    FOR w=1,(L+1)/2 DO num_model=num_model+(4*(L-w))
    print,L,num_model
endfor
end
;-----------------------------------------------------------------

PRO reserve_map

COMMON map_pot
COMMON new_to_old

IF NOT KEYWORD_SET(map_open_flag) THEN BEGIN
  print,'Error: No map file loaded to reserve'
  return
ENDIF

map_data_alt=map_data
reserve_map_flag=1

END

;----------------------------------------------------------------------------------------
	FUNCTION READ_ACE_SWE_CDF,filename,path=path,append=append,quiet=quiet
	
	IF NOT KEYWORD_SET(path) THEN path=getenv('ACE_DATA_PATH')+'swe/h0/'
	open_file=path+'/'+filename+'.cdf'
	
; check on the file type	
	file_type=STRMID(filename,3,1)
	n_records=50000
;Set up the ace_swepam_info structure		
	;IF NOT KEYWORD_SET(append) THEN ace_swepam_info={spacecraft:'',year:-1,month:-1,day:-1,type:'',int:-1.0,rec:-1.}
	
; echo the filename to be opened to screen
    	IF NOT keyword_set(quiet) THEN BEGIN	
	  PRINT,'------------------------------------------'
	  IF KEYWORD_SET(append) THEN PRINT,'Appending file: ',STRTRIM(open_file,2) ELSE $
	  PRINT,'Opening file: ',STRTRIM(open_file,2)
	ENDIF
	IF file_type EQ 'k' THEN BEGIN
		PRINT,'**      Use level 2 (ho) data       **'
		return,0
	ENDIF
		
; Copy the old data to somewhere safe if the append keyword is set
	IF KEYWORD_SET(APPEND) THEN BEGIN
		old_ace_swepam=append
		;old_ace_swepam_info=ace_swepam_info
	ENDIF	

; set up the arrays to read data to
	ace_swepam_time=DBLARR(100000)
	ace_gse_pos=DBLARR(3,100000)
	ace_gse_proton_vel=DBLARR(4,100000)
	ace_gsm_pos=DBLARR(3,100000)
	ace_gsm_proton_vel=DBLARR(4,100000)
	ace_proton_density=DBLARR(100000)
	ace_proton_Pdyn=DBLARR(100000)
	
; Now do the reading
	swepam_cdf=CDF_OPEN(open_file)
	!QUIET=1

; First get the timing information	
		CDF_VARGET,swepam_cdf,'Epoch',epoch,rec_count=n_records,/ZVARIABLE
; Proton density	
		CDF_VARGET,swepam_cdf,'Np',np,rec_count=n_records,/ZVARIABLE		
; Proton speed	
		CDF_VARGET,swepam_cdf,'Vp',vp,rec_count=n_records,/ZVARIABLE			
; GSE positions
		CDF_VARGET,swepam_cdf,'SC_pos_GSE',gse_pos,rec_count=n_records,/ZVARIABLE	
; GSE bulk velocity
		CDF_VARGET,swepam_cdf,'V_GSE',vgse,rec_count=n_records,/ZVARIABLE	
; GSM positions 	
		CDF_VARGET,swepam_cdf,'SC_pos_GSM',gsm_pos,rec_count=n_records,/ZVARIABLE	
; GSM bulk velocity
  		CDF_VARGET,swepam_cdf,'V_GSM',vgsm,rec_count=n_records,/ZVARIABLE
	
	!QUIET=0
	CDF_CLOSE,swepam_cdf
	
	
; Now that the data has been read in from the files, fill the appropriate arrays
	ok_data=WHERE(epoch(1,*) GT 0.0)
	month=0
	    FOR i=0,N_ELEMENTS(ok_data)-1 DO BEGIN
		CDF_EPOCH,epoch(0,ok_data(i)),yr,month,day,hr,min,sec,/break			
		yr=DOUBLE(yr)
		IF (i EQ 0) THEN first_rec_year=yr		
		dy=day_no(yr,month,day)
		dy=DOUBLE(dy)
		hr=DOUBLE(hr)
		min=DOUBLE(min)
		sec=DOUBLE(sec)
; Let's convert this time into a number of seconds since the beginning of the year
; because CDF epoch values are a bit cumbersome
; Also need a check to see if number of data points is greater than the arraysize
		IF (i GT N_ELEMENTS(ace_swepam_time)-1) THEN BEGIN
			  PRINT,'ERROR: Data file is larger than array.'
			  PRINT,'Try resizing the arrays using the REMAKE command.'
			  PRINT,'---> EXITING'
			  RETURN,0
		ENDIF
		
		ace_swepam_time(i)=((dy-1)*86400.0)+(hr*3600.0)+(min*60.0)+sec
		ace_gse_pos(0,i)=gse_pos(0,ok_data(i))
		ace_gse_pos(1,i)=gse_pos(1,ok_data(i))
		ace_gse_pos(2,i)=gse_pos(2,ok_data(i))
		ace_gse_proton_vel(0,i)=VGSE(0,ok_data(i))
		ace_gse_proton_vel(1,i)=VGSE(1,ok_data(i))
		ace_gse_proton_vel(2,i)=VGSE(2,ok_data(i))
		ace_gse_proton_vel(3,i)=vp(0,ok_data(i))
		ace_proton_density(i)=np(0,ok_data(i))
		ace_gsm_pos(0,i)=gsm_pos(0,ok_data(i))
		ace_gsm_pos(1,i)=gsm_pos(1,ok_data(i))
		ace_gsm_pos(2,i)=gsm_pos(2,ok_data(i))
		ace_gsm_proton_vel(0,i)=VGSM(0,ok_data(i))
		ace_gsm_proton_vel(1,i)=VGSM(1,ok_data(i))
		ace_gsm_proton_vel(2,i)=VGSM(2,ok_data(i))
		ace_gsm_proton_vel(3,i)=vp(0,ok_data(i))
; Need to calculate dynamic pressure in nPa p=density*v^2=number density*mass of proton*v^2	
; Actually use a proton mass that assumes 5% of "protons" are actually aplha particles		
			
		IF ace_proton_density(i) GE 0.0 AND ace_gse_proton_vel(3,i) GE 0.0 THEN $
		ace_proton_Pdyn(i)=((ace_proton_density(i)*1E6)*(1.92E-27)*(ace_gse_proton_vel(3,i)*1000.0)^2)/1E-9 ELSE $
		ace_proton_Pdyn(i)=-1E31
	    ENDFOR
; clear the remaining elements of the array (if any exist). This prevents confusion
; if the previous file contained more records than the current one
	IF (N_ELEMENTS(ace_swepam_time) GT i) THEN BEGIN
		ace_swepam_time(i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_gse_pos(*,i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_gse_proton_vel(*,i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_gsm_pos(*,i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_gsm_proton_vel(*,i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_proton_density(i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
		ace_proton_Pdyn(i:N_ELEMENTS(ace_swepam_time)-1)=-1.0e+31
	ENDIF
	
;  Store some information about the file into the fgm_info structure
;	ace_swepam_info.spacecraft='ACE'
;	status = cnv_sec_mdhms(first_rec_year, month, day, hour, minute, sec, ace_swepam_time(0))
;	ace_swepam_info.year=first_rec_year
;	ace_swepam_info.month=month
;	ace_swepam_info.day=day
;	ace_swepam_info.type=STRMID(filename,3,2)
; create the final data array
	array_size=i
;	ace_swepam_info.rec=array_size	
;	IF STRMID(filename,3,2) EQ 'k1' THEN ace_swepam_info.int=3600.0
;	IF STRMID(filename,3,2) EQ 'k0' THEN ace_swepam_info.int=300.0
;	IF STRMID(filename,3,2) EQ 'h0' THEN ace_swepam_info.int=64.0
	ace_swepam=DBLARR(17,array_size)
	ace_swepam(*,*)=-1.0e+31
	
; At what number of seconds do the records occur i.e. 0,10,20,30 or 5,15,25?
; Find out the seconds of the first record and assume the timing remains
; constant throughout the file. Not perfect but good enough for now.
	start_sec=ace_swepam_time(0)
	
; Fill the array
	FOR j=0.0,array_size-1 DO BEGIN
		ace_swepam(0,j)=start_sec+64.0*j
	ENDFOR

		
; Transfer data from the cluster_time,cluster_gse_pos and cluster_gse_fgm arrays into the c array
	FOR k=0.0,i-1 DO BEGIN
; Calculate the index of the current cluster_time record in the new c array
		loc=(ace_swepam_time(k) MOD start_sec)/64.0
		ace_swepam(1,loc)=ace_gse_pos(0,k)
		ace_swepam(2,loc)=ace_gse_pos(1,k)
		ace_swepam(3,loc)=ace_gse_pos(2,k)
		ace_swepam(4:7,loc)=ace_gse_proton_vel(*,k)
		ace_swepam(8:10,loc)=ace_gsm_pos(*,k)
		ace_swepam(11:14,loc)=ace_gsm_proton_vel(*,k)
		ace_swepam(15,loc)=ace_proton_density(k)
		ace_swepam(16,loc)=ace_proton_Pdyn(k)
	ENDFOR	
	

	IF KEYWORD_SET(APPEND) THEN BEGIN	
	    	result=where(append(0,*) GE 0.0,count)
		array_size=count+array_size
		new_ace=DBLARR(17,array_size)
		new_ace(*,0:count-1)=old_ace_swepam
		new_ace(*,count:array_size-1)=ace_swepam
		ace_swepam=new_ace
		;ace_swepam_info=old_ace_swepam_info
		;ace_swepam_info.rec=array_size	
	ENDIF	
	RETURN,ace_swepam
	END



;------------------------------------------------------------------------------	

FUNCTION get_dst_mp,year,month,day,hour
	
	dst_path='/data/ag27/dst/'
	open_file=dst_path+'dst'+STRING(year,FORMAT='(I4.4)')+'.txt'

	OPENR, inunit, open_file, /GET_LUN

; initialise a couple of variables
	line = ''

	i=0
; Now begin reading in in lines from the file until EOF is reached	
	WHILE NOT EOF(inunit) DO BEGIN
		READF, inunit, line
		hourly_values=INTARR(24)
		READS,line,yr,mn,dy,hourly_values,FORMAT='(3x,I2,I2,1x,I2,10x,24I4)'
		IF mn EQ month AND dy EQ day THEN dst=hourly_values(hour)
	ENDWHILE
	FREE_LUN, inunit
	
	RETURN,dst	
		
END

;------------------------------------------------------------------------------

FUNCTION get_hmb,hmb,yr,mo,dy,shr,smin

 latref=59
 nlat=26
 latmin=1.0*hmb
 bfac=(90-latmin)/(90-latref)
 del_L=bfac*5.5
 latx1=fltarr(2,73)

 FOR i=0,72 DO BEGIN
   lon=i*5.
   mlt1=MLT(yr,cnvtime(yr,mo,dy,shr,smin,0),lon)
   latx=latmin
   IF mlt1 GE 11 AND mlt1 LE 19 THEN latx=latmin+del_L*(1+cos((!PI/8)*(mlt1-11))) $
   ELSE IF mlt1 LT 11 AND mlt1 GE 5 THEN latx=latmin+del_L*(1+cos((!PI/6)*(11-mlt1)))
   latx1(*,i)=[lon,latx]
 ENDFOR

RETURN,latx1
END

;----------------------------------------------------------------------------------

FUNCTION get_current_latmins

COMMON map_pot
COMMON new_to_old

hmbnd={year:0, month:0, day:0, hour:0, min:0, sec:0, lat:0 }
hmbnd=replicate(hmbnd,n_elements(map_data))
for i=0,n_elements(map_data) -1 DO hmbnd(i)={year:map_data(i).syr, month:map_data(i).smo, day:map_data(i).sdy, $
   hour:map_data(i).shr, min:map_data(i).smt, sec:map_data(i).ssc, lat: fix(map_data(i).latmin) }
   
RETURN,hmbnd
END

;----------------------------------------------------------------------------------

FUNCTION constrain_latmin,initial

;given a set of latmins, this function returns a new set in which the latitude change between successive 
;maps is limited by the value of Vpc

COMMON map_pot
COMMON new_to_old

IF n_params() EQ 0 THEN initial=0

hmbnd={year:0, month:0, day:0, hour:0, min:0, sec:0, lat:0.0 }
hmbnd=replicate(hmbnd,n_elements(map_data))
hmbnd(initial)={year:map_data(initial).syr, month:map_data(initial).smo, day:map_data(initial).sdy, $
                hour:map_data(initial).shr,min:map_data(initial).smt, sec:map_data(initial).ssc, $
                lat: map_data(initial).latmin }
            
Re = 6375.0d*1000.
Altitude = 300.0*1000.0
bpolar = 6.2e-5
done=0
n=[initial,n_elements(map_data)]
jj=[0,initial]
direc=[-1,1]
m=[1,2]
FOR d=0,1 DO BEGIN
 IF initial EQ 0 AND d EQ 0 THEN CONTINUE
 FOR h=jj(d),n(d) - m(d) DO BEGIN
  IF d EQ 0 THEN j=initial-h ELSE j=h
  k=j+direc(d)
  r1=re*cos(hmbnd(j).lat*!dtor)
  a1=!pi*r1^2
  vpc=map_data(j).pot_drop
  phi=90-hmbnd(j).lat
  Bi=bpolar*(1.0-3.0*Altitude/Re)*sqrt(3.0*cos(phi*!dtor)^2+1.0)/2.0
  i=(60*map_data(j).emt+map_data(j).esc)-(60*map_data(j).smt+map_data(j).ssc)
  IF i LT 0 THEN i=i+3600
  a2=[a1+i*vpc/Bi,a1-i*vpc/Bi]
  r2=sqrt(a2/!pi)
  latmin_range=acos(r2/re)/!dtor
  latmin_range=[hmbnd(j).lat-0.5,hmbnd(j).lat+0.5]
  IF map_data(k).latmin LT latmin_range(0) THEN new_latmin=latmin_range(0) ELSE $
  IF map_data(k).latmin GT latmin_range(1) THEN new_latmin=latmin_range(1) ELSE $
                                                new_latmin=map_data(k).latmin
  hmbnd(k)={year:map_data(k).syr, month:map_data(k).smo, day:map_data(k).sdy, hour:map_data(k).shr, $
            min:map_data(k).smt, sec:map_data(k).ssc, lat: float(new_latmin) }
  
  ;print,d,h,j,k,map_data(j).latmin,map_data(k).latmin,new_latmin
 ENDFOR
ENDFOR

RETURN,hmbnd
END

;----------------------------------------------------------------------------------
