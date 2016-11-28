
PRO shockpars_ip,v_shock=v_shock, t1=t1, b1=b1, usw=usw, rho1=rho1
 
 ; Recoded 08/07/2016
 ; shock is always x-planar, propagating radially

 ; This scripts assumes tangential components are in z-direction only. 

 ; Shock propagation speed is an input parameter (in SNIF) - this gives Vx
 ; solar wind speed Usw is also given in SNIF, and is radial. Thus,
 ; we simply deduct this from V_shock to find the shock-normal inflow speed.
 ; Vz is solved via velocity transform to dHT
 ;
 ; Shock-normal angle is calculated via solar wind speed as per Parker
 ; spiral winding, at 1 AU.

 r_s = 6.96e8
 AU = 215.*r_s
 c  = 2.99792458e8
 m_p = 1.6726e-27
 e_p = 1.6021773e-19
 k_b = 1.3807e-23
 adiab = 5./3.
 mu0 = 1.25664e-6
 avogadro = 6.0221367e23
 Rydb = k_b*avogadro
; OmegaSun = 2.972109871e-6 ;in radians, based on equatorial
; differential rotation speed
 OmegaSun = 2.6394e-6 ;in radians, based on 27 days / rotation

 IF ~keyword_set(v_shock) THEN v_shock=7.50e5 ;default 750 km/s
 IF ~keyword_set(t1) THEN t1=1.0e5 ;default 0.1 MK 
 IF ~keyword_set(b1) THEN b1=5.e-9 ;default 5 nT
 IF ~keyword_set(usw) THEN usw=5.00e5 ;default 500 km/s
 IF ~keyword_set(rho1) THEN rho1=5.0e6 ;default 5 particles cm^-3

 thetrad = atan(OmegaSun*(AU-r_s)/usw)
 theta = (180./!pi)*thetrad
 cos11 = cos(thetrad)
 cos12 = cos11^2
 sin11 = sin(thetrad)
 sin12 = sin11^2
 ; angle, sine and cosine are non-negative

 ;Solve DHT frame values
 bn1 = b1*cos11
 bt1 = b1*sin11
 ; field is pointing outwards (positive in both x and z)
 un1 = -v_shock+usw
 ut1 = un1*bt1/bn1
 uht1 = -sqrt(un1^2+ut1^2)
 ; now B and U are antiparallel

 print, 'Angle ',theta
 print, 'dHT transformation speed ',ut1
 
 va1 = sqrt(b1*b1/(mu0*rho1*m_p))
 ;Use: pressure = n k T (does not include electron pressure)
 pressure1 = rho1 * k_B * t1
 vsound1 = sqrt( adiab*pressure1 /(m_p *rho1) )
 
 u12 = uht1^2
 va12 = va1^2
 vs12 = vsound1^2
 MA2 = u12/va12

 print,' va1 ',va1,' vsound1 ',vsound1, ' uht1 ',uht1,' MA ',sqrt(MA2)
; print,' un1 ',un1,' ut1 ',ut1, ' sqrt(un1^2+ut1^2) ',sqrt(un1^2+ut1^2)

 ; Initialize compression ratio for looping
 compr=0

 ; Calculate beta
 beta1 = (2./adiab)* vs12/va12
 calctemp1 = 1.D + 0.5*adiab*beta1

 ; Calculate plasma compression ratio through iterative Newton method
 Ztry = ((0.5D/cos12)*(calctemp1 + sqrt(calctemp1^2 - 2.*adiab*beta1*cos12)) -1.) > $     ; First (root for M^2) -1
    ((0.5D/cos12)*(calctemp1 - sqrt(calctemp1^2 - 2.*adiab*beta1*cos12)) -1.) > 0.        ; Second and third (root for M^2) -1

 fform = (1.D +Ztry)*((Ztry^2)*8*cos12 +(3.D -5.D*Ztry)*sin12) -(Ztry^2)*5*beta1
 gform = (1.D +Ztry)*((Ztry^2)*2*cos12 +(3.D +Ztry)*sin12)
 Rtry = fform/gform

 M2try = (1.D +Ztry)*Rtry

 rstep = 1.0
 WHILE ((rstep Ge 0.0001) && (compr LT 0.001))  DO BEGIN
    Ztry = ((0.5D/cos12)*(calctemp1 + sqrt(calctemp1^2 - 2.*adiab*beta1*cos12)) -1.) > $        ; First (root for M^2) -1
       ((0.5D/cos12)*(calctemp1 - sqrt(calctemp1^2 - 2.*adiab*beta1*cos12)) -1.) > 0.           ; Second and third (root for M^2) -1
    
    fform = (1.D +Ztry)*((Ztry^2)*8*cos12 +(3.D -5.D*Ztry)*sin12) -(Ztry^2)*5*beta1
    gform = (1.D +Ztry)*((Ztry^2)*2*cos12 +(3.D +Ztry)*sin12)
    Rtry = fform/gform
    M2try = (1.D +Ztry)*Rtry
    
    WHILE abs(M2try - MA2) GT 0.0001 DO BEGIN
       fderi = (Ztry^2)*8.D*cos12 +(3.D -8.D*Ztry)*sin12 -10.D*Ztry*beta1 +(1.D +Ztry)*(16.D*Ztry*cos12 -5.D*sin12)
       gderi = (Ztry^2)*2.D*cos12 +(3.D +Ztry)*sin12 +(1.D +Ztry)*(4.D*Ztry*cos12 +sin12)
       rderi = (gform*fderi-fform*gderi)/(gform^2)
       m2deri = (1.D +Ztry)*rderi + Rtry
       
       ; Newton step forward
       Ztry = Ztry + (MA2 - M2try)/m2deri * 0.5*rstep
       
       ; Calculate new Rtry and M2try
       fform = (1.D +Ztry)*((Ztry^2)*8*cos12 +(3.D -5.D*Ztry)*sin12) -(Ztry^2)*5*beta1
       gform = (1.D +Ztry)*((Ztry^2)*2*cos12 +(3.D +Ztry)*sin12)
       Rtry = fform/gform
       M2try = (1.D +Ztry)*Rtry
       
    ENDWHILE
;    print, 'Ztry ',Ztry,' Rtry ',Rtry,' M2try ',M2try,' MA2 ',MA2
    IF (Rtry LE MA2) THEN compr = Rtry
    rstep = rstep * 0.1
 ENDWHILE

 ; Now solve magnetic compression ratio
 comprb = sqrt( cos12 + (1.-cos12)*( ( compr * (u12-va12)/(u12-va12*compr))^2 ) )

 un2 = un1/compr ;still negative
 bn2 = bn1 ;still positive
 bt2 = sqrt(comprb^2 -cos12)*b1 ;still positive

 theta2 = atan(bt2/bn2) ;positive
 cos21 = cos(theta2)
 sin21 = sin(theta2)
 cos22 = cos21^2
 sin22 = sin21^2

 uht2 = un2 / cos21
 ut2 = uht2 * sin21
 rho2 = rho1*compr

 pressure2 = pressure1*compr + pressure1*(adiab-1)*compr*(uht1*uht1-uht2*uht2)/(2.*vs12)

 ;Use: pressure = n k T (does not include electron pressure)
 t2 = pressure2/(rho2*k_b)

 ; Print values to be copy-pasted into IPShock.cfg
 print,'VX0u = ',un1
 print,'VY0u = ',0.0
 print,'VZ0u = ',ut1
 print,'BX0u = ',bn1
 print,'BY0u = ',0.0
 print,'BZ0u = ',bt1
 print,'rhou = ',rho1
 print,'Temperatureu = ',t1

 print,'VX0d = ',un2
 print,'VY0d = ',0.0
 print,'VZ0d = ',ut2
 print,'BX0d = ',bn2
 print,'BY0d = ',0.0
 print,'BZ0d = ',bt2
 print,'rhod = ',rho2
 print,'Temperatured = ',t2

 ; Print values to input into upstream.dat (and, if wanted, downstream.dat)
 print,''
 print,'upstream.dat'
 print, 0.0, rho1, t1, un1, 0.0, ut1, bn1, 0.0, bt1

 print,''
 print,'downstream.dat'
 print, 0.0, rho2, t2, un2, 0.0, ut2, bn2, 0.0, bt2

 print,''
 print,'Upstream particle parameters'
 gyrofreq = e_p*b1/m_p
 v_th = sqrt(2*t1*k_b/m_p)
 print,'  gyrofreq ',gyrofreq,' gyrotime ',2.*!pi/gyrofreq 
 print,'  v_th ',v_th,' ion inertial range ',v_th/gyrofreq
 print,'  bulk larmor radius ',abs(uht1)/gyrofreq,' ion Alfvenic range ',va1/gyrofreq

END
