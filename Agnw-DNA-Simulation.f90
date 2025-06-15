implicit double precision(a-h,o-z)
 character(80)  :: filename = ' '
integer        :: i = 0
 REAL :: r(1)

real:: xpol(1000), ypol(1000), zpol(1000)		!Npol

real:: xw(100), yw(100), zw(100)


open(19,file='input.txt',status='unknown')  !results
open(20,file='output.txt',status='unknown')  !results


read(19,*) ndrawscreen
read(19,*) ndrawsystem



!constants and values
pi=3.141592
epsi0=8.85e-12
epsir=82.
xkB=1.38e-23
temp=273.15+15.
xNA=6.022e23
e=1.6e-19

nwire=10    !number of points for Agnw sampling


read(19,*) nwirea
read(19,*) conc
read(19,*) xL
read(19,*) xrAgbulk
read(19,*) xNmon
read(19,*) xdmon
read(19,*) xMwmon
read(19,*) xLp
read(19,*) xV
read(19,*) sig



!constantes del potencial Likos
epsipot=1.87
sigmapot=1.08




!__________________________________________________________________________________
!CALCULATION OF MAXIMUM APPROACH DISTANCE
read(19,*) xmM
rhos=xmM*xNA
read(19,*) zetapot
zetapot=zetapot/1000.
read(19,*) xchargemon


xRc=xrAgbulk*1e-9    !Agnw radius
xNbp=2*xLp/xdmon     !base pairs in a Kuhn segment
Zbp=xchargemon*xNmon  !charge of a polymer chain


!Characteristic lengths
xlb=e**2/(4*pi*epsi0*epsir*xkB*temp)    !Bjerrum length
xkappa=(8*pi*xlb*rhos)**0.5	 	!Inverse of the Debye length



!limits
r0=xkappa*xRc
rf=50*xkappa*xRc
itot=100000
Dr=(rf-r0)/real(itot)
xInt=0

!K0ini(xkappa*xRc)
arg=r0
CALL CALCK0(ARG,RESULT,JINT)
BESEK0ini = RESULT

A=zetapot*Zbp*e/(BESEK0ini*xkB*temp)

JINT = 1   !pide la funcion de Bessel

do i=1,itot
arg=xkappa*xRc+i*Dr
xr=arg/xkappa
CALL CALCK0(ARG,RESULT,JINT)
BESEK0 = RESULT
xInt=xInt+(1-exp(-A*besek0))*(xr*1e9)*(Dr*1e9/xkappa)    !nm2
enddo


xReff=((xRc*1e9)**2+2*xInt)**0.5

if(ndrawscreen.ne.0)then
write(6,*) '______________________________________'
write(6,*) 'Maximum approach distance calculation'
write(6,*) 'xReff=', xReff, 'nm', '         xReff-xRc=', xReff-xRc*1e9, 'nm'
write(6,*) '______________________________________'
write(6,*) '           '
endif

xrAg=xReff



!__________________________________________________________________________________
!MAIN QUANTITIES IN CALCULATIONS

if(ndrawscreen.ne.0)then
write(6,*) 'Input data'
write(6,*) 'xrAg', xrAg, 'nm'
write(6,*) 'xL', xL, 'nm'
write(6,*) 'sigma', sig
write(6,*) 'xNmon', xNmon
write(6,*) 'xLp', xLp, 'nm'
write(6,*) 'xV', xV, 'nm'
write(6,*) 'zetapot', zetapot*1000, 'mV'
write(6,*) 'xmM', xmM, 'mM'
write(6,*) 'monomer charge', xchargemon, 'e units'
write(6,*) '______________________________________'
write(6,*) '       '
endif


xLc=xNmon*xdmon   !Longitud de contorno en nm
xLk=2*xLp     !Longitud de Kuhn en nm
xMw=xNmon*xMwmon   !Peso molecular
Rg=xLp*sqrt(xLc/(3*xLp)-1+2*xLp/xLc-2*(xLp/xLc)**2*(1-exp(-xLc/xLp))) !en nm

Nk=xLc/(2*xLp)




!_________________________________------
!CICLE ON EXPERIMENTS
do nexp=12,16  !cicle on experiments

conc=10*nexp   !en mg/L
npol=conc/xMw*xNa*1e-18*(xV/1000.)**3


!__________________________________________________________________________________
!ANISOTROPIC POLYMER BUILD

elon=0.8

apol=Rg+256*(1-exp(-(elon/1.22)**0.81))  !Elongated radius (um)
bpol=Rg-275*(1-exp(-(elon/2.64)**1.))


!___________________________________________________
!POLYMER PLACING
!Colocamos el polímero con solapamiento restringido


!colocamos la primera madeja
xpol(1)=xV*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
ypol(1)=xV*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
zpol(1)=xV*(r(1)-0.5)




!Colocamos la madeja n, que solapa con las anteriores según una probabilidad xP

do n=2,npol  !ciclo polímeros

nciclos=0

do while (ncondition.lt.(n-1))   !Ciclo condition

ncondition=0


CALL RANDOM_NUMBER(r)   
xpol(n)=xV*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
ypol(n)=xV*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
zpol(n)=xV*(r(1)-0.5)

do j=1,n-1    !ciclo en polímeros anteriores



xpar=((zpol(n)-zpol(j))**2/(sigmapot*apol)**2+(xpol(n)-xpol(j))**2/(sigmapot*bpol)**2+(ypol(n)-ypol(j))**2/(sigmapot*bpol)**2)
xUnorm=epsipot*(exp(-xpar))*5
xProb=exp(-xUnorm)


CALL RANDOM_NUMBER(r)
xRan=r(1)


if(xRan.le.xProb)then   !AQUI
!if(xUnorm.le.xRan)then
ncondition=ncondition+1
else
nciclos=nciclos+1
endif


enddo   !ciclo en polímeros anteriores
enddo !ciclo condition 
enddo!ciclo polímeros

!_________________________________________________________________
!DRAW ELLIPSOIDS


if(ndrawsystem.ne.0)then !condición dibujo elipsoides

write(filename,'(a,i0,a)')'mesh-ellip-',nexp,'.txt'
open(16,file=trim(filename),status='unknown') !mesh ellipsoids


nite2=20

do i=1,npol

!En cada punto del mesh un elipsoide
do n=1, nite2
do m=1, nite2

u=real(m-1)/real(nite2-1)*2*pi-pi
v=real(n-1)/real(nite2-1)*pi-pi/2


x=bpol*sin(u)*cos(v)+xpol(i)
y=bpol*sin(u)*sin(v)+ypol(i)
z=apol*cos(u)+zpol(i)


write(16,*) x,y,z

u=real(m+1-1)/real(nite2-1)*2*pi-pi    !writing twice for better image
v=real(n+1-1)/real(nite2-1)*pi-pi/2


x=bpol*sin(u)*cos(v)+xpol(i)
y=bpol*sin(u)*sin(v)+ypol(i)
z=apol*cos(u)+zpol(i)

write(16,*) x,y,z

enddo

write(16,*) ' '  !Este espacio para gnuplot
write(16,*) ' '  !Este espacio para gnuplot
write(16,*) ' '  !Este espacio para gnuplot

enddo

write(16,*) ' '  !Este espacio para gnuplot
write(16,*) ' '  !Este espacio para gnuplot
write(16,*) ' '  !Este espacio para gnuplot

enddo




endif   !condición dibujo elipsoides







!__________________________________________________________________________________
!SILVER NANOWIRE POSITIONS

P2Ag=0
naccepted=0
nrejected=0

nlimit=0

do while (naccepted.lt.nwirea)   !Cicle on Agnws

noverlap=0

CALL RANDOM_NUMBER(r)   !Initial position for Agnw CM
x0=(xV-xL/2)*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
y0=(xV-xL/2)*(r(1)-0.5)
CALL RANDOM_NUMBER(r) 
z0=(xV-xL/2)*(r(1)-0.5)
CALL RANDOM_NUMBER(r)  !position jump (angle1) in cosine
q1=acos(2*(r(1)-0.5))
CALL RANDOM_NUMBER(r) !position jump (angle2)
q2=2*pi*r(1)



!Obstacles
noverlap=0
do n=1,nwire  !Sampling the wire length
xw(n)=x0-(xL/2.)*sin(q1)*cos(q2)+real(n-1)*xL/real(nwire-1)*sin(q1)*cos(q2)
yw(n)=y0-(xL/2.)*sin(q1)*sin(q2)+real(n-1)*xL/real(nwire-1)*sin(q1)*sin(q2)
zw(n)=z0-(xL/2.)*cos(q1)+real(n-1)*xL/real(nwire-1)*cos(q1)
do i=1,npol
if(((zw(n)-zpol(i))**2)/apol**2+((xw(n)-xpol(i))**2)/bpol**2+((yw(n)-ypol(i))**2)/bpol**2.le.1)then     !Ellipsoid field direction -> z
noverlap=noverlap+1
endif
enddo
enddo


!________________________________________
!IF CICLE FOR NO OVERLAP

if(ndrawsystem.ne.0)then
write(filename,'(a,i0,a)')'Agnws-accepted-',nexp,'.txt'
open(12,file=trim(filename),status='unknown') !AgNW positions and orientations
write(filename,'(a,i0,a)')'Agnws-rejected-',nexp,'.txt'
open(13,file=trim(filename),status='unknown') !AgNW positions and orientations
endif

if(noverlap.eq.0)then
P2Ag=P2Ag+0.5*(3*cos(q1)**2-1)
naccepted=naccepted+1
if(ndrawsystem.ne.0)then
xf=(xL/2.)*sin(q1)*cos(q2)
yf=(xL/2.)*sin(q1)*sin(q2)
zf=(xL/2.)*cos(q1)
write(12,*)  (x0-xf),(y0-yf),(z0-zf), (2*xf),(2*yf), (2*zf)
endif
endif


!________________________________________
!IF CICLE FOR OVERLAP
if(noverlap.ne.0)then
nrejected=nrejected+1
if(ndrawsystem.ne.0)then
xf=(xL/2.)*sin(q1)*cos(q2)
yf=(xL/2.)*sin(q1)*sin(q2)
zf=(xL/2.)*cos(q1)
write(13,*)  (x0-xf),(y0-yf),(z0-zf), (2*xf),(2*yf), (2*zf)
endif
endif

nlimit=nlimit+1

if(nlimit.gt.nwirea*1000)then
write(6,*) 'OVERTIME'
go to 1
endif

enddo  !Cicle on Agnw positions

1 continue


P2Ag=P2Ag/real(naccepted)
write(20,*) conc, P2Ag

if(ndrawscreen.ne.0)then
write(6,*) '    nciclos=  ', nciclos
write(6,*) '    accepted        ', 'rejected      ', 'percentage   '
write(6,*) naccepted, nrejected, naccepted/real(naccepted+nrejected)*100.
write(6,*) '    '
write(6,*) '    conc        ', 'P2Ag      ', 'npol'
write(6,*) conc, P2Ag, npol
endif

enddo   !cicle of experiments

 close(19)
 close(20)



end program





!------------------------------------------------------------------


SUBROUTINE CALCK0(ARG,RESULT,JINT)
  implicit none

  integer ( kind = 4 ) i,jint
  real ( kind = 8 ) &
      arg,f,g,one,p,pp,q,qq,result,sumf,sumg,sump,sumq,temp, &
      x,xinf,xmax,xsmall,xx,zero
  dimension p(6),q(2),pp(10),qq(10),f(4),g(3)
!
!  Mathematical constants
!
  data one/1.0d0/,zero/0.0d0/
!
!  Machine-dependent constants
!
  data xsmall/1.11d-16/,xinf/1.79d+308/,xmax/705.342d0/
!
!  Coefficients for XSMALL <=  ARG  <= 1.0
!
  data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
        1.1999463724910714109d+01, 4.6850901201934832188d+02, &
        5.9169059852270512312d+03, 2.4708152720399552679d+03/
  data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
  data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
       -1.7733784684952985886d+04,-4.0320340761145482298d+05/
  data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
       -1.6128136304458193998d+06/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
        3.1075408980684392399d+04, 1.0577068948034021957d+05, &
        1.7398867902565686251d+05, 1.5097646353289914539d+05, &
        7.1557062783764037541d+04, 1.8321525870183537725d+04, &
        2.3444738764199315021d+03, 1.1600249425076035558d+02/
  data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
        3.1474655750295278825d+04, 9.7418829762268075784d+04, &
        1.5144644673520157801d+05, 1.2689839587977598727d+05, &
        5.8824616785857027752d+04, 1.4847228371802360957d+04, &
        1.8821890840982713696d+03, 9.2556599177304839811d+01/

  x = arg

  if (x > zero) then

    if (x <= one) then
!
!  0.0 <  ARG  <= 1.0
!
      temp = log(x)

      if (x < xsmall) then
!
!  Return for small ARG
!
        result = p(6)/q(2) - temp

      else

        xx = x * x
        sump = ((((p(1)*xx + p(2))*xx + p(3))*xx + &
          p(4))*xx + p(5))*xx + p(6)
        sumq = (xx + q(1))*xx + q(2)
        sumf = ((f(1)*xx + f(2))*xx + f(3))*xx + f(4)
        sumg = ((xx + g(1))*xx + g(2))*xx + g(3)
        result = sump/sumq - xx*sumf*temp/sumg - temp
        if (jint == 2) result = result * exp(x)

      end if

    else if ((jint == 1) .and. (x > xmax)) then
!
!  Error return for ARG > XMAX
!
      result = zero

    else
!
!  1.0 < ARG
!
      xx = one / x
      sump = pp(1)
      do i = 2, 10
        sump = sump*xx + pp(i)
      end do
      sumq = xx
      do i = 1, 9
        sumq = (sumq + qq(i))*xx
      end do
      sumq = sumq + qq(10)
      result = sump / sumq / sqrt(x)
      if (jint == 1) result = result * exp(-x)

    end if

  else
!
!  Error return for ARG <= 0.0
!
    result = xinf

  end if

  return
end
