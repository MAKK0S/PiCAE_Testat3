c      Titel: Isotrope, viskoelastische UMAT im Zeitbereich
c      zu beachten: gleiche Relaxationszeiten für die
c                   unterschiedlichen Materialkennwerte müssen
c                   verwendet werden
c                   Es werden 6+6 Zustandsvariablen benötigt
c
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,knc)
c------Deklaration ABAQUS
      implicit none
      integer :: kspt,layer,npt,noel,nprops,nstatv,ntens,
     1 nshr,ndi,knc
      double precision :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,
     1 pnewdt,celent,dfgrd0(3,3),dfgrd1(3,3),time(2),stress(ntens),
     2 statev(nstatv),ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),predef(1),dpred(1),props(nprops),
     4 coords(3),drot(3,3)
      integer :: jstep(4)
      character*80 cmname
c------Lokale Deklarationen
      integer :: i, j
      real :: stressm(ntens),stressE(ntens),
     1 zero,one,two,three,four,half,ten,expTerm
      real :: G0,
     1 k0,G,k,Gm,km,f0(ntens,ntens),
     2 f(ntens,ntens)
c
c-----Nuetzliche Zahlen
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0,
     1  half=0.5d0,ten=10.d0)
c
c-----Initialisieren mit null
c     gi=zero; ki=zero; taui=zero
C-----Anzahl a der Parameter bestimmen
      integer a=1
      real b=props(a)
      write (*,*) b,props(1000)
      do while (b/=zero)
          b=props(a)
          props(a)
          a=a+1
      enddo 

c-----Überprüfen ob Rest a/5 zero
      
c-----Deklaration für Prony-Reihe
      real E(a/5),nu(a/5),gi(a/5),taui(a/5),ki(a/5)
      
c-----Materialparameter aus props einlesen
      do i=0,(a/5)
          E(i)=props(i+1); nu(i)=props(i+2)
          gi(i)=props(i+3)
          ki(i)=props(i+4)
          taui(i)=props(i+5)
      end do
c
c-----viskoelastischer Schritt
c-----Zuweisen der Materialparamter
c-----Schleife
      G0=E(0)/(two*(one+nu(0))); k0=E(0)/(three*(one-two*nu(0)))

      G=G0-G0*gi
      k=k0-k0*ki
c-----elastische Steifigkeitsmatrix
c-----Schleife?
      f0=zero
      do i=1,3
        do j=1,3
          f0(i,j)=k-two/three*G
        end do
        f0(i,i)=k+four/three*G
      end do
      do i=4,6
        f0(i,i)=G
      end do
      km=zero;Gm=zero;f=zero
      Gm=gi*G0
      km=ki*k0
c-----Komponenten der Steifigkeitsmatrix für jedes Prony-Element
      do i=1,3
        do j=1,3
          f(i,j)=km-two/three*Gm
        end do
        f(i,i)=km+four/three*Gm
      end do
      do i=4,6
        f(i,i)=Gm
      end do
c
c-----Berechnung des elastischen Spannungsanteils
c-----Schleife (bis a/5)
      stressE=zero
      do i=1,3
        stressE(i)=statev(i)
        do j=1,3
          stressE(i)=stressE(i)+f0(i,j)*dstran(j)
        end do
      end do
      do j=4,ntens
        stressE(j)=statev(j)+f0(j,j)*dstran(j)
      end do
c-----Relaxieren der Spannung
      stressm=zero
      do j=1,ntens
        stressm(j)=exp(-dtime/taui)*statev(6+j)
      end do
c-----Berechnung der Spannung am Ende des Inkrementes
      expTerm=taui*(one-exp(-dtime/taui))
c
      do i=1,3
        do j=1,3
          stressm(i)=stressm(i)+f(i,j)*dstran(j)/dtime*expTerm
        end do
      end do
      do j=4,ntens
        stressm(j)=stressm(j)+f(j,j)*dstran(j)/dtime*expTerm
      end do
c-----Berechnung der Gesamtspannung am Ende des Inkrementes
c-----schleife (schauen wie lange schleife laufen muss)
      stress=zero
      do i=1,6
        stress(i)=stressE(i)+stressm(i)
      end do
c
c-----Jacobi-Matrix
      ddsdde=zero
      ddsdde(1,1)=(k+four/three*G)
      ddsdde(1,2)=(k-two/three*G)
      ddsdde(4,4)=G

      ddsdde(1,1)=ddsdde(1,1)+(km+four/three*Gm)*
     1  taui/dtime*(one-exp(-(dtime/taui)))
      ddsdde(1,2)=ddsdde(1,2)+(km-two/three*Gm)*taui/dtime*
     1  (one-exp(-(dtime/taui)))
      ddsdde(4,4)=ddsdde(4,4)+Gm*taui/dtime*
     1  (one-exp(-(dtime/taui)))

      ddsdde(2,2)=ddsdde(1,1)
      ddsdde(3,3)=ddsdde(1,1)
      ddsdde(1,3)=ddsdde(1,2)
      ddsdde(2,1)=ddsdde(1,2)
      ddsdde(2,3)=ddsdde(1,2)
      ddsdde(3,1)=ddsdde(1,2)
      ddsdde(3,2)=ddsdde(1,2)
      ddsdde(5,5)=ddsdde(4,4)
      ddsdde(6,6)=ddsdde(4,4)
c
c-----Aktualisieren der gespeicherten Zustandsvariablen
      statev(1:6)=stressE(1:6)
      do i=1,ntens
        statev(6+i)=stressm(i)
      end do
c
      end subroutine umat
