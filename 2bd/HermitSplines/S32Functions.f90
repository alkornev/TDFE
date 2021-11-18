!>>>>>>>>> This module is not yet implemented!!!!!<<<<<<<<<<<<<<<<<<<<
module HermitSplines
!  This module takes care of the system discretization. 
!  It supports some basic functionality for S_32 Hermit splines and orthogonal collocations.
!
  implicit none
  integer :: Nx,Ncx
  integer :: BCLeft=0,BCRight=0
  double precision :: Xl,Xr
  double precision, allocatable :: XG(:),XGc(:)
  double precision, allocatable :: LeftBounds(:),RightBounds(:)
  logical :: ifAllocated=.false.
  contains
!
  subroutine InitDiscretization(NPoints,leftBound,rightBound,specialPoint, &
      &                         densityLeft,densityRight,iBCLeft,iBCRight)
    integer, intent(in) :: iBCLeft,iBCRight
    integer, intent(in) :: NPoints
    double precision :: leftBound,rightBound,specialPoint
    double precision :: densityLeft,densityRight
    Xl=leftBound
    Xr=rightBound
    BCLeft=iBCLeft
    BCRight=iBCRight
    NX=NPoints
    Ncx=2*Nx
    if (abs(BCRight).gt.1) then
      stop 'Error initializing the grid, BCRight can be only 0,1 or -1'
    end if
    if (abs(BCLeft).gt.1) then
      stop 'Error initializing the grid, BCleft can be only 0,1 or -1'
    end if
    if (BCRight.eq.-1) Ncx=Ncx+1
    if (BCLeft.eq.-1) Ncx=Ncx+1
    if (.not.ifAllocated) then
      allocate(Xg(0:Nx),XGc(Ncx),LeftBounds(Ncx),RightBounds(Ncx))
      ifAllocated=.true.
    else if (size(XGc).ne.Ncx) then
        deallocate(Xg,XGc,LeftBounds,RightBounds)
        allocate(Xg(0:Nx),XGc(Ncx),LeftBounds(Ncx),RightBounds(Ncx))
    end if
    if ((specialPoint.gt.rightBound) .or. (specialPoint.lt.leftBound)) then
        print *,'Error in defining the grid: special point outside the interval', leftBound,specialPoint,rightBound
        stop 'Error initializing the grid'
    end if
    call initGrid(specialPoint,densityLeft,densityRight)
    call setCOLLOC
    call SetBounds
  end subroutine InitDiscretization
!
 subroutine SetBounds
 integer :: i,j
 double precision :: tl,tr,tstl,tstr,dt
 do i=1,Ncx
     do j=0,Nx-1
         dt=(XG(j+1)-XG(j))*0.1d0
         tl=XG(j)-dt
         tr=XG(j)+dt
         tstl=abs(HSplineBC(i-1,tl))
         tstr=abs(HSplineBC(i-1,tr))
         !print *,'??? Left:',i,j,LeftBounds(i),tl,tstl,tr,tstr
         if ((tstl.lt.1.0d-10).and.(tstr.gt.0.0d0)) then
             LeftBounds(i)=XG(j)
             !print *,'ACCEPTED'
             !print *,'Left:',i,j,LeftBounds(i),tl,tstl,tr,tstr
         end if
         tl=XG(j+1)-dt
         tr=XG(j+1)+dt
         tstl=abs(HSplineBC(i-1,tl))
         tstr=abs(HSplineBC(i-1,tr))
         !print *,'???right:',i,j+1,RightBounds(i),tl,tstl,tr,tstr
         if ((tstr.lt.1.0d-10).and.(tstl.gt.0.0d0)) then
             RightBounds(i)=XG(j+1)
             !print *,'ACCEPTED'
             !print *,'right:',i,j+1,RightBounds(i),tl,tstl,tr,tstr
         end if
     end do
 end do
 end subroutine SetBounds
!
 double precision function xmap(t,c0,c1,cs,ts)
 double precision :: t,c0,c1,cs,ts
 double precision :: res
   res=0.0
    if (t.lt.ts) then
      res=(-c0+cs)*t**2/ts/2+c0*t
    else
      res=(-c0+cs)*ts/2+c0*ts
      res=res+(-cs+c1)*(t**2-ts**2)/(1-ts)/2+cs*(t-ts)/(1-ts)-c1*(t-ts)*ts/(1-ts)
    end if
   xmap=res
 end function xmap
!
  subroutine initGrid(specialPoint,densityLeft,densityRight)
!    character, intent(in) :: confName
    double precision :: specialPoint,densityLeft,densityRight
    double precision :: t,dt,a0,a1
    double precision :: c0,c1,cs,ts,sp,s
    integer :: i
    a0=1.0d0/densityLeft
    a1=1.0d0/densityRight
    sp=specialPoint-Xl
    s=Xr-Xl
    cs=2*(s+((a1-a0)/(a0+1))*sp)/(1+a1)
    c0=a0*cs
    c1=a1*cs
    ts=2*sp/cs/(a0+1)
    dt=1.0d0/Nx
    do i=0,Nx
      t=i*dt
      XG(i)=xmap(t,c0,c1,cs,ts)+Xl
    end do
  end subroutine initGrid

  subroutine resetGrid(xg1,nx1)
      integer :: nx1
      double precision :: xg1(0:nx1)
      integer :: i
      if (nx1.ne.nx) stop 'Sent grind size is inconsistent with initialized size'
      do i=0,nx
          xg(i)=xg1(i)
      end do
      call setColloc
  end subroutine resetGrid
!
      SUBROUTINE setCOLLOC
      INTEGER :: I
      integer :: ir
      double precision :: SH1,SH2
      double precision :: sh31,sh32,sh33
      sh31=(1.0d0-sqrt(0.6d0))*0.5d0
      sh32=0.5d0
      sh33=(1.0d0+sqrt(0.6d0))*0.5d0
      SH1=(1.0D0-1.0D0/DSQRT(3.0D0))*0.5D0
      SH2=(1.0D0+1.0D0/DSQRT(3.0D0))*0.5D0
      if (BCLeft.ne.-1) then
          DO I=1,NX
             XGC(2*I-1)=(XG(I)-XG(I-1))*SH1+XG(I-1)
             XGC(2*I)=(XG(I)-XG(I-1))*SH2+XG(I-1)
          END DO
          ir=2*NX-1
      else
          XGC(1)=(XG(1)-XG(0))*SH31+XG(0)
          XGC(2)=(XG(1)-XG(0))*SH32+XG(0)
          XGC(3)=(XG(1)-XG(0))*SH33+XG(0)
          DO I=2,NX
             XGC(2*I)=(XG(I)-XG(I-1))*SH1+XG(I-1)
             XGC(2*I+1)=(XG(I)-XG(I-1))*SH2+XG(I-1)
          END DO
          ir=2*NX
      end if
      if (BCRight.eq.-1) then
          XGC(ir)=(XG(Nx)-XG(Nx-1))*SH31+XG(Nx-1)
          XGC(ir+1)=(XG(Nx)-XG(Nx-1))*SH32+XG(Nx-1)
          XGC(ir+2)=(XG(Nx)-XG(Nx-1))*SH33+XG(Nx-1)
      end if
      END SUBROUTINE setCOLLOC

!**************************************************************
      DOUBLE PRECISION FUNCTION HSPlineBC(I,X)
      integer, intent(in) :: I
      integer :: i1
      double precision, intent(in) :: x
      ! this function takes care of appropriate python idexing i->-1 of the basis
      ! for the given boundary conditions.
      IF (BCLeft.eq. -1) then
        i1=i
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      Else if (BCLeft.eq.0) then
        i1=i+1
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      ELse if (BCLeft.eq.1) then
        i1=i+1
        if (i.eq.0) i1=0
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=I1+1
      end if
      HSPlineBC=HSPline(I1,X)
      END FUNCTION HSPlineBC

      DOUBLE PRECISION FUNCTION HSPlineBC1(I,X)
      integer, intent(in) :: I
      integer :: i1
      double precision, intent(in) :: x
      ! this function takes care of appropriate python idexing i->-1 of the basis
      ! for the given boundary conditions.
      IF (BCLeft.eq. -1) then
        i1=i
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      Else if (BCLeft.eq.0) then
        i1=i+1
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      ELse if (BCLeft.eq.1) then
        i1=i+1
        if (i.eq.0) i1=0
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=I1+1
      end if
      HSPlineBC1=HSPline1(I1,X)
      END FUNCTION HSPlineBC1

      DOUBLE PRECISION FUNCTION HSPlineBC2(I,X)
      integer, intent(in) :: I
      integer :: i1
      double precision, intent(in) :: x
      ! this function takes care of appropriate python idexing i->-1 of the basis
      ! for the given boundary conditions.
      IF (BCLeft.eq. -1) then
        i1=i
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      Else if (BCLeft.eq.0) then
        i1=i+1
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=i1+1
      ELse if (BCLeft.eq.1) then
        i1=i+1
        if (i.eq.0) i1=0
        if ((BCRight.eq.0).and.(I.eq.Ncx-1)) i1=I1+1
      end if
      HSPlineBC2=HSPline2(I1,X)
      END FUNCTION HSPlineBC2

!-----------------------------------------------------------------!
      DOUBLE PRECISION FUNCTION HSPline(I,X)
      integer, intent(in) :: I
      double precision, intent(in) :: x
      double precision :: Spl
      integer NX2, IPAR
      double precision X1I,XI,XI1
      integer :: i1i,ii,ii1
!      write(01,*) 'PARAMETERS OF HSPLIN'
!      write(01,*) 'X,I,NX=',X,I,NX
      NX2=NX*2
      IPAR=I/2
      IPAR=I-IPAR*2
      IF (IPAR.EQ.0) THEN
         ii=I/2
         if (ii.gt.NX) ii=NX
         i1i=ii-1
         ii1=ii+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
         X1I=XG(i1i)
         XI=XG(ii)
         XI1=XG(ii1)
         SPL=PHI(X,X1I,XI,XI1)
      ELSE
        i1i=(I-1)/2-1
        ii=(I-1)/2
        ii1=(I-1)/2+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
        X1I=XG(i1i)
        XI=XG(ii)
        XI1=XG(ii1)
        SPL=PSI(X,X1I,XI,XI1)
      ENDIF
      HSpline=SPL
      END FUNCTION HSpline

      DOUBLE PRECISION FUNCTION HSPline1(I,X)
      integer, intent(in) :: I
      double precision, intent(in) :: x
      double precision :: Spl1
      integer NX2, IPAR
      double precision X1I,XI,XI1
      integer :: i1i,ii,ii1
!      write(01,*) 'PARAMETERS OF HSPLIN'
!      write(01,*) 'X,I,NX=',X,I,NX
      NX2=NX*2
      IPAR=I/2
      IPAR=I-IPAR*2
      IF (IPAR.EQ.0) THEN
         ii=I/2
         if (ii.gt.NX) ii=NX
         i1i=ii-1
         ii1=ii+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
         X1I=XG(i1i)
         XI=XG(ii)
         XI1=XG(ii1)
         SPL1=PHI1(X,X1I,XI,XI1)
      ELSE
        i1i=(I-1)/2-1
        ii=(I-1)/2
        ii1=(I-1)/2+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
        X1I=XG(i1i)
        XI=XG(ii)
        XI1=XG(ii1)
        SPL1=PSI1(X,X1I,XI,XI1)
      ENDIF
      HSpline1=SPL1
      END FUNCTION HSpline1

      DOUBLE PRECISION FUNCTION HSPline2(I,X)
      integer, intent(in) :: I
      double precision, intent(in) :: x
      double precision :: Spl2
      integer NX2, IPAR
      double precision X1I,XI,XI1
      integer :: i1i,ii,ii1
!      write(01,*) 'PARAMETERS OF HSPLIN'
!      write(01,*) 'X,I,NX=',X,I,NX
      NX2=NX*2
      IPAR=I/2
      IPAR=I-IPAR*2
      IF (IPAR.EQ.0) THEN
         ii=I/2
         if (ii.gt.NX) ii=NX
         i1i=ii-1
         ii1=ii+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
         X1I=XG(i1i)
         XI=XG(ii)
         XI1=XG(ii1)
         SPL2=PHI2(X,X1I,XI,XI1)
      ELSE
        i1i=(I-1)/2-1
        ii=(I-1)/2
        ii1=(I-1)/2+1
        if (i1i.lt.0) i1i=0
        if (ii1.gt.nx) ii1=nx
        X1I=XG(i1i)
        XI=XG(ii)
        XI1=XG(ii1)
        SPL2=PSI2(X,X1I,XI,XI1)
      ENDIF
      HSpline2=SPL2
      END FUNCTION HSpline2
!-------------------------------------------------------------------!
      real*8 FUNCTION PHI(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      IF (X.LT.X1I) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.EQ.XI) THEN
        PHI=1.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
        IF (HK.EQ.0.0D0) THEN
          HK=XI1-XI
          RES=2*(X-XI+0.5D0*HK)*(X-XI1)*(X-XI1)/(HK*HK*HK)
        ELSE
          RES=-2*(X-X1I)*(X-X1I)*(X-XI-0.5D0*HK)/(HK*HK*HK)
        END IF
      ELSE
        HK=XI1-XI
        IF (HK.EQ.0.0D0) THEN
          HK=XI-X1I
          RES=-2*(X-X1I)*(X-X1I)*(X-XI-0.5D0*HK)/(HK*HK*HK)
        ELSE
          RES=2*(X-XI+0.5D0*HK)*(X-XI1)*(X-XI1)/(HK*HK*HK)
        END IF
      END IF
      PHI=RES
      RETURN
      END FUNCTION PHI
      real*8 FUNCTION PSI(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      IF (X.LT.X1I) THEN
        PSI=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PSI=0.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
        IF (HK.EQ.0.0D0) THEN
          HK=XI1-XI
          RES=(X-XI)*(X-XI1)*(X-XI1)/(HK*HK)
        ELSE
          RES=(X-X1I)*(X-X1I)*(X-XI)/(HK*HK)
        END IF
      ELSE
        HK=XI1-XI
        IF (HK.EQ.0.0D0) THEN
          HK=XI-X1I
          RES=(X-X1I)*(X-X1I)*(X-XI)/(HK*HK)
        ELSE
          RES=(X-XI)*(X-XI1)*(X-XI1)/(HK*HK)
        END IF
      END IF
      PSI=RES
      RETURN
      END FUNCTION PSI
      
      double precision function phi1(x, x1i, xi, xi1)
      double precision :: x, x1i,xi,xi1
      double precision :: res, HK, DL, DR
      RES=0.0d0
      if (x.lt.x1i) then
        res=0.0d0
      else if (x.gt.xi1) then
        res=0.0d0
      else if (x.eq.xi ) then
        res=0.0d0
      else if (x.lt.xi) then
        HK=xi-x1i
        DL=x-x1i
        DR=x-xi
        res=-2*(DL*(DR-0.5*HK)+DL*(DR-0.5*HK)+DL*DL)/(HK*HK*HK)
      else
        HK=xi1-xi
        DL=x-xi
        DR=x-xi1
        res=2*(DR*DR+(DL+0.5*HK)*DR+(DL+0.5*HK)*DR)/(HK*HK*HK);
      end if
      phi1=res  
      end function phi1

      double precision function psi1(x, x1i, xi, xi1)
      double precision :: x, x1i,xi,xi1
      double precision :: res, HK, DL, DR
      res=0.0d0
      if (x.lt.x1i) then
        res=0.0d0
      else if (x.gt.xi1) then
        res=0.0d0
      else if (x.eq.xi ) then
        res=1.0d0
      else if (x.lt.xi) then
        HK=xi-x1i
        DL=x-x1i
        DR=x-xi
        res=(DL*DR+DL*DR+DL*DL)/(HK*HK)
      else
        HK=xi1-xi
        DL=x-xi
        DR=x-xi1
        res=(DR*DR+DL*DR+DL*DR)/(HK*HK)
      end if
      psi1=res
      end function psi1
      
      real*8 FUNCTION PHI2(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      RES=0.0d0
      IF (X.LT.X1I) THEN
        PHI2=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PHI2=0.0D0
        RETURN
      END IF
      IF (X.LE.XI) THEN
        HK=XI-X1I
        IF (HK.NE.0.0D0) RES=-4*(X-XI-0.5D0*HK+2*(X-X1I))/(HK*HK*HK)
      ELSE
        HK=XI1-XI
        IF (HK.NE.0.0D0) RES=4*(X-XI+0.5D0*HK+2*(X-XI1))/(HK*HK*HK)
      END IF
      PHI2=RES
      RETURN
      END FUNCTION PHI2
      real*8 FUNCTION PSI2(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      RES=0.0d0
      IF (X.LT.X1I) THEN
        PSI2=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PSI2=0.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
        IF (HK.NE.0.0D0) RES=2*(2*(X-X1I)+(X-XI))/(HK*HK)
      ELSE
        HK=XI1-XI
        IF (HK.NE.0.0D0) RES=2*((X-XI)+2*(X-XI1))/(HK*HK)
      END IF
      PSI2=RES
      RETURN
      END  FUNCTION PSI2

      SUBROUTINE LOCAT(X,I1,I2)
      ! returns indexes of the first and the last non-zero spline functions 
      ! for a given point x
      integer i1,i2
      double precision X
      integer i0
      DO I0=1,NX
        IF ((X.GT.XG(I0-1)).AND.(X.LT.XG(I0))) GOTO 3
      END DO
 3    I1=I0*2-2
      I2=I0*2+1
      IF (I1.LT.1) THEN
        I1=1
        I2=3
      END IF
      IF (I2.GT.NX*2) THEN 
        I2=NX*2
        I1=NX*2-2
      END IF
      RETURN
      END  SUBROUTINE LOCAT
    
end module HermitSplines
