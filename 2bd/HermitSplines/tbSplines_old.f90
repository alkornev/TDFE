!
module tbSplines
!
!  This module takes care of the system discretization. 
!  We employ tb splines and orthogonal collocations.
!  We assume the same grid for all the channels.
!  the preliminary version is 1-channel.
!
  implicit none
  integer :: Nx,Ncx
  integer :: BCLeft=0,BCRight=0
  integer :: SplineType=0
  double precision :: Xl,Xr
  double precision, allocatable :: XG(:),XGc(:)
  double precision, allocatable :: LeftBounds(:),RightBounds(:)
  logical :: ifAllocated=.false.
  contains
!
  subroutine InitDiscretization(NPoints,Grid,iBCLeft,iBCRight, iSplineType)
    integer, intent(in) :: iBCLeft,iBCRight,iSplineType
    integer, intent(in) :: NPoints
    double precision, intent(in) :: Grid(0:NPoints)
!    Xl=leftBound
!    Xr=rightBound
    BCLeft=iBCLeft
    BCRight=iBCRight
    SplineType=iSplineType
    NX=NPoints
    if(SplineType.eq.0) then
      Ncx=2*Nx
    else if(SplineType.eq.1) then
      Ncx=3*Nx
    else 
      stop 'SplineType can be only 0 or 1'
    end if
    if (abs(BCRight).gt.1) then
      stop 'Error initializing the grid, BCRight can be only 0,1 or -1'
    end if
    if (abs(BCLeft).gt.1) then
      stop 'Error initializing the grid, BCleft can be only 0,1 or -1'
    end if
    if (BCRight.eq.-1) Ncx=Ncx+1
    if (BCLeft.eq.-1) then
	if (SplineType.eq.1) then
	  Ncx=Ncx+2
	else if(SplineType.eq.0) then
	  Ncx=Ncx+1	
	end if
    end if
    if (.not.ifAllocated) then
      allocate(Xg(0:Nx),XGc(Ncx),LeftBounds(Ncx),RightBounds(Ncx))
      ifAllocated=.true.
    else if (size(XGc).ne.Ncx) then
        deallocate(Xg,XGc,LeftBounds,RightBounds)
        allocate(Xg(0:Nx),XGc(Ncx),LeftBounds(Ncx),RightBounds(Ncx))
    end if
    call resetGrid(Grid,NPoints)
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
         if ((tstl.lt.1.0d-10).and.(tstr.gt.0.0d0)) then
             LeftBounds(i)=XG(j)
         end if
         tl=XG(j+1)-dt
         tr=XG(j+1)+dt
         tstl=abs(HSplineBC(i-1,tl))
         tstr=abs(HSplineBC(i-1,tr))
         if ((tstr.lt.1.0d-10).and.(tstl.gt.0.0d0)) then
             RightBounds(i)=XG(j+1)
         end if
     end do
 end do
 end subroutine SetBounds
!
!

  subroutine resetGrid(xg1,nx1)
      integer :: nx1
      double precision :: xg1(0:nx1)
      integer :: i
      if (nx1.ne.nx) stop 'Sent grid size is inconsistent with initialized size'
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
      double precision :: sh41, sh42, sh43, sh44
      double precision :: sh51, sh52, sh53, sh54, sh55
      double precision :: sh61, sh62, sh63, sh64, sh65, sh66

      SH1= (1.0D0-1.0D0/DSQRT(3.0D0))*0.5D0!0.21132487d0
      SH2=(1.0D0+1.0D0/DSQRT(3.0D0))*0.5D0!0.78867513d0

      sh31=(1.0d0-sqrt(0.6d0))*0.5d0
      sh32=0.5d0
      sh33=(1.0d0+sqrt(0.6d0))*0.5d0

      sh41=(1.0d0-0.86113631d0)*0.5d0!1594053d0
      sh42=(1.0d0-0.33998104d0)*0.5d0!3584856d0
      sh43=(1.0d0+0.33998104d0)*0.5d0!3584856d0
      sh44=(1.0d0+0.86113631d0)*0.5d0!1594053d0

      sh51=(1.0d0-0.90617984d0)*0.5d0!5938664d0
      sh52=(1.0d0-0.53846931d0)*0.5d0!0105683d0
      sh53=0.5d0
      sh54=(1.0d0+0.53846931d0)*0.5d0!0105683d0)
      sh55=(1.0d0+0.90617984d0)*0.5d0!5938664d0

      sh61=-0.93246951d0!4203152d0
      sh62=-0.66120938d0!6466265d0
      sh63=-0.23861918d0!6083197d0
      sh64=0.23861918d0!6083197d0
      sh65=0.66120938d0!6466265d0
      sh66=0.93246951d0!4203152d0
      sh61=(1.0d0+sh61)*0.5d0
      sh62=(1.0d0+sh62)*0.5d0
      sh63=(1.0d0+sh63)*0.5d0
      sh64=(1.0d0+sh64)*0.5d0
      sh65=(1.0d0+sh65)*0.5d0
      sh66=(1.0d0+sh66)*0.5d0

      if (SplineType.eq.1) then
      	if (BCLeft.ne.-1) then
          DO I=1,NX
             XGC(3*I-2)=(XG(I)-XG(I-1))*SH31+XG(I-1)
             XGC(3*I-1)=(XG(I)-XG(I-1))*SH32+XG(I-1)
             XGC(3*I)=(XG(I)-XG(I-1))*SH33+XG(I-1)
          END DO
          ir=3*NX-2
        else
            XGC(1)=(XG(1)-XG(0))*SH51+XG(0)
            XGC(2)=(XG(1)-XG(0))*SH52+XG(0)
            XGC(3)=(XG(1)-XG(0))*SH53+XG(0)
            XGC(4)=(XG(1)-XG(0))*SH54+XG(0)
            XGC(5)=(XG(1)-XG(0))*SH55+XG(0)
            DO I=2,NX-1
                 XGC(3*I)=(XG(I)-XG(I-1))*SH31+XG(I-1)
                 XGC(3*I+1)=(XG(I)-XG(I-1))*SH32+XG(I-1)
                 XGC(3*I+2)=(XG(I)-XG(I-1))*SH33+XG(I-1)
            END DO
            ir=3*NX
        end if
        if (BCRight.eq.-1) then
              XGC(ir)=(XG(Nx)-XG(Nx-1))*SH41+XG(Nx-1)
              XGC(ir+1)=(XG(Nx)-XG(Nx-1))*SH42+XG(Nx-1)
              XGC(ir+2)=(XG(Nx)-XG(Nx-1))*SH43+XG(Nx-1)
              XGC(ir+3)=(XG(Nx)-XG(Nx-1))*SH44+XG(Nx-1)
        end if
        if (BCLeft.eq.-1.and.BCRight.eq.-1.and.NX.eq.1) then
              XGC(1)=(XG(1)-XG(0))*SH61+XG(0)
              XGC(2)=(XG(1)-XG(0))*SH62+XG(0)
              XGC(3)=(XG(1)-XG(0))*SH63+XG(0)
              XGC(4)=(XG(1)-XG(0))*SH64+XG(0)
              XGC(5)=(XG(1)-XG(0))*SH65+XG(0)
              XGC(6)=(XG(1)-XG(0))*SH66+XG(0)
        end if
      else if (SplineType.eq.0) then
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
	NX2=NX
      	IPAR=I/1
      	IPAR=I-IPAR*1
       	ii=I/1
           if (ii.gt.NX) ii=NX
           i1i=ii-1
           ii1=ii+1
           if (i1i.lt.0) i1i=0
           if (ii1.gt.nx) ii1=nx
           X1I=XG(i1i)
           XI=XG(ii)
           XI1=XG(ii1)
           SPL=PHI(X,X1I,XI,XI1)
        END IF	
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
	NX2=NX
      	IPAR=I/1
      	IPAR=I-IPAR*1
       	ii=I/1
           if (ii.gt.NX) ii=NX
           i1i=ii-1
           ii1=ii+1
           if (i1i.lt.0) i1i=0
           if (ii1.gt.nx) ii1=nx
           X1I=XG(i1i)
           XI=XG(ii)
           XI1=XG(ii1)
           SPL=TPHI1(X,X1I,XI,XI1)
        END IF	
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
      if(SplineType.eq.0) then
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
          SPL2=TPHI2(X,X1I,XI,XI1)
        END IF

      HSpline2=SPL2
      END FUNCTION HSpline2
!-------------------------------------------------------------------!
      real*8 FUNCTION TPHI(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision :: X, X1I,XI,XI1
      double precision :: RES, HK, DL, DR
      IF (X.LT.X1I) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.EQ.XI) THEN
        PHI=2.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
	DL=X-X1I
	DR=X-XI
	RES=SIN(DL/2.0d0)/SIN(HK/2.0d0)
      ELSE
	HK=XI1-XI
	DL=X-XI
	DR=X-XI1	
        RES=SIN(-DR/2.0d0)/SIN(HK/2.0d0)
      END IF
      PHI = RES
      RETURN
      END FUNCTION PHI

      real*8 FUNCTION TPHI1(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision :: X, X1I,XI,XI1
      double precision :: RES, HK, DL, DR
      IF (X.LT.X1I) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.EQ.XI) THEN
        PHI=2.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
	DL=X-X1I
	DR=X-XI
	RES=0.5d0*COS(DL/2.0d0)/SIN(HK/2.0d0)
      ELSE
	HK=XI1-XI
	DL=X-XI
	DR=X-XI1	
        RES=-0.5d0*COS(-DR/2.0d0)/SIN(HK/2.0d0)
      END IF
      PHI = RES
      RETURN
      END FUNCTION PHI

      real*8 FUNCTION TPHI2(X,X1I,XI,XI1)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision :: X, X1I,XI,XI1
      double precision :: RES, HK, DL, DR
      IF (X.LT.X1I) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.GT.XI1) THEN
        PHI=0.0D0
        RETURN
      END IF
      IF (X.EQ.XI) THEN
        PHI=2.0D0
        RETURN
      END IF
      IF (X.LT.XI) THEN
        HK=XI-X1I
	DL=X-X1I
	DR=X-XI
	RES=-0.25d0*SIN(DL/2.0d0)/SIN(HK/2.0d0)
      ELSE
	HK=XI1-XI
	DL=X-XI
	DR=X-XI1	
        RES=0.25d0*SIN(-DR/2.0d0)/SIN(HK/2.0d0)
      END IF
      PHI = RES
      RETURN
      END FUNCTION PHI

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
    
end module tbSplines
