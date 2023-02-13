!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   jrnmm : Jansen-Rit Neural Mass Model for a cortical neuron.
!       Here we analyze the homogeneous component of properly normalized
!       network, which boils down to study a self-coupled JR-NMM    
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION y0,y1,y2,y3,y4,y5,y6,a1,a0,b1,b0,c1,c2,c3,c4,p
      DOUBLE PRECISION e0,r,v0,c,eps,sigm,dsigm,s,x1,x2
      INTEGER i,j

       dsigm(s)=2*e0*r*exp(r*(v0-s))/((1+exp(r*(v0-s)))**2)
       sigm(s)=2*e0/(1+exp(r*(v0-s)))
       y0=U(1)
       y1=U(2)
       y2=U(3)
       y3=U(4)
       y4=U(5)
       y5=U(6)

       x1=U(7)
       x2=U(8)
        
       C=135.0

        a1 = 3.25
        a0 = 100.0
        b1 = 22.0
        b0 = 50.0
        c1= c
        c2= 0.8*c
        c3= 0.25*c
        c4= 0.25*c

        p =PAR(1) 
        eps = PAR(2)

        e0=2.5
        r=0.56
        v0=6.0


       F(1)= y3
       F(2)= y4
       F(3)= y5
       F(4)= a1*a0*sigm(y1-y2+x1) - 2*a0*y3 - a0**2*y0
       F(5)= a1*a0*(p + c2*sigm(C1*y0)) - 2*a0*y4 - a0**2*y1
       F(6)= b1*b0*c4*sigm(c3*y0+x1) -2*b0*y5 - b0**2*y2

       F(7)= x2 ! \dot x1
       F(8)= a1*a0*eps*sigm(y1-y2) - 2*a0*x2 - a0**2*x1 ! \dot x2


      IF(IJAC.EQ.0)RETURN

        !do i = 1, 6
        !        do j = 1, 6
        !                DFDU(i,j)=0;
        !        end do
        !end do

        !DFDU(1,4)=1;
        !DFDU(2,5)=1;
        !DFDU(3,6)=1;

        !DFDU(4,4)=-2*a0;
        !DFDU(5,5)=-2*a0;
        !DFDU(6,6)=-2*b0;

        !DFDU(4,1)=-a0**2;
        !DFDU(5,2)=-a0**2;
        !DFDU(6,3)=-b0**2;

        !DFDU(4,2)= a0*a1*dsigm(y1-y2);
        !DFDU(4,3)=-a0*a1*dsigm(y1-y2);
        !DFDU(5,1)= a0*a1*c1*c2*dsigm(c1*y0);
        !DFDU(6,1)= b0*b1*c3*c4*dsigm(c3*y0);

        !DFDU(5,2) =DFDU(5,2) + a1*a0*eps*dsigm(y1-y2);
        !DFDU(5,3) =DFDU(5,3) - a1*a0*eps*dsigm(y1-y2);


      IF(IJAC.EQ.1)RETURN

       ! do i = 1, 6
       !         do j = 1, 2
       !                 DFDP(i,j)=0;
       !         end do
       ! end do

       !DFDP(5,1)=a1*a0
       !DFDP(5,2)=a1*a0*sigm(y1-y2)

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=10.0 ! P
      PAR(2)=20.0 ! EPSILON

      U(1)=7.7E-02
      U(2)=-6.05E-02
      U(3)=0.01
      U(4)=0.01
      U(5)=0.01
      U(6)=0.01
      U(7)=0.01
      U(8)=0.01


      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MAX(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: PCLU_Y(0:NTST*NCOL)

      PCLU_Y=U(2,:)-U(3,:)

      GETUY_MAX=MAXVAL(PCLU_Y)

      END FUNCTION GETUY_MAX
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: PCLU_Y(0:NTST*NCOL)

      PCLU_Y=U(2,:)-U(3,:)

      GETUY_MIN=MINVAL(PCLU_Y)

      END FUNCTION GETUY_MIN
!-----------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUY_MAX,GETUY_MIN
      INTEGER NDX,NCOL,NTST

      NDX=NINT(GETP('NDX',0,U))
      NTST=NINT(GETP('NTST',0,U))
      NCOL=NINT(GETP('NCOL',0,U))

      PAR(3)=GETUY_MIN(U,NDX,NTST,NCOL)
      PAR(4)=GETUY_MAX(U,NDX,NTST,NCOL)

      END SUBROUTINE PVLS
