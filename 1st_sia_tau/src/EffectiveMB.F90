
       FUNCTION EffectiveMB(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut
       REAL(kind=dp) :: rhowater,rhoice

          rhowater = 1.028
          rhoice   = 0.917

          IF ((VarIn(1).GT.-100.0)) THEN
             VarOut=(rhowater/rhoice*VarIn(1)-VarIn(2))*(+1.0)
!+1.8675
          ELSE
             VarOut=0._dp
          END IF

        End FUNCTION EffectiveMB

!#
       FUNCTION HPassive(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

          IF ((VarIn(1).GT.-100.0).AND.(VarIn(2).GT.-100.0)) THEN
             VarOut=-1._dp
          ELSE
             VarOut=+1._dp
          END IF

        End FUNCTION HPassive
!#
       FUNCTION Velocity(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(3),VarOut
       REAL(kind=dp) :: eps=1.0e+01,Vlim=5.0

          IF (VarIn(3).GT.Vlim) THEN
             VarOut=VarIn(1)
          ELSE
             VarOut=-eps*VarIn(2)
          END IF

        End FUNCTION Velocity
       FUNCTION Vdir(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(5),VarOut
       REAL(kind=dp) :: eps=1.0e+01,Vlim=5.0
       REAL(kind=dp) :: u,norm

          IF (VarIn(5).GT.Vlim) THEN
             !norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             norm=1._dp
             u=VarIn(1)
          ELSE
             norm=sqrt(VarIn(3)*VarIn(3)+VarIn(4)*VarIn(4))
             u=-eps*VarIn(3)
          END IF

             if (norm.GT.1.0e-6) Then
               VarOut=u/norm
             else
             Varout = 0._dp
             endif

        End FUNCTION Vdir
!#
       FUNCTION Variance(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut
       REAL(kind=dp) :: Vlim=25.0

!          IF (4.0/5.0*VarIn.GT.4.0/5.0*Vlim) THEN
          IF (VarIn.GT.Vlim) THEN
             VarOut=1000.0_dp
          ELSE
             VarOut=1000.0_dp
          END IF

        End FUNCTION Variance
!#
       FUNCTION VarianceSMB(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

!          IF (VarIn.GT.Vlim) THEN
             VarOut=2.0_dp
!          ELSE
!             VarOut=1000.0_dp
!          END IF
          !IF (VarIn.GT.0) THEN
          !    VarOut = VarOut*100.
          !ENDIF

        End FUNCTION VarianceSMB
!#
       FUNCTION Vdirx(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(5),VarOut
       REAL(kind=dp) :: eps=0.5e+01,Vlim=5.0
       REAL(kind=dp) :: u,norm,a1,a2
       REAL(kind=dp) :: aa,bb
       REAL(kind=dp) :: Arate,rhoi,nflow
       REAL(kind=dp) :: secperyear,grav
 
       secperyear = 60*60*24*365.25
       rhoi       = 917.0
       grav       = 9.81
       nflow      = 3.0
       Arate      = 100.0*1.0e-25

          IF (VarIn(3).GT.Vlim.and.VarIn(3).lt.365.and.VarIn(5).lt.800) THEN
             !norm=1._dp
             !u=abs(VarIn(3))*sin((VarIn(4)+0.0)*3.14159265358979/180.0)
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             u=-abs(VarIn(3))*VarIn(1)
          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
!             a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
!             a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
             u=-eps*VarIn(1)
          END IF

          aa = 300.0/2.0
          bb = aa/630.0
          norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
          if(norm.eq.0)then
          norm=1.0_dp
          u = abs(aa - bb*VarIn(5))*sin((VarIn(4)+0.0)*3.14159265358979/180.0)
          else
          u = (-1)*abs(aa - bb*VarIn(5))*VarIn(1)
          endif

          VarOut = u/norm
!             if (norm.GT.1.0e-6) Then
!             VarOut=4.0/5.0*u/norm
!!             VarOut=u/norm
!             else
!             VarOut = sign(1.0,VarIn(1))*4.0/5.0/sqrt(2.0)
!             !Varout = 4.0/5.0/sqrt(2)
!             endif

        End FUNCTION Vdirx
!#
       FUNCTION Vdiry(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(5),VarOut
       REAL(kind=dp) :: eps=0.5e+01,Vlim=5.0
       REAL(kind=dp) :: u,norm,a1,a2
       REAL(kind=dp) :: aa,bb
       REAL(kind=dp) :: Arate,rhoi,nflow
       REAL(kind=dp) :: secperyear,grav
 
       secperyear = 60*60*24*365.25
       rhoi       = 917.0
       grav       = 9.81
       nflow      = 3.0
       Arate      = 100.0*1.0e-25

          IF (VarIn(3).GT.Vlim.and.VarIn(3).lt.365.and.VarIn(5).lt.800) THEN
             !norm=1._dp
             !u=abs(VarIn(3))*cos((VarIn(4)+0.0)*3.14159265358979/180.0)
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             u=-abs(VarIn(3))*VarIn(1)
          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
!             a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
!             a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
             u=-eps*VarIn(1)
          END IF

          aa = 300.0/2.0
          bb = aa/630.0
          norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
          if(norm.eq.0)then
          norm=1.0_dp
          u = abs(aa - bb*VarIn(5))*cos((VarIn(4)+0.0)*3.14159265358979/180.0)
          else
          u = (-1)*abs(aa - bb*VarIn(5))*VarIn(1)
          endif

          VarOut = u/norm

!             if (norm.GT.1.0e-6) Then
!             VarOut=4.0/5.0*u/norm
!!             VarOut=u/norm
!             else
!             VarOut = sign(1.0,VarIn(1))*4.0/5.0/sqrt(2.0)
!             endif

        End FUNCTION Vdiry
!#
       FUNCTION computeSIAthickness(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut
       REAL(kind=dp) :: Arate,rhoice,nflow
       REAL(kind=dp) :: secperyear,grav
       REAL(kind=dp) :: vmag,slope,hsia
       REAL(kind=dp) :: a1,a2,a3,b1,b2

       secperyear = 60*60*24*365.25
       !rhowater   = 1028.
       rhoice     = 917.0
       grav       = 9.81
       nflow      = 3.0
       Arate      = 100.0*1.0e-25

       vmag  = sqrt(VarIn(3)**2 + VarIn(4)**2) 
       slope = sqrt(VarIn(1)**2 + VarIn(2)**2)

       IF(slope.gt.0.002)THEN 
       a1    = vmag/secperyear/(slope**nflow)
       a2    = (nflow+1)/(2*Arate)*1.0/((rhoice*grav)**nflow)
       a3    = 1/(nflow+1)
       b1    = a1**a3
       b2    = a2**a3
       Hsia  = b1*b2
       ELSE
       Hsia  = -1.0
       ENDIF

       VarOut = Hsia

       End FUNCTION computeSIAthickness
!#
!#
       FUNCTION computeSMB(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(1),VarOut
       REAL(kind=dp) :: aa,bb

       aa = 0.1/300.0
       bb = 300.0

       VarOut = aa*(VarIn(1)-bb)

       VarOut = VarIn(1)

       End FUNCTION computeSMB
!#
       FUNCTION vmag(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(1),VarOut
       REAL(kind=dp) :: aa,bb

       aa = 300.0
       bb = aa/630.0

       VarOut = abs(aa - bb*VarIn(1))



       End FUNCTION vmag
!#



