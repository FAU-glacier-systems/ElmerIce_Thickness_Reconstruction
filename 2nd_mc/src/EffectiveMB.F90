
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

          !IF ((VarIn(1).GT.-100.0).AND.(VarIn(2).GT.-100.0)) THEN
          !   !VarOut=rhowater/rhoice*VarIn(1)
          !   VarOut=(rhowater/rhoice*VarIn(1)-VarIn(2))*(+1.0)
          !   !VarOut=(+1.0/1.0)*((rhowater/rhoice*VarIn(1)-0.05*VarIn(2))+0.0)
!         !    VarOut=VarIn(1)-VarIn(2)
          !ELSE
          !   VarOut=0._dp
          !END IF

! Schaefer et al. (2013)
          IF (VarIn(1).lt.1200)THEN
          VarOut = (0.013*VarIn(1) - 16.2) - VarIn(2)
          ELSE
          VarOut = (0.0084*VarIn(1) - 9.8) - VarIn(2)
          ENDIF
          ! 1200
          !IF (VarIn(1).lt.1200)THEN
          !VarOut = 0.013*(VarIn(1) - 1200.0)
          !ELSE
          !VarOut = 0.0084*(VarIn(1) - 1200.0)
          !ENDIF
          !VarOut = VarOut*1.2 -VarIn(2)

! Koppes et al. (2011)
!          IF (VarIn(1).lt.1800)THEN
!          VarOut = (0.013*VarIn(1) - 16.2) - VarIn(2)
!          ELSE
!          VarOut = ((-1.0)*0.0022*VarIn(1) + 11.7) - VarIn(2)
!          ENDIF
! Collao et al. (2017): best fit Rafael
!          IF (VarIn(1).lt.1250)THEN
!          VarOut = (0.013*VarIn(1) - 16.2) - VarIn(2)
!          ELSE
!          VarOut = (5.95/2750.0*VarIn(1) + (6.0-4000.0/2750.0*5.95)) - VarIn(2)
!          ENDIF


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
       REAL(kind=dp) :: VarIn,VarOut
       REAL(kind=dp) :: eps=1.0e+01,Vlim=5.0

          !IF (VarIn(3).GT.Vlim) THEN
             VarOut=VarIn
          !ELSE
          !   VarOut=-eps*VarIn(2)
          !END IF

        End FUNCTION Velocity
!
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

        End FUNCTION VarianceSMB
!#
       FUNCTION VarianceH(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(3),VarOut

          IF (VarIn(1).LT.100.or.VarIn(2).lt.50) THEN
!          IF (VarIn(1).LT.50.or.VarIn(2).lt.1.0/4.0*VarIn(3)) THEN
             VarOut=2.0_dp
          ELSE
             VarOut=1000.0_dp
          END IF

        End FUNCTION VarianceH
!#
       FUNCTION Vdirx(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(5),VarOut
       REAL(kind=dp) :: eps=0.5e+01,Vlim=100.0
       REAL(kind=dp) :: u,norm,a1,a2
       REAL(kind=dp) :: aa,bb
       REAL(kind=dp) :: Arate,rhoi,nflow
       REAL(kind=dp) :: secperyear,grav
 
       secperyear = 60*60*24*365.25
       rhoi       = 917.0
       grav       = 9.81
       nflow      = 3.0
       Arate      = 100.0*1.0e-25

!!          IF (4.0/5.0*VarIn(3).GT.4.0/5.0*Vlim) THEN
!          IF (VarIn(3).GT.Vlim) THEN
!             norm=1._dp
!             u=(+1.0)*abs(VarIn(3))*sin((VarIn(4)+20.0)*3.14159265358979/180.0)
             !norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
!          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             if(norm.gt.0)then
!        !     a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
!        !     a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
!            ! u=-eps*VarIn(1)
!             if(norm.lt.1e-2)then
!             norm=1e-2
!             endif
             u=abs(VarIn(3))*(-1.0)*VarIn(1)
             else
             norm=1._dp
             u=(+1.0)*abs(VarIn(3))*sin((VarIn(4)+0.0)*3.14159265358979/180.0)
             endif
!          END IF

             !if (norm.GT.1.0e-6) Then
       !      VarOut=4.0/5.0*u/norm
             VarOut=u/norm
!            ! VarOut=u/norm
             !else
             !VarOut = 0._dp
             !endif

        End FUNCTION Vdirx
!#
!#
       FUNCTION Vdirx2(Model,nodenumber,VarIn) RESULT(VarOut)

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

!          IF (4.0/5.0*VarIn(3).GT.4.0/5.0*Vlim) THEN
          IF (VarIn(3).GT.Vlim) THEN
             norm=1._dp
             u=(+1.0)*abs(VarIn(3))*sin((VarIn(4)+0.0)*3.14159265358979/180.0)
             !norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
        !     a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
        !     a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
             u=(-1.0)*eps*VarIn(1)
             !if(norm.lt.1e-2)then
             !norm=1._dp
             !endif
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
          END IF

             if (norm.GT.1.0e-6) Then
!             VarOut=4.0/5.0*u/norm
             VarOut=(+1.0)*u/norm
             else
             VarOut = 0._dp
             endif

        End FUNCTION Vdirx2
!#
       FUNCTION Vdiry(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut
       REAL(kind=dp) :: eps=0.5e+01,Vlim=100.0
       REAL(kind=dp) :: u,norm,a1,a2
       REAL(kind=dp) :: aa,bb
       REAL(kind=dp) :: Arate,rhoi,nflow
       REAL(kind=dp) :: secperyear,grav
 
       secperyear = 60*60*24*365.25
       rhoi       = 917.0
       grav       = 9.81
       nflow      = 3.0
       Arate      = 100.0*1.0e-25

!          IF (4.0/5.0*VarIn(3).GT.4.0/5.0*Vlim) THEN
!          IF (VarIn(3).GT.Vlim) THEN
!             norm=1._dp
!             u=(+1.0)*abs(VarIn(3))*cos((VarIn(4)+0.0)*3.14159265358979/180.0)
             !norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
!          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             if(norm.gt.0)then
!        !     a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
!        !     a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
!           !  u=-eps*VarIn(1)
!             if(norm.lt.1e-2)then
!             norm=1e-2
!             endif
             u=abs(VarIn(3))*(-1.0)*VarIn(1)
             else
             norm=1._dp
             u=(+1.0)*abs(VarIn(3))*cos((VarIn(4)+0.0)*3.14159265358979/180.0)
             endif
!          END IF

             !if (norm.GT.1.0e-6) Then
      !       VarOut=4.0/5.0*u/norm
             VarOut=u/norm
!            ! VarOut=u/norm
             !else
             !VarOut = 0._dp
             !endif

        End FUNCTION Vdiry
!#
       FUNCTION Vdiry2(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut
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

!          IF (4.0/5.0*VarIn(3).GT.4.0/5.0*Vlim) THEN
          IF (VarIn(3).GT.Vlim.and.VarIn(1).ne.0) THEN
             norm=1._dp
             u=(+1.0)*abs(VarIn(3))*cos((VarIn(4)+0.0)*3.14159265358979/180.0)
             !norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
          ELSE
             norm=sqrt(VarIn(1)*VarIn(1)+VarIn(2)*VarIn(2))
        !     a1=(norm**nflow)*VarIn(5)**(nflow+1.0)
        !     a2=(2*Arate)/(nflow+1.0)*((rhoi*grav)**nflow)
             u=(-1.0)*eps*VarIn(1)
             !if(norm.lt.1e-2)then
             !norm=1._dp
             !endif
             !u=abs(VarIn(3))*(-1.0)*VarIn(1)
          END IF

             if (norm.GT.1.0e-6) Then
!             VarOut=4.0/5.0*u/norm
             VarOut=(+1.0)*u/norm
             else
             VarOut = 0._dp
             endif

        End FUNCTION Vdiry2
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



