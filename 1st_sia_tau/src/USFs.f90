!#
       FUNCTION largerABSOLUTE(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut
 
       IF (abs(VarIn(1)).GT.abs(VarIn(2)))THEN
          VarOut = abs(VarIn(1))
       ELSE
          VarOut = abs(VarIn(2))
       ENDIF
   
        End FUNCTION largerABSOLUTE
!#
       FUNCTION AtimesC_BtimesC(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(3),VarOut

       VarOut = VarIn(1)*VarIn(3)+VarIn(2)*VarIn(3)

        End FUNCTION AtimesC_BtimesC
!#
       FUNCTION smb(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

       VarOut = VarIn/917

        End FUNCTION smb
!#
       FUNCTION AoverB(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

       VarOut = VarIn(1)/VarIn(2)

        End FUNCTION AoverB

!#
       FUNCTION TenPowerA(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn,VarOut

       VarOut = 10._dp**(VarIn)

       End FUNCTION TenPowerA
!#
       FUNCTION Zb(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

       if (VarIn(2).GT.(VarIn(1)-10.0)) then 
          VarOut=VarIn(1)-10.0 
       else 
          VarOut=VarIn(2)
       endif

       End FUNCTION Zb
!#
       FUNCTION Dalpha(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut

       VarOut = VarIn(1)*(10.0**(VarIn(2)))*log(10.0)

       End FUNCTION Dalpha
!#
!#
       FUNCTION Cost(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut

       VarOut = 0._dp
       if (abs(VarIn(2)).lt.2.0d05) &
            VarOut = 0.5*(VarIn(1)-VarIn(2))*(VarIn(1)-VarIn(2)) 
       if (abs(VarIn(4)).lt.2.0d05) &
       VarOut = VarOut + 0.5*(VarIn(3)-VarIn(4))*(VarIn(3)-VarIn(4))

       End FUNCTION Cost
!#
       FUNCTION DCost(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut
       
       VarOut = 0._dp
       if (abs(VarIn(2)).lt.2.0d05) &
          VarOut = (VarIn(1)-VarIn(2))

       End FUNCTION DCost
!#
       FUNCTION fill_ux(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut
       !-----------------
       ! VarIn(1) - speed
       ! VarIn(2) - azimuth
       ! VarIn(4) - dsdx
       ! VarIn(5) - dsdy

          IF ((VarIn(1).LE.10)) THEN
             VarOut=-VarIn(3)/sqrt(VarIn(3)**2+VarIn(4)**2)*abs(VarIn(1))
          ELSE
             VarOut=abs(VarIn(1))*sin((VarIn(2)+0.0)*3.14159265358979/180.0)
          END IF

        End FUNCTION fill_ux
!#
       FUNCTION fill_uy(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(4),VarOut
       !-----------------
       ! VarIn(1) - speed
       ! VarIn(2) - azimuth
       ! VarIn(3) - dsdx
       ! VarIn(4) - dsdy


          IF ((VarIn(1).LE.10)) THEN
             VarOut=-VarIn(4)/sqrt(VarIn(3)**2+VarIn(4)**2)*abs(VarIn(1))
!             VarOut=abs(VarIn(3))
          ELSE
             VarOut=-abs(VarIn(1))*cos((VarIn(2)+0.0)*3.14159265358979/180.0)
          END IF

        End FUNCTION fill_uy
!#
       FUNCTION thickness_limit(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(1),VarOut
       !-----------------
       ! VarIn(1) - thi

          IF ((VarIn(1).LE.5)) THEN
             VarOut=5.0
          ELSE
             VarOut=VarIn(1)
          END IF

        End FUNCTION thickness_limit
!#
       FUNCTION calc_ux(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut
       !-----------------
       ! VarIn(1) - speed
       ! VarIn(2) - azimuth
       ! VarIn(3) - angle shift in azimuth field

       VarOut = abs(VarIn(1))*sin((VarIn(2)+0.0)*3.14159265358979/180.0)

       End FUNCTION calc_ux
!#
       FUNCTION calc_uy(Model,nodenumber,VarIn) RESULT(VarOut)

       USE DefUtils

       implicit none
       !-----------------
       TYPE(Model_t) :: Model
       INTEGER :: nodenumber
       REAL(kind=dp) :: VarIn(2),VarOut
       !-----------------
       ! VarIn(1) - speed
       ! VarIn(2) - azimuth
       ! VarIn(3) - angle shift in azimuth field

       VarOut = -abs(VarIn(1))*cos((VarIn(2)+0.0)*3.14159265358979/180.0)

       End FUNCTION calc_uy
!#
