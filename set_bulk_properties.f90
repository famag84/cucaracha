! Here all variables related with bulk

subroutine set_bulk_properties(pHbulk, csaltbulk, cpolbulk)
#   include "control_run.h"
    use globales
    use csys
    implicit none
    real(KIND=8), intent(in) :: pHbulk, csaltbulk, cpolbulk

! Volume fraction salt in mol/l 
    xsalt=(csaltbulk*Na/(1.0d24))*(vsalt*vsol)   

cHplus = 10**(-pHbulk) ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
! Relacion que sale del Atkins íntrinseco a las definiciones de equilibrio con pH
pOHbulk= pKw -pHbulk 
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
!------- Volume fractions
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
if(pHbulk.le.7) then  ! pH<= 7
    xposbulk= xsalt/zpos
 !xnegbulk=-xsalt/zneg -xsalt2/zneg +(xHplusbulk -xOHminbulk)*vsalt ! NaCl+ HCl
    xnegbulk=-xsalt/zneg +(xHplusbulk -xOHminbulk)*vsalt ! NaCl+ HCl
else                  ! pH >7
    xposbulk= xsalt/zpos +(xOHminbulk -xHplusbulk)*vsalt ! NaCl+ NaOH
    !xnegbulk=-xsalt/zneg -xsalt2/zneg
    xnegbulk=-xsalt/zneg
endif

# if CHAIN == 1 && MUPOL == 1
! Volume fraction pol
    xpolbulk=(cpolbulk)*(Na/1.0d24)*(vpol*vsol)
    xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk ! - xpolbulk
# else
! En bulk no hay polimero.
         xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk !- xposbulk2

#endif

! Ojo Kw y Ka estan en mol/l mientras que vsol y xsolbulk estan en nm3
! 1.0d24 factor de conversion de volumen de litro a nm3  1nm3 = 10^-24 l = 10^-24 dm3

! remember that vsalt, vHplus and vOHmin are in vsol unit!
           expmupos = xposbulk/vsalt / xsolbulk**vsalt 
           expmuneg = xnegbulk/vsalt / xsolbulk**vsalt 
!ojo! expmuHplus ya tiene el signo menos!!
         expmuHplus = xHplusbulk / xsolbulk**vHplus   ! vsol = vHplus
         expmuOHmin = xOHminbulk / xsolbulk**vOHmin  ! vsol = vOHmin

! intrinstic equilibruim constants in mol/nm3
            Ka0 = (Ka*vsol/xsolbulk)*(Na/1.0d24)
            Kb0 = (Kb*vsol/xsolbulk)*(Na/1.0d24)
         Kwall0 = (Kwall*vsol/xsolbulk)*(Na/1.0d24)
#ifdef PAHCL
          K_Cl0 = (K_Cl*vsol/xsolbulk)*(Na/1.0d24)
#endif 

! Polymer Bulk Properties
# if CHAIN != 0 
#   if POL == 0 /* PAH */
#       ifdef PAHCL
#       else
        fdisbulk = Ka0 / (( xOHminbulk/xsolbulk)  + Ka0   ) ! PAH using fdiswall symmetry
#       endif
#   elif POL == 1 /* PMEP */
#   elif POL == 2 /*Neutral Polymer*/
#   endif /* POL */
# endif /* CHAIN */

!print*, fdisbulk, ":fdisbulk"
!Q_bulk = cuantas * xsolbulk**(long*vpol) /(1.0-fdisbulk)**long
!print*, "Q_bulk , cuantas *xsolbulk**(long*vpol), (1-fdisbulk)**long, log(1-fdisbulk)", &
!         Q_bulk , cuantas *xsolbulk**(long*vpol), (1-fdisbulk)**long,long*log(1.0-fdisbulk)
!expmupol = xpolbulk/vpol / Q_bulk

expmupol = xpolbulk/(long*vpol) / xsolbulk**(long*vpol) /cuantas * (fdisbulk)**long

! std_mupol : bulk equation
  std_mupol = dlog(xpolbulk/(long*vpol) ) - log(1.0*cuantas) - vpol*long*dlog(xsolbulk)+long*dlog(fdisbulk)

!print*, "std_mupol, expmupol: ", std_mupol, expmupol

end subroutine set_bulk_properties
