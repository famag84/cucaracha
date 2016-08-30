function fchem_eq()
#   include "control_run.h"
    use globales, only: dimR, delta, vpol, vsol, radio
    use csys, only: expmuHplus, expmuOHmin, expmuneg, K_Cl0,Ka0, Kb0
!    use FreeEnergy, only: checknumber
# if POL == 0 /* PAH */
#ifdef PAHCL
    use pore, only: avpol, fdis , fdis2
#else
    use pore, only: avpol, fdis !, fdis2
#endif
# elif POL == 2 /* Neutral Polymer */
    use pore, only: avpol, fdis !, fdis2
# elif POL == 1 /* PMEP */
    use pore, only: avpol, fdis, fdis2
# endif
    implicit none
    real(kind=8) :: fchem_eq
    integer :: i
    
!    print*,  "WARNING! fchem_eq Calc. Be sure that neither f's is equal to zero!"
!    print*,  "Check fdis(), fdis2() and 1-f()-fdis2()" 
!    call printstate("fchem L14") ! Report of State
    
    fchem_eq = 0.0
    
    do i = 1, dimR
# if POL == 0 /* PAH */
#ifdef PAHCL
! Nueva expresion para el equilibrio quimico PAH + CL
      fchem_eq = fchem_eq + ( fdis(i) *dlog( fdis(i) ) &
                          +   fdis2(i) *dlog( fdis2(i)*Ka0  ) &
                          + (1-fdis(i)-fdis2(i))*dlog( (1-fdis(i)-fdis2(i) )*K_Cl0 ) )&
                          *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio
      fchem_eq = fchem_eq + fdis(i)*dlog(expmuneg) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio ! ojo! 
      fchem_eq = fchem_eq + fdis2(i)*(dlog(expmuneg)-dlog(expmuOHmin)) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio ! 
!! Mantiene la coherencia, expmuneg cambia con el pH! por lo tanto es una constante valiosa.
      fchem_eq = fchem_eq - dlog(expmuneg) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio !!
#else
! Nueva expresion para el equilibrio quimico
!      fchem_eq = fchem_eq + ( fdis(i) *dlog( fdis(i)/Ka0 ) &
!                          + (1-fdis(i))*dlog( (1-fdis(i)) ) )  *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio
!      fchem_eq = fchem_eq + fdis(i)*dlog(expmuOHmin) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio ! ojo! expmuHplus ya tiene el signo menos!!




! Propuesta Mario
      fchem_eq = fchem_eq + ( fdis(i) *dlog( fdis(i)) &
                          + (1-fdis(i))*dlog( (1-fdis(i))*Ka0 ) )  *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio
      fchem_eq = fchem_eq - (1-fdis(i))*dlog(expmuOHmin) *(avpol(i)/(vpol*vsol)) *delta*(dfloat(i)-0.5)*delta/Radio ! ojo! expmuHplus ya tiene el signo menos!!

! Original Mario
!      fchem_eq = fchem_eq + fdis(i)*dlog(fdis(i))            *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!      fchem_eq = fchem_eq + (1.0-fdis(i))*dlog(1.0-fdis(i))  *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!                   
!      fchem_eq = fchem_eq + (1.0-fdis(i))*dlog(Ka0)            *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
!      fchem_eq = fchem_eq + (1.0-fdis(i))*(-dlog(expmuOHmin)) *avpol(i)/vpol *(dfloat(i)-0.5)*delta/Radio
#endif

# elif POL == 1
! Here PMEP - TO CHECK!
      fchem_eq = fchem_eq + ( fdis(i)*dlog(fdis(i)/Ka0) &
                            + fdis2(i)*dlog(fdis2(i)/(Ka0*Kb0)) &
                            + (1.0-fdis(i)-fdis2(i))*dlog(1.0-fdis(i)-fdis2(i)) &
                            ) *(avpol(i)/(vpol*vsol)) *delta *(dfloat(i)-0.5)*delta/Radio
!                            ) *avpol(i)/vpol*(dfloat(i)-0.5)*delta/Radio

      fchem_eq = fchem_eq + ( fdis(i)*(dlog(expmuHplus)) & ! ojo! expmuHplus ya tiene el signo menos!!
                            + 2*fdis2(i)*(dlog(expmuHplus)) &
                            ) *(avpol(i)/(vpol*vsol)) *delta *(dfloat(i)-0.5)*delta/Radio
!                            ) *avpol(i)/vpol*(dfloat(i)-0.5)*delta/Radio

# endif

    enddo

!    print*, "fchem_eq llama checknumber: fchem_eq", fchem_eq
!    call checknumber(fchem_eq,'fchem_eq')
    
    return
    contains 
!        function checknumber(var, arg) result(bool)
!            implicit none
!            real(kind=8), intent(in) :: var
!            character(len=*), intent(in) :: arg
!            logical :: bool 
!            ! Check if Not a Number
!            if ( var /= var ) then
!                print*, arg, " real number is NaN"
!                call printstate('checknumber real NaN')
!                bool=.true. 
!            endif
!        
!            if ( var-1 == var ) then
!                print*, arg, " real number is infinity"
!                call printstate('checknumber real infinity')
!                bool=.true.
!            endif
!        end function checknumber
!        subroutine checknumber(var, arg)
!            implicit none
!            real(kind=8), intent(in) :: var
!            character(len=*), intent(in) :: arg
!            ! Check if Not a Number
!            if ( var /= var ) then
!                print*, arg, " real number is NaN"
!                call printstate('checknumber real NaN')
!                stop
!            endif
!        
!            if ( var-1 == var ) then
!                print*, arg, " real number is infinity"
!                call printstate('checknumber real infinity')
!                stop
!            endif
!        endsubroutine checknumber
end function 
