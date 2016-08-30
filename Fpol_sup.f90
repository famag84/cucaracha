function fpol_sup()
    use globales, only: delta, radio, dimR, vpol, vsol, long
    use csys, only: eps1, eps
    use pore, only: avpol
!    use FreeEnergy, only: checknumber
    implicit none
    real(kind=8) :: fpol_sup
    integer :: iR
!  2015-04-20 Interaccion polimero superficie (electrostatica(?)<-nop it could be chemist or entropic interaction)
    fpol_sup=0.0

!    do iR = dimR-1, dimR
!    iR=dimR
            !fpol_sup = fpol_sup - eps1*avpol(iR)*delta*(dfloat(iR)-0.5)*delta/radio
!            fpol_sup = fpol_sup - eps(iR)*avpol(iR)*delta*(dfloat(iR)-0.5)*delta/radio
!    enddo
        ! Check avpol(iR) expression, here we need rho_pol
! Por compatibilidad con la discretización es necesario poner el jacobiano
        fpol_sup = eps(dimR)*avpol(dimR)/vpol/vsol*delta*(dfloat(dimR)-0.5)*delta/radio

    print*, "fpol_sup L16: " , fpol_sup 
!    print*, "f_vdW llama checknumber: f_vdW", f_vdW
    call checknumber(fpol_sup,'fpol_sup')
    
    return
    contains 
        subroutine checknumber(var, arg)
            implicit none
            real(kind=8), intent(in) :: var
            character(len=*), intent(in) :: arg
            ! Check if Not a Number
            if ( var /= var ) then
                print*, arg, " real number is NaN"
                call printstate('checknumber real NaN')
                stop
            endif
        
            if ( var-1 == var ) then
                print*, arg, " real number is infinity"
                call printstate('checknumber real infinity')
                stop
            endif
        endsubroutine checknumber
end function 
