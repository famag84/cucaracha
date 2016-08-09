!--------------------------------------------------------------------
!        Calcula Energia libre
! 2015-04-20 F.M.G. : Energía de Pong super chequeada!
!--------------------------------------------------------------------

subroutine calc_energy(pHbulk)
#   include "control_run.h"
    use globales, only: delta, vsol, vsalt, radio, dimR , pi, lb, zwall
    use pore, only: qtot, psi, fdiswall
    use csys, only: eps1, log_q, sigma, sigmaq!, constq
    use FreeEnergy
    use mpmodule
    implicit none
    real(kind=8), intent(in) :: pHbulk
    real(kind=8) :: suma_pong, std_mupol!, aux_mp
!    real(kind=8) :: F_Mix_pos, F_Mix_neg, F_Mix_Hplus, F_Mix_OHmin, &
!                    F_Conf, F_Eq, F_Eq_wall, F_vdW, F_electro, F_eps
    integer :: iR,j

    print*, "Calculating energies..." 
    Free_Energy = 0.0
    Free_Energy2 = 0.0
    std_mupol= 0.0

! ******************************************************************
! ENTROPIAS DE MEZCLA (independiente de la presencia del polimero!)
! ******************************************************************
    F_Mix = 0.0
    F_Mix = fmix() ! Calcula las entropías de mezcla por separado y devuelve el total
                   ! (cada contribucion es guardada en una variable individual)
    Free_Energy =  Free_Energy + F_Mix

! ********************* POLYMER ************************************
! 6. Entropia interna polimero + Eq. Quimico (cambiar?)
! ******************************************************************
    F_Conf = 0.0
# if CHAIN != 0 
    F_Conf = fconf_pol() *(sigma*delta/vsol)
# endif
    Free_Energy = Free_Energy + F_Conf
!    print*, "E + F_Fconf" , Free_energy
 
    !std_mupol = std_mupol + F_conf ! just the sum: P*log(P) (viejo)
    std_mupol = -log_q

! 7. Chemical Equilibria
! ******************************************************************
    F_Eq = 0.0 
    F_Eq_wall = 0.0
# if CHAIN != 0 && POL != 2 
! Funciona bien
      F_Eq = fchem_eq() ! Solo calcula el ChemEq. en el polimero
# endif
# if fsigmaq == 1
      F_Eq_wall = fchem_eq_wall()
# endif
      Free_Energy = Free_Energy + F_Eq + F_Eq_wall
!    print*, "E + F_Eq" , Free_energy

! 8.vdW ! Ojo, los  son negativos => atraccion
    F_vdW = 0.0
# ifdef VDW
    F_vdW = fvdW()
    Free_Energy = Free_Energy + F_vdW
    print*, "E + F_vdW" , Free_energy
# else
    print*, "F_vdW = 0.0 "
# endif

! 9. Electrostatic - Esta calculada en las unidades correctas
F_electro = 0.0    
do iR  = 1, dimR
    F_electro = F_electro &
                + psi(iR)*(qtot(iR))/2.0 *(delta)*(dfloat(iR)-0.5)*delta/Radio
enddo
! La siguiente opcion no tiene el mismo sigmaq 
! que esta en la condicion de borde de psi(dimR+1)
    F_electro = F_electro + sigmaq*(delta/vsol)*zwall*fdiswall*psi(dimR)/2.0 
! La siguiente opcion es coherente con la definicion de psi(dimR+1)
!    F_electro = F_electro + sigmaq*(delta/vsol)*constq*zwall*fdiswall*psi(dimR)/2.0 
    
    Free_Energy = Free_Energy + F_electro
!    print*, "E + F_electro" , Free_energy
!*******************************************************************
! Observacion:
! La energia correspondiente a la distribucion de carga superficial 
! sigmaq es tenida en cuenta en ambas energias Pong y acá.
! (ver pong_energy.f90) 
!*******************************************************************

! 10. Pol-Sup
    F_eps = 0.0
#if CHAIN != 0
    F_eps = fpol_sup() 
!    if(eps1.ne.0.0) then
!        print*, 'EPS should be 0 or correct free energy!'
!        stop
!    endif
#endif
    Free_Energy = Free_Energy + F_eps
 ! 11. Osmotic-Pressure!? should be zero. this is a check of the packing constraint
    F_ospi =0.0 
!    F_ospi = fospi()
!    F_ospi = F_ospi * delta
!    Free_Energy = Free_Energy - F_ospi
!    print*, "Energia osmotic pressure: " , F_ospi 

    print*, "E " , Free_energy
!*******************************************************************
! Observacion:
! La energia correspondiente a la distribucion de carga superficial sigmaq 
! es tenida en cuenta dentro de pong_energy. Esto es asi por que 
! la integral corespondiente a la energia electrostatica en Free_energy 
! llega hasta R-delta/2
!*******************************************************************
    Free_Energy2 = 0.0
    suma_pong = pong_energy() ! pong_energy(): considera la carga superficial sigmaq
    Free_Energy2 = suma_pong - (delta/vsol)*sigma*log_q - F_vdW !&
!*********************************************************************************
! La siguiente definicion no coincide con la condicion de contorno psi(dimR+1)
!    Free_Energy2 = Free_Energy2 + sigmaq*zwall*fdiswall*(delta/vsol)*psi(dimR)/2.0 & 
    
    Free_Energy2 = Free_Energy2 + sigmaq*zwall*fdiswall*(delta/vsol)*psi(dimR)/2.0 &
                                + F_Eq_wall
!
! Aunque la expresion para la carga superficial aparece en todas mis deducciones 
! aparentemente no hace falta. Es un termino constante. 
! Si se considera, debe ser considerada en las dos energias al mismo tiempo.
! (ver 9. Electrostatic en calc_energy)
! Expresion para la energia en la superficie:
!*********************************************************************************

! Printing output
    Print*, "E2 ", Free_Energy2
    Print*, "diff", Free_Energy -  Free_Energy2

! Guarda energia libre
         write(311,*) sigma*delta/vsol, F_electro
         write(312,*) sigma*delta/vsol, Free_energy2, suma_pong
         write(301,*) sigma*delta/vsol, Free_energy
         write(302,*) sigma*delta/vsol, F_Mix_s 
         write(303,*) sigma*delta/vsol, F_Mix_pos
         write(304,*) sigma*delta/vsol, F_Mix_neg
         write(305,*) sigma*delta/vsol, F_Mix_Hplus
         write(306,*) sigma*delta/vsol, F_Mix_OHmin
         write(308,*) sigma*delta/vsol, F_Eq
         write(319,*) sigma*delta/vsol, F_Eq_wall
         write(202,*) sigma*delta/vsol, std_mupol
#ifdef VDW
         write(309,*) sigma*delta/vsol, F_vdW
#endif
# if CHAIN != 0
         write(307,*) sigma*delta/vsol, F_Conf
         write(310,*) sigma*delta/vsol, F_eps
# endif
         write(201,*) sigma*delta/vsol, F_ospi

!c--------------------------- FIN DE ENERGIA LIBRE -----------------

end subroutine calc_energy
