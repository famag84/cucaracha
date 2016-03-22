subroutine set_pore_distrib
! To control the code to be included in the compilation
#   include "control_run.h" 
use mpmodule
use globales
use csys
use pore
implicit none
interface    
    function factorcurv(i)
        integer, intent(in) :: i
        double precision :: factorcurv
    end function factorcurv
end interface 
integer :: iR, aR, i,j 
real(kind=8) :: shift, denon ! it is  multiplicative factor in pro(cuantas)
type (mp_real) q ! este valor puede ser MUY grande
!real(kind=8) :: aux=0.0

# ifdef fdebug    
    call printstate("set_pore_distrib L19")
# endif
!**************************************************************
! Volume fractions and degrees of dissociations: fdis
do iR=1,dimR
! Mis expresiones:
! expmu**** = son las exponenciales con los pot. quimicos calculadas 
!             con concentraciones de BULK,
! * xh : concentracion de solvente, tiene una expresion en terminos de pi(r).
!   xpos(iR) = expmupos*(xh(iR)**vsalt) *dexp(-psi(iR)*zpos) ! ion plus vol. fraction
   xpos(iR) = vsalt*expmupos*(xh(iR)**vsalt) *dexp(-psi(iR)) ! ion plus vol. fraction

!   xneg(iR) = expmuneg*(xh(iR)**vsalt) *dexp(-psi(iR)*zneg) ! ion neg vol. fraction
   xneg(iR) = vsalt*expmuneg*(xh(iR)**vsalt) *dexp(psi(iR))  ! ion neg vol. fraction

!  xHplus(iR) = expmuHplus*(xh(iR)**vHplus) *dexp(-psi(iR)*zH) ! H+ vol. fraction
  xHplus(iR) = expmuHplus*xh(iR) *dexp(-psi(iR))               ! H+ vol. fraction

!  xOHmin(iR) = expmuOHmin*(xh(iR)**vOHmin) *dexp(-psi(iR)*zOH) ! OH- volume fraction
  xOHmin(iR) = expmuOHmin*xh(iR) *dexp(psi(iR))           ! OH-  volume fraction

!  fdis(iR) = dissos_degree(1,iR)
# if CHAIN != 0 
#   if POL == 0 /* PAH */
!**************************************************************
! New symmetric equations!
! Facundo like expression
 !fdis(iR) = 1.0 / (1.0 + expmuOHmin*dexp(psi(iR)*zpol)/Ka0 ) ! FACUNDO funciona(!)
 !fdis(iR) = Ka0 / (expmuOHmin*dexp( psi(iR)*zpos) + Ka0 ) ! using fdiswall symmetry 

fdis(iR) = Ka0 / ( xOHmin(iR)/xh(iR)  + Ka0   ) ! using fdiswall symmetry 

! fdis(iR) = 1.0 / (1.0 + xOHmin(iR)/(xh(iR)*Ka0) ) ! MARIO EXPRESSION ! Funciona 
! fdis(iR) = 0.9
! fdis(iR) = 1.0d0/(1.0d0 + (Ka0*dexp(psi(iR)*zpos)/expmuHplus))! from the theory
#   elif POL == 1 /* PMEP */
! Para estas expresiones es necesario definir las constantes K1 y K2

     denon = ( expmuHplus*dexp(psi(iR)*zneg) )**2 &
             + expmuHplus*dexp(psi(iR)*zneg) *Ka0 &
             + Ka0*Kb0 ! Kb0 = Ka20 for dobule protonation
    fdis(iR) = expmuHplus*dexp(psi(iR)*zneg)*Ka0 / denon 
    fdis2(iR) = Ka0*Kb0 / denon ! Kb0 = Ka20 for dobule protonation

!!     denon = ( expmuOHmin*dexp(psi(iR)*zneg) )**2 &
!!             + expmuOHmin*dexp(psi(iR)*zneg) *Ka0 &
!!             + Ka0*Kb0 ! Kb0 = Ka20 for dobule protonation
!!    fdis(iR) = expmuOHmin*dexp(psi(iR)*zneg)*Ka0 / denon 
!!    fdis2(iR) = Ka0*Kb0 / denon ! Kb0 = Ka20 for dobule protonation
!   fdis2(iR) = 1.0d0/(1.0d0 + (dexp(psi(iR)*zpos)*expmuOHmin/Kb0) )
# elif POL == 2 /*Neutral Polymer*/
    fdis(iR) = 0.0 ! Neutral polymer
#   endif /* POL */
# endif /* CHAIN */
end do
!**************************************************************
! Dissociation in the inner wall of the pore 
! + interaction with the electrostatic Potential    
!    fdiswall = Kwall0 / (expmuHplus*dexp(psi(dimR)*zwall) + Kwall0) 
!    fdiswall = 1.0 / (expmuHplus*dexp(psi(dimR)*zwall)/Kwall0 + 1.0 ) 
# if fsigmaq == 1 /*  Si hay regulacion en la pared */
    fdiswall = 1.0 / ( xHplus(dimR)/xh(dimR) /Kwall0 + 1.0 ) 
# elif fsigmaq == 0 /* No hay regulacion en la pared */
    fdiswall = 1.0  
# endif
!**************************************************************

!**************************************************************
! Boundary Conditions: Electrical Potential
!! Estas ecuaciones son suplementarias a la eq. de Poisson discretizada
!   *   The derivate at r=0 should be zero. 
!       Se utiliza orden 1 en la derivada en r=0.
     psi(0) = psi(1) ! La derivada en r = 0 es igual a cero
!   *   The derivate at r=R takes into account the superficial charge: sigmaq. 
!       Se utiliza orden 2 en la derivada.
     psi(dimR+1) = psi(dimR) + delta*(constq)*(sigmaq*delta/vsol)*zwall*fdiswall 
! factor in poisson eq:          constq=(4*pi*lb) 

! working on: Boundary condition of a dipole layer (elefante)
!     psi(dimR+1) = psi(dimR) - (2*sigmaq/constq) 
! La derivada en r=R es el salto de la carga superficial, ver notas
!**************************************************************

# if CHAIN != 0
!!!!!! AQUI FALTA ACTUALIZAR xH! debe tomar el valor de  x1
do iR = 1, dimR
! (xh(iR)**vpol): viene de reemplazar la presion osmotica por 
!                 la expresion para el solvent
#   if POL == 0 /* PAH */
    ! Para polimero con regulacion de carga (cargado positivamente) elefante mayor
    xpot(iR) = (xh(iR)**vpol) / (1-fdis(iR))
    !xpot(iR) = (xh(iR)**vpol)*exp(-psi(iR)*zpol) /(1-fdis(iR)-fdis2(iR) )
!    xpot(iR) = xpot(iR) *exp(-psi(iR)*zpol) *expmuOHmin *Ka0*(1.0-fdis(iR))/fdis(iR)
#   elif POL == 1 /* PMEP */
    xpot(iR) = (xh(iR)**vpol)/(1.0-fdis(iR)-fdis2(iR) ) 
#   elif POL == 2 /* Neutral Polymer */
    xpot(iR) = (xh(iR)**vpol) ! Polimero neutro ! For Neutral Polymers OK!
#   endif /* PAH || PMEP */
# ifdef fdebug_set_pore_distrib
print*, "set_pore_distrib L115 iR, xpot(iR), fdis(iR), xh(iR): ", iR, xpot(iR), fdis(iR), xh(iR)
#endif
# ifdef VDW
!**************************************************************
! To calculate VdW force coefficient
!**************************************************************
! Taking in to account the probability with the vdW force coefficient for every layer
   do j = 1, dimR
!   xpot(iR) = xpot(iR) * dexp(vsol*vpol*Xu(iR,j)*xtotal(j+1))!Porque j+1? (elefante)
       xpot(iR) = xpot(iR) * dexp(vsol*vpol* Xu(iR,j)*xtotal(j))
   end do
#endif /* VDW */
enddo
# ifdef fdebug_set_pore_distrib
        print*, "set_pore_distrib L133: loop end"
#endif

!      q = 0.0 ! Normalization of P(alpha)
      q = '0.0' ! Normalization of P(alpha) remember that: type (mp_real) q !!
# if POL == 2 /* Neutral Polymers */
  shift = 1.0 ! (Should be "shift" an input parameter in fort.8?) 
# else
  shift = 1.0d-50 ! (Should be "shift" an input parameter in fort.8?) 
# endif
avpol(:)= 0.0 ! line important to probability calculus

!**************************************************************
! Calculo de la densidad de probabilidad de cada configuracion
! Estos Do's estan bien, siempre dejar la coma a la izquierda.
do i=1,cuantas ! i enumerate configurations (configurations ensamble)
    pro(i)=1.0*shift
!      do j=1,long, 2 ! (long=28) ! Here you choose the type of the first segment
    do j=1,long ! (long=28)
        aR = pR(i, j) 
!            pR()'s output is the layer where is the segment j of configuration i.
!            bR = pR(i, j+1)   
! The configuration's probability is the product of layer's probability
! many times as particles in that layer.
            pro(i) = pro(i) * xpot(aR)
! ATENTTION: Now you insert structural detail about alternating type of 
! monomers, A - B - A - B - etc
!            pro(i) = pro(i) * xpot_neg(aR) * xpot_pos(bR)
# ifdef fdebug_set_pore_distrib
        print*, "set_pore_distrib L154, i, j, aR=pR(i,j), pro(i): ", i, j, aR, pro(i)
#endif
    enddo
! q es la suma de todas las probabilidades
    q=q+pro(i)
    do j = 1,long
        aR = pR(i, j) 
! pR devuelve en que layer se encuentra el monómero j del polímero i.
     ! OJO aca no es sigma el factor multiplicativo! elefante!
       avpol(aR) = avpol(aR) + pro(i)*sigma*vpol*factorcurv(aR) ! cilindro, ver notas
! ver eq:factorcurv in mis_apuntes.lyx
!       bR = pR(i, j+1)
!       avpol(aR) = avpol(aR) + pro(i)*sigma*vpol ! plano
! sigma*vpol*Factorcurv(aR) es el factor en la densidad de polimero.
! cilindro, ver notas
!      avpoln(aR) = avpoln(aR) + pro(i)*sigma*vpol*Factorcurv(aR) 
!      avpolp(bR) = avpolp(bR) + pro(i)*sigma*vpol*Factorcurv(bR)
# ifdef fdebug_set_pore_distrib
        print*, "set_pore_distrib L172, aR=pR(i,j), avpol(aR): ", aR, avpol(aR)
#endif
    end do
enddo ! End loop over chains/configurations
# ifdef fdebug_set_pore_distrib
        print*, "set_pore_distrib L177: End Loop"
#endif

    call mpwrite(11,q/shift)
    log_q = log(q/shift) ! Variable clave en el calculo de energias! MPLOG
! log(q) is nepperian log of q. has 15 digits of precision.
    
    do i=1,cuantas 
        pro(i) = pro(i)/q
    enddo
!    print*, "pro(:): ", pro(:)
!    call printstate("L105 set_pore_distrib")
    do iR=1, dimR            ! norma avpol
        avpol(iR) = avpol(iR)/ q
!        avpoln(iR) = avpoln(iR) / q
!        avpolp(iR) = avpolp(iR) / q 
!         avpol(iR) = avpoln(iR) + avpolp(iR)
    enddo
! Chequeo avpol
!!    do iR=1,dimR
!!        temp2= temp2+avpol(iR)
!!    enddo
!!    print*, "suma avpol", sigma, temp2

! Chequeo q
!    call checknum(q,'q en set_pore_distrib')
!    call printstate("L115 set_pore_distrib")
# endif

! Local charge
    do iR=1,dimR  

        qtot(iR) = (zpos*xpos(iR)+zneg*xneg(iR))/(vsalt*vsol) &
!                  + zpos2*xpos2(iR)/vsalt2 &
# if CHAIN == 1
#   if POL == 0 /* PAH */
            + (fdis(iR) )*zpol*avpol(iR)/(vpol*vsol) &
#   elif POL == 1 /* PMEP */
            + (fdis(iR) + 2.0*fdis2(iR))*zpol*avpol(iR)/(vpol*vsol) &
!            + (avpolp(iR)*fdis2(iR)*zpos + avpoln(iR)*fdis(iR)*zneg)/vpol
#   elif POL == 2 /* Neutral Polymers */
! Neutral Polymers do not add local charge!
#   endif /* PAH || PMEP */
# endif
            + xHplus(iR)/vsol -xOHmin(iR)/vsol 
    enddo

end subroutine set_pore_distrib
