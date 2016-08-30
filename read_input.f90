subroutine read_input
#   include "control_run.h"
    use globales, only: cuantas, dimR, long, radio, delta, neq
    use csys
    use pore, only: shift_f
    implicit none
    integer :: i
    character :: basura

    open(8,status='old', action='read' )
    print*, 'Reading input parameters'
! Input variables from fort.8
      read(8, *), basura
      read(8, *), nsigma   ! Surface coverage
! OJO! esta memoria debe ser liberada en el main! becareful!
      allocate(vsigma(nsigma)) 
      read(8, *), (vsigma(i), i=1, nsigma) ![#chains/nm^2]

      read(8, *), basura
      read(8, *), sigmaq   ! Surface chargeable sites [e/nm^2]

      read(8, *), basura
      read(8, *), eps1     ! Interaccion Pol-Sup  [J?]

! Salt concentration
      read(8, *), basura
      read(8, *), ncsalt   ! Salt Concentration [M]
      allocate(vcsalt(ncsalt)) 
      read(8, *), (vcsalt(i), i=1, ncsalt)
        csalt=vcsalt(1) ! Me quedo con la primeer sal!
! polymer concentration
      read(8, *), basura
      read(8, *), ncpol   ! Monomer of pol Concentration
      allocate(vcpol(ncpol)) 
      read(8, *), (vcpol(i), i=1, ncpol)

      read(8, *), basura
      read(8, *), pKawall ! Eq. Const. wall charge

      read(8, *), basura
      read(8, *), npKa
      read(8, *), pKa, pKb    ! Eq. Const. Polymer
     
      read(8, *), basura
#ifdef PAHCL
      read(8, *), pK_Cl   ! Eq. Const. Polymer
#else
      read(8, *), basura  ! Eq. Const. Polymer
#endif
 
      read(8, *), basura
      read(8, *), npH
! OJO! esta memoria debe ser liberada en el main! becareful!
      allocate(pHs(npH)) 
      read(8, *), (pHs(i), i=1, npH)

      read(8, *), basura
      read(8, *), infile
! st parameter is for vdW interaction on/off ?
      read(8, *), basura
!      read(8, *), st
      read(8, *), cuantas, dimR, long, shift_f  ! radio size
        radio=dimR*delta ! Radio in [nm], dimR = Diameter in [nm] = number of layers
        neq=2*dimR

      read(8, *), basura
      read(8, *), movpos, movneg, movHplus, movOHmin, longporo ! [(m^2)*mS/M], longporo in micrometers

      read(8, *), basura
      read(8, *), seed
        print*, seed, "<- Seed"
    close(unit=8)

# if CHAIN != 0
    print*, "Program compiled to include chains in the system!"
    if (cuantas*sigma*long <= 0.0 ) then
        print*, "No polymer chains in this run."
        chainpol = .False.
    endif
#endif

end subroutine 
