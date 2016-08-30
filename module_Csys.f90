! ###
!
! El objetivo es poner aca todas las variables del sistema fisicoquimico
! NO poner variables del solver
! NO poner constantes fisicas o geometricas (ver globales)
!
! ##

module csys
!  use mpmodule
    logical :: chainpol
    integer, parameter :: sp = kind(1.0),    &
                          dp = selected_real_kind(2*sp),   &
                          qp = selected_real_kind(2*dp)   
! Variables Fisicoquimicas
!    use globales !  implicit none esta en globales

! System Variables inputs
    real(kind=dp) :: Q_bulk
    real(kind=8) :: sigma, sigmaq, constq, &
        cpol, xpolbulk, fdisbulk, &
        csalt, csalt2=0.0, xsalt, xsalt2,& ! Salt concentration
        pKawall, kwall, Kwall0, pKa, Ka, Ka0, pKb, Kb, Kb0, pK_Cl, K_Cl, K_Cl0,& ! chemical equilibrium
        eps1, & ! Electrostatic potential
        longporo, &
        movpos, movneg, movHplus, movOHmin,&  ! Mobilities
        st
! Bulk values
    real(kind=8) :: cHplus, cOHmin, &
                    pOHbulk, xposbulk, xnegbulk, xsolbulk, &
                    xHplusbulk, xOHminbulk, expmupol,&
                    std_mupol, expmupos, expmuneg, expmuHplus, expmuOHmin
    real(kind=8), dimension(:), allocatable :: eps
    real(kind=8), dimension(:), allocatable :: pHs, vsigma, vcsalt, vcpol ! list of bulk pHs
    real(kind=8), dimension(:,:), allocatable :: Xu
    real(kind=8), dimension(:,:,:), allocatable :: in1 ! guarda las configuraciones de cadena 
    real(kind=8), dimension(:), allocatable :: x1, xg1, xflag
    real(kind=dp),dimension(:), allocatable ::  pro !  list of probability 
!    type(mp_real),dimension(:), allocatable ::  pro !  list of probability 
    integer, dimension(:,:), allocatable :: pR ! pR stores the position of the segment j from conf. i . 
!    integer, dimension(:,:)dimension(cuantas,long)  :: pR ! pR stores the position of the segment j from conf. i . 
!    real(kind=8), dimension(cuantas,long,3) :: in1 ! guarda las configuraciones de cadena 
!    real(kind=8), dimension(dimR,dimR) :: Xu
    real(kind=16) :: log_q=0 !
!    real(kind=8) :: norma
    integer :: infile, iter, npH, ncpol, ncsalt, nsigma, seed, npKa=2 ! ojo! seed cambia a lo largo del programa!
    integer :: ps_i, flag ! used in check_run subroutine
    integer :: ipH, icsalt, icpol, isigma ! Contador en el main
    ! ps_i : printstate counter
end module csys
