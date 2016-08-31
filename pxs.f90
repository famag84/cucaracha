!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

subroutine pxs
#   include "control_run.h"
!  Esta subrutina se encarga de chequear que  todos los segmentos esten dentro del slab
    use globales
    use csys
    implicit none

    REAL(KIND=8) :: vect
    INTEGER :: i,j
#ifdef fdebug_pxs
    print* , "long, cuantas, chaintot: ", long, cuantas, chaintot
#endif
    chaintot=0
    do i=1,cuantas
        chaintot =chaintot+1
        do j=1,long ! 1 = y, 2 = x, 3 = z... z no se usa...)
            vect = sqrt((in1(i,j,1))**2 + in1(i,j,2)**2) ! distancia del centro al segmento
# if MUPOL == 1
!            print*, i,j,vect, (int(vect/delta)+1), dimR
            if( (int(vect/delta)+1).gt.dimR) then
                chaintot=chaintot-1
                exit
            else
                pR(chaintot,j)=int(vect/delta)+1 ! pR tiene la celda donde cae el segmento j de la conf. i
            endif
# else
! **************************************************************************************
! ATENTO: el segundo parametro de int(expression, KIND) fija la cantidad de memoria asignada para el entero
! KIND=1 => 1 byte! = 2^8 = 256 enteros
! viejo       pR(i,j)=int(vect/delta,1)+1 ! pR tiene la celda donde cae el segmento j de la conf. i
! **************************************************************************************
            pR(i,j)=int(vect/delta)+1 ! pR tiene la celda donde cae el segmento j de la conf. i
#       ifdef fdebug_pxs
            print*, "pxs L22: i, j, pR(i,j), vect, in1(i,j,1), in1(i,j,2) ", i, j, pR(i,j), vect, in1(i,j,1), in1(i,j,2)
#       endif
            if(pR(i, j).gt.dimR) then
                Print*,'ERROR grave en pxs!!!!', pR(i, j)
                stop
            endif         
# endif
        enddo
!         print*, j, pR(10, j)
    enddo
    ! Pongo la ultima conformacion en el segmento pegado a la pared
    pR(chaintot,:) = dimR
    
#ifdef fdebug_pxs
    print* , "long, cuantas, chaintot: ", long, cuantas, chaintot
#endif
end subroutine pxs
