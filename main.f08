module input
    implicit none
    
    double precision, parameter ::      pi = 2.d0*acos(0.d0)
    complex(kind(0.d0)), parameter ::   im = cmplx(0.d0, 1.d0, kind(0.d0))

    double precision, parameter ::      m  = 1.d0
    double precision, parameter ::      hb = 1.d0
    double precision, parameter ::      e  = 1.d0

end module input





module md_h1
    use input
	implicit none 

	integer, parameter ::              n   = 2**4 !joutai kazu

    double precision, parameter ::     x_a = 1.d0 !hako no ookisa 0~a
	integer, parameter ::              x_n = 2**6
    double precision, parameter ::     x_d = x_a/dble(x_n)
    double precision, save ::          x(0:x_n)

    double precision, parameter ::     k_a = 2.d0*pi/x_a*(dble(x_n)/2.d0) !hasuu kuukann -a~a
    integer, parameter ::              k_n = x_n
    double precision, parameter ::     k_d = 2.d0*pi/x_a
    double precision, save ::          k(0:k_n)

    complex(kind(0.d0)), save ::       psi_a(1:n) 
    double precision, save ::          psi_steady(1:n, 0:x_n)
    complex(kind(0.d0)), save ::       psi_x(0:x_n), psi_k(0:k_n) !psi_x = psi_a * psi_steady
    double precision, save ::          ope_1(1:n), ope_2(1:n, 1:n)

contains



subroutine hamilton
    integer :: i, j, l
    
    x(:) = 0.d0
    do i = 0, x_n
        x(i) = x_d*dble(i)
    enddo

    k(:) = 0.d0
    do i = -k_n/2, k_n/2
        k(k_n/2 +i) = k_d*dble(i)
    enddo

    psi_a(:) = 0.d0
    psi_steady(:, :) = 0.d0
    do i = 0, x_n
        do j = 1, n
            psi_steady(j, i) = (2.d0/x_a)**0.5d0 &
                                    *sin(dble(j)*pi*x(i)/x_a) &
                                    *x_d**0.5d0
        enddo
    enddo

    ope_1(:) = 0.d0
    ope_2(:, :) = 0.d0
    do i = 1, n
        ope_1(i) = -(dble(n)*pi/x_a)**2.d0
        do j = 1, n
            do l = 0, x_n
                ope_2(j, i) = ope_2(j, i) &
                                +psi_steady(j, l)*x(l)*psi_steady(i, l)
            enddo
        enddo
    enddo

end subroutine hamilton
end module md_h1










module md_h2
    use input
    implicit none

    double precision, parameter ::  t_a     = 2.d0 !jikan hatten no ookisa
    integer, parameter ::           t_n     = 2**18
    double precision, parameter ::  t_d     = t_a/dble(t_n)
    double precision, save ::       t(0:t_n)

!     double precision, parameter ::  denba_a = 10.d-4 !denba parusu no nagasa
    integer, parameter ::           denba_n = 2**8
    double precision, parameter ::  denba_e = 1.d0 !denba no tuyosa
    double precision, parameter ::  denba_w = 1.d0 !denba no sindousuu tani
    double precision, save ::       denba(0:t_n)

contains



function denba_pulse(t) result(f)
    double precision, intent(in) :: t
    double precision :: f

!         if(0 <= t .and. t <= denba_a) then
!             f = sin(pi/denba_a*t)**2.d0
!         else 
!             f = 0.d0
!         endif
    f = 1.d0 !not pulse

end function denba_pulse



function denba_total(t) result(f)
    double precision, intent(in) :: t
    double precision :: f, simple
    integer :: i

    simple = 0.d0
    do i = 1, denba_n
    	simple = simple &
        			+sin(dble(i)*denba_w*t) &
                    /dble(denba_n)
    enddo
    f = denba_e*denba_pulse(t)*simple

end function denba_total



subroutine time
    integer :: i
    
    t(:) = 0.d0
    denba(:) = 0.d0
    do i = 0, t_n
        t(i) = t_d*dble(i)
        denba(i) = denba_total(t(i))
    enddo

end subroutine time
end module md_h2










module md_eq
    use input
    use md_h1
    use md_h2
!     use func
    implicit none 

contains



function eq(x, a) result(f) !equation d/dx a(x) = f(x, a)
    double precision, intent(in) :: x
    complex(kind(0.d0)), intent(in) :: a(1:)
    complex(kind(0.d0)) :: f(1:size(a(1:)))
    integer :: i, j

    f(:) = 0.d0
    do i = 1, size(a(1:))
        f(i) = f(i) +a(i) &
                *im/hb*hb**2.d0/(2.d0*m) &
                *ope_1(i)
        do j = 1, size(a(1:))
            f(i) = f(i) +a(j) &
!                     *im/hb*e*f_11(t(0:), denba(0:), x) & !func hituyou
                    *im/hb*e*denba_total(x) &
                    *ope_2(i, j)
        enddo
    enddo 

end function eq
end module md_eq






























program main 
    use md_h1
    use md_h2
    use md_eq
    use sde
    use fftsg
    implicit none

    integer, parameter :: file_psi_x  = 11
    integer, parameter :: file_psi_k  = 12
    integer, parameter :: file_psi_a  = 13
    integer, parameter :: file_denba  = 21
    integer, parameter :: file_energy = 22

    integer, parameter :: t_plot      = 100
    integer, parameter :: t_pp        = t_n/t_plot  

    double precision, save :: norm, energy
    integer :: i, j, l

    call hamilton
    call time

    open(file_psi_x, file = 'output/psi_x.d')
    open(file_psi_k, file = 'output/psi_k.d')
    open(file_psi_a, file = 'output/psi_a.d')
    open(file_denba, file = 'output/denba.d')
    open(file_energy, file = 'output/energy.d')



    psi_a(1) = 1.d0

    do i = 0, t_n
        if(i /= 0) call rk8_vec(eq, t_d, t(i -1), psi_a(:))

        psi_x(0:x_n) = 0.d0
        do j = 0, x_n
            do l = 1, n
                psi_x(j) = psi_x(j) &
                            +psi_a(l)*psi_steady(l, j)
            enddo
            if(mod(i, t_pp) == 0) write(file_psi_x, *) t(i), x(j), dble(conjg(psi_x(j))*psi_x(j))/x_d
        enddo
        if(mod(i, t_pp) == 0) write(file_psi_x, *)
        
        norm = 0.d0
        do j = 1, n
            norm = norm &
                    +dble(conjg(psi_a(j))*psi_a(j))
            if(mod(i, t_pp) == 0) write(file_psi_a, *) t(i), j, dble(conjg(psi_a(j))*psi_a(j))
        enddo
        if(mod(i, t_pp) == 0) write(file_psi_a, *)

        psi_k(:) = psi_x(:)
        call fft(-1, psi_k(1:))
        do j = 0, k_n
            if(mod(i, t_pp) == 0) write(file_psi_k, *) t(i), k(j), dble(conjg(psi_k(j))*psi_k(j))/k_d
        enddo
!         call fft(+1, psi_k(1:)) !fft test
!         do j = 0, x_n
!             if(mod(i, t_pp) == 0) then
!                 write(file_psi_k, *) t(i), x(j), dble(conjg(psi_k(j))*psi_k(j))/x_d
!             endif
!         enddo
        if(mod(i, t_pp) == 0) write(file_psi_k, *)

        if(mod(i, t_pp) == 0) write(file_denba, *) t(i), denba(i)

        energy = 0.d0
        do j = 1, n
            energy = energy &
                        +dble(conjg(psi_a(j))*psi_a(j)) &
                        	*(pi*hb*dble(j)/x_a)**2.d0/(2.d0*m)
		enddo
        if(mod(i, t_pp) == 0) write(file_energy, *) t(i), energy

		write(*, *) dble(i)/dble(t_n)*100.d0, norm
    enddo




    close(file_psi_x)
    close(file_psi_a)
    close(file_psi_k)
    close(file_denba)
    close(file_energy)

end program main