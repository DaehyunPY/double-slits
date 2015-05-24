module input
    use fftsg
	implicit none 

    double precision, parameter :: pi = 2.d0*acos(0.d0)
    complex(kind(0.d0)), parameter :: im = cmplx(0.d0, 1.d0, kind(0.d0))

    double precision, parameter :: m = 1.d0
    double precision, parameter :: hb = 1.d0
    double precision, parameter :: e = 1.d0
    double precision, parameter :: phase_shift = 0.d0
    double precision, parameter :: pimentum_shift = 5.d0



    double precision, parameter :: x1_a = 1.d3 !-a~a 
	integer, parameter :: x1_n = 2**8
    double precision, parameter :: x1_d = (2.d0*x1_a)/dble(x1_n)
    double precision, save :: x1(0:x1_n)

    double precision, parameter :: x2_a = x1_a*50.d-2 !% 0~a
    integer, parameter :: x2_n = x1_n
    double precision, parameter :: x2_d = x2_a/dble(x2_n)
    double precision, save :: x2(0:x2_n)



    double precision, parameter :: x1_slit = x1_a*1.d-2 !% double slit
    double precision, parameter :: x1_sigma = x1_a*.3d-2 !% double slit
    double precision, parameter :: x2_sigma = x1_a*10.d-2 !% double slit

    double precision, parameter :: x2_ob = x2_a*100.d-2 !% ob
    integer, parameter :: x2_ob_n = int(x2_ob/x2_a*dble(x2_n))



    double precision, parameter :: k1_a = (2.d0*pi)/(2.d0*x1_a)*(dble(x1_n)/2.d0) !-a~a
    integer, parameter :: k1_n = x1_n
    double precision, parameter :: k1_d = (2.d0*pi)/(2.d0*x1_a)
    double precision, save :: k1(0:k1_n)

    double precision, parameter :: k2_a = (2.d0*pi)/(2.d0*x2_a)*(dble(x2_n)/2.d0) !-a~a
    integer, parameter :: k2_n = x2_n
    double precision, parameter :: k2_d = (2.d0*pi)/(2.d0*x2_a)
    double precision, save :: k2(0:k2_n)



    double precision, parameter :: t_a = 3.0d3 !jikan hatten no ookisa
    integer, parameter :: t_n = 2**8
    double precision, parameter :: t_d = t_a/dble(t_n)
    double precision, save :: t(0:t_n)


 
    complex(kind(0.d0)), save :: psi_x(0:x1_n, 0:x2_n)
    complex(kind(0.d0)), save :: psi_k(0:k1_n, 0:k2_n)
    double precision, save :: psi_steady(0:x2_n, 0:k2_n)

contains





function gauss(x, mu, sigma) result(f)
    double precision, intent(in) :: x, mu, sigma
    double precision :: f

        f = 1.d0/(sigma*(2.d0*pi)**0.5d0) &
                *exp(-(x -mu)**2.d0/(2.d0*sigma**2.d0))        

end function gauss





function steady(x, n) result(f)
    double precision, intent(in) :: x
    integer, intent(in) :: n
    double precision :: f

!     f = (1.d0/(2.d0*x2_a))**0.5d0 &
!             *sin((dble(n)*pi)*(x +x2_a)/(2.d0*x2_a)) !jiyuudan-koteidan

    f = (2.d0/x2_a)**0.5d0 &
            *sin(dble(n)*pi*x/x2_a) !koteidan-koteidan

!     f = (2.d0/(x2_a*2.d0))**0.5d0 &
!             *sin(dble(n)*pi*x/(x2_a*2.d0)) !koteidan-jiyuudan

end function steady





subroutine hamilton
    double precision :: norm
    integer :: i, j
    
    x1(:) = 0.d0
    do i = -x1_n/2, x1_n/2
        x1(x1_n/2 +i) = x1_d*dble(i)
    enddo

    x2(:) = 0.d0
    do i = 0, x2_n
        x2(i) = x2_d*dble(i)
    enddo

    k1(:) = 0.d0
    do i = -k1_n/2, k1_n/2
        k1(k1_n/2 +i) = k1_d*dble(i)
    enddo

    k2(:) = 0.d0
    do i = -k2_n/2, k2_n/2
        k2(k2_n/2 +i) = k2_d*dble(i)
    enddo

    psi_steady(:, :) = 0.d0
    do i = 0, k2_n
        do j = 0, x2_n
            psi_steady(j, i) = steady(x2(j), i)*x2_d**0.5d0
        enddo
    enddo

    t(:) = 0.d0
    do i = 0, t_n
        t(i) = t_d*dble(i)
    enddo





    psi_x(:, :) = 0.d0
    do i = 0, x1_n
        do j = 0, x2_n
            psi_x(i, j) = 0.5d0*(gauss(x1(i), x1_slit, x1_sigma)**0.5d0*x1_d**0.5d0 &
                                    *gauss(x2(j), 0.d0, x2_sigma)**0.5d0*x2_d**0.5d0 &
                                    +gauss(x1(i), -x1_slit, x1_sigma)**0.5d0*x1_d**0.5d0 &
                                    *gauss(x2(j), 0.d0, x2_sigma)**0.5d0*x2_d**0.5d0) &
                            *exp(-im*pimentum_shift/hb*x2(j) -im*phase_shift/hb)
        enddo
    enddo

    norm = 0.d0
    do i = 0, x1_n
        do j = 0, x2_n
            norm = norm &
                    +dble(conjg(psi_x(i, j))*psi_x(i, j))
        enddo
    enddo
    write(*, *) int(0), dble(0.d0), 'x', norm
    psi_x(:, :) = psi_x(:, :)/norm**0.5d0

    norm = 0.d0
    do i = 0, x1_n
        do j = 0, x2_n
            norm = norm &
                    +dble(conjg(psi_x(i, j))*psi_x(i, j))
        enddo
    enddo
    write(*, *) int(0), dble(0.d0), 'x', norm



    psi_k(:, :) = 0.d0
    do i = 0, x1_n
        call fftad(-1, psi_x(i, 0:), psi_k(i, 0:))
    enddo
    do i = 0, k2_n
        call fft1d(-1, psi_k(1:, i))
    enddo

    norm = 0.d0
    do i = 0, k1_n
        do j = 0, k2_n
            norm = norm &
                    +dble(conjg(psi_k(i, j))*psi_k(i, j))
        enddo
    enddo
    write(*, *) int(0), dble(0.d0), 'k', norm
!     psi_k(:, :) = psi_k(:, :)/norm**0.5d0

end subroutine hamilton





subroutine fftad(type, psi_x_input, psi_a_input)
    integer, intent(in) :: type
    complex(kind(0.d0)), intent(inout) :: psi_x_input(0:)
    complex(kind(0.d0)), intent(inout) :: psi_a_input(0:)
    integer :: i, j, n_x_input, n_a_input

    n_x_input = size(psi_x_input(0:))
    n_a_input = size(psi_a_input(0:))



    if(type == -1) then
        psi_a_input(:) = 0.d0
        do i = 0, n_a_input -1
            do j = 0, n_x_input -1
                psi_a_input(i) = psi_a_input(i) &
                                    +conjg(psi_steady(j, i)*psi_x_input(j))
            enddo
        enddo
    else if(type == +1) then
        psi_x_input(:) = 0.d0
        do i = 0, n_x_input -1
            do j = 0, n_a_input -1
                psi_x_input(i) = psi_x_input(i) &
                                    +psi_a_input(j)*psi_steady(i, j)
            enddo
        enddo
    endif

end subroutine fftad
end module input


















































program main 
    use input
    use fftsg
    implicit none

    integer, parameter :: x1_plot = 120
    integer, parameter :: x1_pp = x1_n/x1_plot  
    integer, parameter :: x2_plot = 50
    integer, parameter :: x2_pp = x2_n/x2_plot  
    integer, parameter :: k1_plot = 50
    integer, parameter :: k1_pp = k1_n/k1_plot  
    integer, parameter :: k2_plot = 50
    integer, parameter :: k2_pp = k2_n/k2_plot  
    integer, parameter :: t_plot = 10
    integer, parameter :: t_pp = t_n/t_plot  

    double precision :: norm
    integer :: i, j, l 
    character(3) :: ch3

    call hamilton





    do i = 0, t_n
        if(mod(i, t_pp) == 0) then 
            open(1000 +int(i/t_pp), file = 'output/psi_x_'//ch3(int(i/t_pp))//'.d')
            open(2000 +int(i/t_pp), file = 'output/psi_k_'//ch3(int(i/t_pp))//'.d')
        endif



        if(i /= 0) then
            do j = 0, k2_n
                psi_k(:, j) = exp(-im/hb &
                                    *(dble(j)*pi/(2.d0*x2_a))**2.d0 &
                                    *hb**2.d0/(2.d0*m) &
                                    *t_d/2.d0) &
                                *psi_k(:, j)
            enddo

            do j = 0, k1_n
                psi_k(j, :) = exp(-im/hb &
                                    *k1(j)**2.d0 &
                                    *hb**2.d0/(2.d0*m) &
                                    *t_d) &
                                *psi_k(j, :)
            enddo

            do j = 0, k2_n
                psi_k(:, j) = exp(-im/hb &
                                    *(dble(j)*pi/(2.d0*x2_a))**2.d0 &
                                    *hb**2.d0/(2.d0*m) &
                                    *t_d/2.d0) &
                                *psi_k(:, j)
            enddo

            norm = 0.d0
            do j = 0, k1_n
                do l = 0, k2_n
                    norm = norm &
                            +dble(conjg(psi_k(j, l))*psi_k(j, l))
                enddo
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'k', norm
!             psi_k(:, :) = psi_k(:, :)/norm**0.5d0



            psi_x(:, :) = 0.d0
            do j = 0, x1_n
                call fftad(+1, psi_x(j, 0:), psi_k(j, 0:))
            enddo
            do j = 0, k2_n
                call fft1d(+1, psi_x(1:, j))
            enddo

            norm = 0.d0
            do j = 0, k1_n
                do l = 0, k2_n
                    norm = norm &
                            +dble(conjg(psi_x(j, l))*psi_x(j, l))
                enddo
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'x', norm
!             psi_x(:, :) = psi_x(:, :)/norm**0.5d0
        endif



        do j = 1, x1_n
            do l = 1, x2_n
                if(mod(i, t_pp) == 0 .and. mod(j, x1_pp) == 0 .and. mod(l, x2_pp) == 0) then 
                    write(1000 +int(i/t_pp), *) & 
                        x1(j), x2(l), &
                        dble(conjg(psi_x(j, l))*psi_x(j, l))/(x1_d*x2_d)**0.5d0
                endif
            enddo
            write(1000 +int(i/t_pp), *)
        enddo
        write(1000 +int(i/t_pp), *)

        do j = 1, k1_n
            do l = 1, k2_n
                if(mod(i, t_pp) == 0 .and. mod(j, k1_pp) == 0 .and. mod(l, k2_pp) == 0) then 
                    write(2000 +int(i/t_pp), *) & 
                        k1(j), k2(l), &
                        dble(conjg(psi_k(j, l))*psi_k(j, l))/(k1_d*k2_d)**0.5d0
                endif
            enddo
            write(2000 +int(i/t_pp), *)
        enddo
        write(2000 +int(i/t_pp), *)



        if(mod(i, t_pp) == 0) then 
            close(1000 +int(i/t_pp))
            close(2000 +int(i/t_pp))
        endif
    enddo

end program main