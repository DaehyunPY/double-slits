module input
    implicit none
    
    double precision, parameter :: pi = 2.d0*acos(0.d0)
    complex(kind(0.d0)), parameter :: im = cmplx(0.d0, 1.d0, kind(0.d0))

    double precision, parameter :: m = 1.d0
    double precision, parameter :: hb = 1.d0
    double precision, parameter :: e = 1.d0
    double precision, parameter :: phase_shift = 0.d0
    double precision, parameter :: pimentum_shift = 2.d0

end module input





module md_h1
    use input
    use fftsg
	implicit none 

    integer, parameter :: n = 2**8

    double precision, parameter :: x_a = 1.d2 !-a~a 
	integer, parameter :: x_n = 2**10
    double precision, parameter :: x_d = (2.d0*x_a)/dble(x_n)
    double precision, save :: x(0:x_n)

    double precision, parameter :: x_slit = -x_a*50.d-2 !double slit
    double precision, parameter :: x_sigma = x_a*1.d-2 !double slit

    double precision, parameter :: k_a = (2.d0*pi)/(2.d0*x_a)*(dble(x_n)/2.d0) !-a~a
    integer, parameter :: k_n = x_n
    double precision, parameter :: k_d = (2.d0*pi)/(2.d0*x_a)
    double precision, save :: k(0:k_n)

    double precision, save :: psi_steady(0:x_n, 0:n)
    complex(kind(0.d0)), save :: psi_x(0:x_n), psi_k(0:k_n), psi_a(0:n)

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

    f = (2.d0/(2.d0*x_a))**0.5d0*sin((dble(n)*pi)*(x +x_a)/(2.d0*x_a))

end function steady



subroutine hamilton
    double precision :: norm
    integer :: i, j
    
    x(:) = 0.d0
    do i = -x_n/2, x_n/2
        x(x_n/2 +i) = x_d*dble(i)
    enddo

    k(:) = 0.d0
    do i = -k_n/2, k_n/2
        k(k_n/2 +i) = k_d*dble(i)
    enddo

    psi_steady(:, :) = 0.d0
    do i = 0, n
        do j = 0, x_n
            psi_steady(j, i) = steady(x(j), i)*x_d**0.5d0
        enddo
    enddo

    psi_x(:) = 0.d0
    norm = 0.d0
    do i = 0, x_n
        psi_x(i) = gauss(x(i), x_slit, x_sigma)**0.5d0*x_d**0.5d0 &
                    *exp(-im*pimentum_shift/hb*x(i) -im*phase_shift/hb)
        norm = norm &
                +dble(conjg(psi_x(i))*psi_x(i))
    enddo
    write(*, *) int(0), dble(0.d0), 'x', norm
!     psi_x(:) = psi_x(:)/norm**0.5d0

    psi_k(:) = psi_x(:)
    norm = 0.d0
    call fft1d(-1, psi_k(1:))
    do i = 0, k_n
        norm = norm &
                +dble(conjg(psi_k(i))*psi_k(i))
    enddo
    write(*, *) int(0), dble(0.d0), 'k', norm
!     psi_k(:) = psi_k(:)/norm**0.5d0

    psi_a(:) = 0.d0
    norm = 0.d0
    call fftad(-1, psi_x(0:), psi_a(0:))
    do i = 0, n
        norm = norm &
                +dble(conjg(psi_a(i))*psi_a(i))
    enddo
    write(*, *) int(0), dble(0.d0), 'a', norm
!     psi_a(:) = psi_a(:)/norm**0.5d0

!     norm = 0.d0
!     do i = 0, n
!         norm = norm &
!                 +dble(conjg(psi_a(i))*psi_a(i))
!     enddo
!     write(*, *) int(0), dble(0.d0), 'a1', norm

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
end module md_h1










module md_h2
    use input
    implicit none

    double precision, parameter :: t_a = 1.0d2 !jikan hatten no ookisa
    integer, parameter :: t_n = 2**10
    double precision, parameter :: t_d = t_a/dble(t_n)
    double precision, save :: t(0:t_n)

contains



subroutine time
    integer :: i
    
    t(:) = 0.d0
    do i = 0, t_n
        t(i) = t_d*dble(i)
    enddo

end subroutine time
end module md_h2


















































program main 
    use input
    use md_h1
    use md_h2
    use fftsg
    implicit none

!     character(3) :: ch3
    integer, parameter :: file_x = 11
    integer, parameter :: file_k = 12
    integer, parameter :: file_a = 13
    integer, parameter :: file_energy = 14

    integer, parameter :: x_plot = 100
    integer, parameter :: x_pp = x_n/x_plot  
    integer, parameter :: k_plot = 100
    integer, parameter :: k_pp = k_n/k_plot  
    integer, parameter :: a_plot = 100
    integer, parameter :: a_pp = n/a_plot  
    integer, parameter :: t_plot = 100
    integer, parameter :: t_pp = t_n/t_plot  

    complex(kind(0.d0)) :: psi_tmp(0:n)
    double precision :: norm
    integer :: i, j, l 

    call hamilton
    call time



    open(file_x, file = 'output/psi_x.d')
    open(file_k, file = 'output/psi_k.d')
    open(file_a, file = 'output/psi_a.d')
    open(file_energy, file = 'output/energy.d')



    do i = 0, t_n
        if(i /= 0) then
            psi_tmp(:) = psi_a(:)
            psi_a(:) = 0.d0
            norm = 0.d0
            do j = 0, n
                psi_a(j) = exp(-im/hb &
                                    *((dble(j)*pi)/(2.d0*x_a))**2.d0*(hb**2.d0)/(2.d0*m) &
                                    *t_d) &
                            *psi_tmp(j)
                norm = norm &
                        +dble(conjg(psi_a(j))*psi_a(j))
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'a', norm
!             psi_a(:) = psi_a(:)/norm**0.5d0

            psi_x(:) = 0.d0
            norm = 0.d0
            call fftad(+1, psi_x(0:), psi_a(0:))
            do j = 0, x_n
                norm = norm &
                        +dble(conjg(psi_x(j))*psi_x(j))
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'x', norm
!             psi_x(:) = psi_x(:)/norm**0.5d0

            psi_k(:) = psi_x(:)
            norm = 0.d0
            call fft1d(-1, psi_k(1:))
            do j = 0, k_n
                norm = norm &
                        +dble(conjg(psi_k(j))*psi_k(j))
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'k', norm
!             psi_k(:) = psi_k(:)/norm**0.5d0
        endif



        do j = 1, x_n
            if(mod(i, t_pp) == 0 .and. mod(j, x_pp) == 0) then 
                write(file_x, *) t(i), x(j), dble(conjg(psi_x(j))*psi_x(j)/x_d)
            endif
        enddo
        write(file_x, *)

        do j = 1, k_n
            if(mod(i, t_pp) == 0 .and. mod(j, k_pp) == 0) then 
                write(file_k, *) t(i), k(j), dble(conjg(psi_k(j))*psi_k(j)/x_d)
            endif
        enddo
        write(file_k, *)

        do j = 1, n
            if(mod(i, t_pp) == 0 .and. mod(j, a_pp) == 0) then 
                write(file_a, *) t(i), j, dble(conjg(psi_a(j))*psi_a(j))
            endif
        enddo
        write(file_a, *)
    enddo



    close(file_x)
    close(file_k)
    close(file_a)
end program main