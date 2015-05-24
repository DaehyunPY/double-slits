module input
    use fftsg
	implicit none 

    double precision, parameter :: pi = 2.d0*acos(0.d0)
    complex(kind(0.d0)), parameter :: im = cmplx(0.d0, 1.d0, kind(0.d0))

    double precision, parameter :: m = 1.d0
    double precision, parameter :: hb = 1.d0
    double precision, parameter :: e = 1.d0
    double precision, parameter :: phase_shift = 0.d0
    double precision, parameter :: pimentum_shift = 6.d0





    double precision, parameter :: x1_a = 1.d4 !-a~a
	integer, parameter :: x1_n = 2**8 !-a~a
    double precision, parameter :: x1_d = (2.d0*x1_a)/dble(x1_n -1)
    double precision, save :: x1(1:x1_n)

    double precision, parameter :: x2_a = x1_a*0.5d0 !-2a~2a
    integer, parameter :: x2_n = x1_n !0~a
    double precision, parameter :: x2_d = (4.d0*x2_a)/dble(4*x2_n -1)
    double precision, save :: x2(1:4*x2_n)





    integer, parameter :: k1_n = x1_n
    double precision, parameter :: k1_d = (2.d0*pi)/(2.d0*x1_a)
    double precision, parameter :: k1_a = (2.d0*pi)/(2.d0*x1_a)*dble(k1_n)/2.d0 !-pi~pi(-1)
    double precision, save :: k1(1:k1_n)

    integer, parameter :: k2_n = 4*x2_n
    double precision, parameter :: k2_d = (2.d0*pi)/(4.d0*x2_a)
    double precision, parameter :: k2_a = (2.d0*pi)/(4.d0*x2_a)*dble(k2_n)/2.d0 !-pi~pi(-1)
    double precision, save :: k2(1:k2_n)





    double precision, parameter :: x1_slit1 = x1_a*0.01d0 !% double slit
    double precision, parameter :: x1_slit2 = -x1_slit1 !% double slit
    double precision, parameter :: x2_slit = -x2_a !% double slit
    double precision, parameter :: x1_sigma = x1_a*0.001d0 !% double slit
    double precision, parameter :: x2_sigma = x2_a*0.3d0 !% double slit





    complex(kind(0.d0)), save :: psi_x(1:x1_n, 1:4*x2_n)
    complex(kind(0.d0)), save :: psi_k(1:k1_n, 1:k2_n)

	double precision, parameter :: poten_e = -1.d30
	double precision, parameter :: poten_sigma = x1_a*0.03d0 !% double slit
    double precision, save :: poten(1:x1_n, 1:4*x2_n)
    double precision, save :: kinet(1:k1_n, 1:k2_n)





    double precision, parameter :: t_a = 1.d5
    integer, parameter :: t_n = 2**8
    double precision, parameter :: t_d = t_a/dble(t_n)
    double precision, save :: t(0:t_n)

contains










function gauss(x, mu, sigma) result(f)
    double precision, intent(in) :: x, mu, sigma
    double precision :: f

        f = 1.d0/(sigma*(2.d0*pi)**0.5d0) &
                *exp(-(x -mu)**2.d0/(2.d0*sigma**2.d0))        

end function gauss





subroutine hamilton
    complex(kind(0.d0)) :: psi_tmp(1:x1_n, 1:4*x2_n)
    double precision :: norm
    integer :: i, j
    
    x1(:) = 0.d0
    do i = 1, x1_n
        x1(i) = -x1_a +x1_d*dble(i -1)
    enddo

    x2(:) = 0.d0
    do i = 1, 4*x2_n
        x2(i) = -2.d0*x2_a +x2_d*dble(i -1)
    enddo

    k1(:) = 0.d0
    do i = 1, k1_n
        k1(i) = -k1_a +k1_d*dble(i -1)
    enddo

    k2(:) = 0.d0
    do i = 1, k2_n
        k2(i) = -k2_a +k2_d*dble(i -1)
    enddo

    t(:) = 0.d0
    do i = 0, t_n
        t(i) = t_d*dble(i)
    enddo

    poten(:, :) = 0.d0
    do i = 1, 4*x2_n
    	poten(:, i) = poten_e &
                        *gauss(x2(i), 0.d0, poten_sigma)
    enddo

    kinet(:, :) = 0.d0
    do i = 1, k1_n
        do j = 1, k2_n
            kinet(i, j) = (k1(i)**2.d0 +k2(j)**2.d0) &
                            *hb**2.d0/(2.d0*m)
        enddo
    enddo





    psi_x(:, :) = 0.d0
    do i = 1, x1_n
        do j = 1, 4*x2_n
            psi_x(i, j) = (gauss(x1(i), x1_slit1, x1_sigma)**0.5d0*x1_d**0.5d0 &
                                    *gauss(x2(j), x2_slit, x2_sigma)**0.5d0*x2_d**0.5d0 &
                                    +gauss(x1(i), x1_slit2, x1_sigma)**0.5d0*x1_d**0.5d0 &
                                    *gauss(x2(j), x2_slit, x2_sigma)**0.5d0*x2_d**0.5d0) &
                            *exp(-im*pimentum_shift/hb*x2(j) -im*phase_shift/hb)
        enddo
    enddo

    norm = 0.d0
    do i = 1, x1_n
        do j = 1, 4*x2_n
            norm = norm &
                    +dble(conjg(psi_x(i, j))*psi_x(i, j))
        enddo
    enddo
    write(*, *) int(0), dble(0.d0), 'x', norm
    psi_x(:, :) = psi_x(:, :)/norm**0.5d0

    psi_tmp(:, :) = psi_x(:, :)
    psi_x(:, :) = 0.d0
    do i = 1, x1_n
        do j = 1, 4*x2_n
            psi_x(i, j) = psi_tmp(i, j) +psi_tmp(i, 4*x2_n +1 -j) !jiyuudan
!             psi_x(i, j) = psi_tmp(i, j) -psi_tmp(i, 4*x2_n +1 -j) !koteidan
        enddo
    enddo

    norm = 0.d0
    do i = 1, x1_n
        do j = 1, 4*x2_n
            norm = norm &
                    +dble(conjg(psi_x(i, j))*psi_x(i, j))
        enddo
    enddo
    write(*, *) int(0), dble(0.d0), 'x', norm

end subroutine hamilton
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
    integer, parameter :: t_plot = 20
    integer, parameter :: t_pp = t_n/t_plot  

    double precision :: norm
    integer :: i, j, l 
    character(2) :: ch2

    call hamilton





    open(000, file = 'output/poten.d')
    do i = 1, x1_n
    	do j = 1, 4*x2_n
            if(mod(i, x1_pp) == 0 .and. mod(j, x2_pp) == 0) then 
	    		write(0, *) x1(i), x2(j), poten(i, j)
	    	endif
    	enddo
    	write(0, *)
    enddo
    close(000)





    do i = 0, t_n
        if(mod(i, t_pp) == 0) then 
            open(100 +int(i/t_pp), file = 'output/psi_x_'//ch2(int(i/t_pp))//'.d')
!             open(200 +int(i/t_pp), file = 'output/psi_k_'//ch2(int(i/t_pp))//'.d')
        endif

        if(i /= 0) then
            !poten 1
            do j = 1, x1_n
            	do l = 1, 4*x2_n
            		psi_x(j, l) = exp(-im/hb &
	                                    *poten(j, l) &
	                                    *t_d/2.d0) &
	                                *psi_x(j, l)
				enddo
			enddo

			psi_k(:, :) = psi_x(:, :)
	        call fft2d(-1, psi_k(:, :))

            !kinet norm
            norm = 0.d0
            do j = 1, k1_n
                do l = 1, k2_n
                    norm = norm &
                            +dble(conjg(psi_k(j, l))*psi_k(j, l))
                enddo
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'k', norm

        	!kinet 1, 2
            do j = 1, k1_n
                do l = 1, k2_n
                    psi_k(j, l) = exp(-im/hb &
                                        *kinet(j, l) &
                                        *t_d) &
                                    *psi_k(j, l)
                enddo
            enddo

!             !kinet norm
!             norm = 0.d0
!             do j = 1, k1_n
!                 do l = 1, k2_n
!                     norm = norm &
!                             +dble(conjg(psi_k(j, l))*psi_k(j, l))
!                 enddo
!             enddo
!             write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'k', norm

            psi_x(:, :) = psi_k(:, :)
            call fft2d(+1, psi_x(:, :))

            !poten norm
            norm = 0.d0
            do j = 1, x1_n
                do l = 1, 4*x2_n
                    norm = norm &
                            +dble(conjg(psi_x(j, l))*psi_x(j, l))
                enddo
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'x', norm

			!poten 2
            do j = 1, x1_n
                do l = 1, 4*x2_n
                    psi_x(j, l) = exp(-im/hb &
                                        *poten(j, l) &
                                        *t_d/2.d0) &
                                    *psi_x(j, l)
                enddo
            enddo

            !poten norm
            norm = 0.d0
            do j = 1, x1_n
                do l = 1, 4*x2_n
                    norm = norm &
                            +dble(conjg(psi_x(j, l))*psi_x(j, l))
                enddo
            enddo
            write(*, *) int(i/t_pp), dble(i)/dble(t_n)*100.d0, 'x', norm
        endif





        !output
        do j = 1, x1_n
            do l = 1, 2*x2_n
                if(mod(i, t_pp) == 0 .and. mod(j, x1_pp) == 0 .and. mod(l, x2_pp) == 0) then 
                    write(100 +int(i/t_pp), *) & 
                        x1(j), x2(l), &
                        dble(conjg(psi_x(j, l))*psi_x(j, l))/(x1_d*x2_d)**0.5d0
                endif
            enddo
            write(100 +int(i/t_pp), *)
        enddo
        write(100 +int(i/t_pp), *)

!         do j = 1, k1_n
!             do l = 1, k2_n
!                 if(mod(i, t_pp) == 0 .and. mod(j, k1_pp) == 0 .and. mod(l, k2_pp) == 0) then 
!                     write(200 +int(i/t_pp), *) & 
!                         k1(j), k2(l), &
!                         dble(conjg(psi_k(j, l))*psi_k(j, l))/(k1_d*k2_d)**0.5d0
!                 endif
!             enddo
!             write(200 +int(i/t_pp), *)
!         enddo
!         write(200 +int(i/t_pp), *)



        if(mod(i, t_pp) == 0) then 
            close(100 +int(i/t_pp))
!             close(200 +int(i/t_pp))
        endif
    enddo

end program main