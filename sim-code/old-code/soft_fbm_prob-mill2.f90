PROGRAM soft_fbm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   time-discrete reflected fractional Brownian motion
!   on finite 1d interval
!
!   fbm_intv1   :  21 Dec 2018      first version based on disk program
!   fbm_intv2   :  23 Dcc 2018      determine distribution over time interval not just at the end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preprocessor directives
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define PARALLEL
#define VERSION 'Soft FBM with Probabilistic Force'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      use IEEE_ARITHMETIC
      use iso_fortran_env
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 
      integer,parameter      :: i8b= SELECTED_INT_KIND(18)            ! 8-byte integers 
      integer,parameter      :: r16b = real128

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(i4b), parameter     :: M=26,NT=2**M               ! number of time steps (in which mu const) 
      integer(i4b), parameter     :: NCONF=25000                    ! number of walkers
      real(r8b), parameter        :: GAMMA = 1.0D0              ! FBM correlation exponent 

      real(r8b)                   :: force_weight = -0.25D0        ! multiplied against density gradient (negative as repelled by density)
      real(r16b), parameter       :: p_accept = 0.8D0
      logical, parameter          :: WRITE_OUTPUT = .FALSE.
      logical, parameter          :: SET_SUM_AVG = .TRUE.             ! if true, we average over the sum of walks from a set, rather 
                                                                      ! rather than average over each individual walker in all set 

      real(r8b),parameter         :: L = 10000000.D0                 ! length of interval
      real(r8b),parameter         :: X0= 0.D0                  ! starting point

      real(r8b), parameter        :: STEPSIG=1.D0             ! sigma of individual step
      character(3)                :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      real(r8b),parameter         :: lambda = 0.2D0/STEPSIG
      real(r8b),parameter         :: wall_force = STEPSIG
      character(4)                :: WALL = 'HARD'

      logical,parameter           :: WRITEDISTRIB = .TRUE.        ! write final radial distribution    
      integer(i4b), parameter     :: NBIN =  5000000                   ! number of bins for density distribution
      integer(i4b), parameter     :: NTSTART=0          ! begin and end of measuring distribution
      integer(i4b), parameter     :: NTEND=2**26 

      real(r8b), parameter        :: outtimefac=2**0.25D0          ! factor for consecutive output times  
            
      integer(i4b), parameter     :: IRINIT=1                    ! random number seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Internal constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      real(r8b), parameter        :: LBY2=L/2
      real(r8b), parameter        :: Pi=3.14159265358979323D0
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer(i4b)	    :: num_passed, sum_num_passed
      real(r8b)              :: xx(0:NT)                        ! walker coordinates
      real(r8b)              :: xix(1:2*NT)                       ! increments

      real(r16b)              :: confxx(1:NT),conf2xx(1:NT),conf4xx(1:NT)                    
      real(r16b)              :: sumxx(1:NT),sum2xx(1:NT),sum4xx(1:NT)         ! sums over machines in MPI version
      !real(r8b)              :: auxxx(1:NT),aux2xx(1:NT)

      !real(r8b)              :: local_corr(1:NT)                  ! product of FBM step and gradient step at each time for a conf
      !real(r8b)              :: global_corr, xix_sum, grad_sum                      ! product of sum of FBM steps and sum of gradient steps for a conf
      !real(r8b)              :: sum_local(1:NT),aux_local(1:NT)
      !real(r8b)              :: sum_global, aux_global

      real(r8b)              :: p_accept_conf
      real(r8b)              :: grad                             ! density gradient 
      real(r16b)             :: conf_weight2, sum_weight2
      real(r8b)              :: force_step
      integer(i4b)           :: iconf, it, ibin_curr, ibin_new, w, test_inc                       ! configuration, and time counters   
      real(r16b)		    :: temp_weight
      integer(i4b)           :: totconf                         ! actual number of confs

      real(r8b)              :: config_history(-NBIN:NBIN)       ! denisty histogram used for gradient calculations 

      real(r8b), allocatable   :: xxdis(:)                ! density histogram 
      real(r8b), allocatable   :: sumdis(:)
      !real(r8b), allocatable   :: auxdis(:) 
      real(r8b)                :: PP,PPsym,x                          ! P(x) 
                        
      external               :: kissinit 
      real(r8b),external     :: rkiss05,erfcc 

      integer(i8b)           :: tlast, tnow, tcount ! used to record time of each configuration



      character(15)          :: avxfile = 'avx00000000.dat'
      character(15)          :: disfile = 'dis00000000.dat' 
            
! Now the MPI stuff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      include 'mpif.h'
      integer(i4b)              :: ierr
      integer(i4b)              :: id,myid                  ! process index
      integer(i4b)              :: numprocs              ! total number of processes
      integer(i4b)              :: status(MPI_STATUS_SIZE)      
#endif             
  

! Start of main program !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set up MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PARALLEL
      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      totconf=(NCONF/numprocs)*numprocs
      
      if (myid==0) then
         print *,'Program ',VERSION,' runing on', numprocs, ' processes'
	 print *, 'mapping factor=', mapping_factor, ', cutoff_accept=', cutoff_accept 
	print *, HUGE(p_accept), log(HUGE(p_accept)), log(1.0D50)
	print *,  log(HUGE(p_accept)) - log(1.0D50)
         print *,'--------------------------------------------------'
      endif ! of if (myid==0)
#else
      totconf=NCONF
      print *,'Program ',VERSION,' runing on single processor'
      print *,'--------------------------------------------------'
#endif 

      if(WRITEDISTRIB) then 
            allocate(xxdis(-NBIN:NBIN),sumdis(-NBIN:NBIN))
            xxdis(:) = 0.D0
      endif 

      confxx(:)=0.D0 
      conf2xx(:)=0.D0
      conf4xx(:)=0.D0

      config_history(:) = 0.D0

      sum2xx(:) = 0.D0
      sumxx(:) = 0.D0
      sum4xx(:) = 0.D0

      sum_weight2 = 0.D0
      conf_weight2 = 0.D0

      !global_corr = 0.D0
      !local_corr(:) = 0.D0
       num_passed = 0

! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
      disorder_loop: do iconf=myid+1,totconf,numprocs
!      if (myid==0) print *, 'dis. conf.', iconf
      if (myid==0) then
            if (iconf==1) then 
              call system_clock(tnow,tcount)
              write(*,'(A,I0,A,I0)') 'dis. conf. ', iconf,' of ', totconf
            else
              tlast=tnow
              call system_clock(tnow)
              write(*,'(A,I0,A,I0,A,I0,A,F0.3,A)') 'dis. conf. ', iconf,' of ',totconf,&
                ' (took ',(tnow-tlast)/(60*tcount),' minutes and ',mod(tnow-tlast,60*tcount)/(tcount*1.D0),' seconds)'
            endif
      endif

#else
      disorder_loop: do iconf=1,totconf
!         print *, 'dis. conf.', iconf
#endif 

        call gkissinit(IRINIT+iconf-1)
        call corvec(xix,2*NT,M+1,GAMMA)                         ! create x-increments (correlated Gaussian random numbers)
        
        !if (STEPDIS.eq.'BOX') then
        !  xix(:) = 1 - (0.5D0*erfcc(xix(:)/sqrt(2.0D0)))                       ! map onto unit inteval
        ! xix(:)=xix(:)-0.5D0                                                  ! center around 0
        !endif
        
        if (STEPDIS.eq.'BIN') then 
          xix(:)=sign(1.D0,xix(:))
        endif  
        
        xix(:)=xix(:)*STEPSIG                               ! scale increments to correct sigma 
        
! Time loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            confxx(:)=0.D0 
            conf2xx(:)=0.D0
	      conf4xx(:)=0.D0
            xx(0)=X0
            config_history(:) = 0.D0
            config_history(0) = 1.0D0

            time_loop: do it=1, NT
 
                  ibin_curr =nint( xx(it-1)*NBIN/LBY2 ) ! find walker's starting bin 
                  grad = 0.D0 ! reset grad in case we went past the wall in soft wall case

                  xx(it) = xx(it-1) + xix(it) ! find new position of walker 

                  if (WALL .eq. 'SOFT') then
                        xx(it) = xx(it) + wall_force*exp(-lambda*(xx(it-1)+LBY2)) - wall_force*exp(lambda*(xx(it-1)-LBY2)) 
                  else ! WALL .eq. 'HARD'
                        if ( abs(xx(it)).gt.LBY2 ) then ! stopping boundaries
                              xx(it)=xx(it-1)
                        endif 
                  end if

                  ibin_new = nint( xx(it)*NBIN/LBY2 )

                  ! if we're trying to walk into a higher density bin, go to p_accept 
                  if(config_history(ibin_new) .gt. config_history(ibin_curr)) then 
                        if (rkiss05() > p_accept) then ! if number outside of p_accept is selected
                              xx(it) = xx(it) - xix(it) ! step isn't accepted, don't move anywhere 
                              ibin_new = nint( xx(it)*NBIN/LBY2 ) ! and recalculate new bin
                        end if 
                  end if 
                       
                  confxx(it)=confxx(it) + xx(it)
                  conf2xx(it)=conf2xx(it) + xx(it)*xx(it)
		      conf4xx(it)=conf4xx(it) + xx(it)*xx(it)*xx(it)*xx(it)
                  !xix_sum = xix_sum + xix(it)
                  
                  if (WRITE_OUTPUT) then
                        if (myid==0) then
                              if (iconf==1) then 
                              write(*,'(I0.1, A, F0.7,A,F0.10,A,F0.3,A,F0.3,A,I0.1,A,I0.1,A,I0.1,A)')&
                              it, ' : ',&
                              xx(it), ' = ',&
                              xx(it-1), ' + ',xix(it), ' + ', force_step,&
                              " [", int(config_history(ibin_new-1)),&
                              ",", int(config_history(ibin_new)), &
                              ",", int(config_history(ibin_new+1)), "]"
                              endif
                        endif
                  end if 

                  ! update full time distribution for gradient 
                  config_history(ibin_new) = config_history(ibin_new) + 1.0D0 

                  if ( (ibin_new.ge.-NBIN) .and. (ibin_new.le.NBIN)) then
                        ! record steady-state distribution 
                        if( (it.ge.NTSTART) .and. (it.le.NTEND) .and. WRITEDISTRIB) then
                              xxdis(ibin_new) = xxdis(ibin_new) + 1.0D0
                        end if 
                  end if
            
           end do time_loop

           sumxx(:) = confxx(:) + sumxx(:)
           sum2xx(:) = conf2xx(:) + sum2xx(:)
	     sum4xx(:) = conf4xx(:) + sum4xx(:)   


      end do disorder_loop      ! of do inconf=1,NCONF


! Now collect and analyze data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL
     if (myid.ne.0) then                                                  ! Send data
         call MPI_SEND(sumxx,NT,MPI_REAL16,0,1,MPI_COMM_WORLD,ierr)
         call MPI_SEND(sum2xx,NT, MPI_REAL16,0,2,MPI_COMM_WORLD,ierr)
	

         if(WRITEDISTRIB) then
            call MPI_SEND(xxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
         end if

	call MPI_SEND(sum4xx,NT, MPI_REAL16,0,4,MPI_COMM_WORLD,ierr)

         !call MPI_SEND(global_corr,1,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)
         !call MPI_SEND(local_corr,NT,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
	!print *, "send all data from myid=", myid

      else

         sumdis(:)=xxdis(:)
         xxdis(:)= 0.D0

         conf2xx(:) = 0.D0
         confxx(:) = 0.D0
	   conf4xx(:) = 0.D0

         !sum_global = global_corr
         !sum_local(:) = local_corr(:)

         do id=1,numprocs-1                                                   ! Receive data
            call MPI_RECV(confxx,NT, MPI_REAL16,id,1,MPI_COMM_WORLD,status,ierr)
	    !print *, "received confxx from myid=", id 

            call MPI_RECV(conf2xx,NT, MPI_REAL16,id,2,MPI_COMM_WORLD,status,ierr)
	    !print *, "received conf2xx from myid=", id 

            sumxx(:)=sumxx(:)+confxx(:) 
            sum2xx(:)=sum2xx(:)+conf2xx(:) 

            if(WRITEDISTRIB) then
                  call MPI_RECV(xxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,id,3,MPI_COMM_WORLD,status,ierr)
                  sumdis(:)=sumdis(:)+xxdis(:)
		 !print *, "received xxdis from myid=", id 
            end if 

	   call MPI_RECV(conf4xx, NT, MPI_REAL16, id, 4, MPI_COMM_WORLD, status, ierr)
	   sum4xx(:) = sum4xx(:) + conf4xx(:)

	    print *, "received all data from myid=", id

            !call MPI_RECV(aux_global,1,MPI_DOUBLE_PRECISION,id,4,MPI_COMM_WORLD,status,ierr)
            !call MPI_RECV(aux_local,NT,MPI_DOUBLE_PRECISION,id,5,MPI_COMM_WORLD,status,ierr)
            !sum_global = sum_global + aux_global
            !sum_local(:) = sum_local(:) + aux_local(:)
	    
         enddo
      endif        
#else          
      sumxx(:)=confxx(:)
      sum2xx(:)=conf2xx(:)
      sumdis(:)=xxdis(:)
#endif
         
#ifdef PARALLEL
      if (myid==0) then
#endif       
        write(avxfile(4:11),'(I8.8)') NT       
        open(2,file=avxfile,status='replace')

        write(2,*) 'Program ', VERSION
        write(2,*) 'average displacement of reflected FBM'
        write(2,*) 'weight on force', force_weight
        write(2,*) 'step distribution ',STEPDIS
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'L=', L
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=', totconf
        write(2,*) 'p_accept=', p_accept
      write(2,*) 'SET_SUM_AVG', SET_SUM_AVG
       ! write(2,*) '<sum(grad)*sum(xix)>', sum_global/totconf
        write (2,*)'=================================='
        write(2,*) '   time         <r>         <r^2>'      
        it=1
        do while(it.le.NT)  
          Write(2,'(1X,I8,6(2X,E13.6))')  it, sumxx(it)/sum_weight2, sum2xx(it)/sum_weight2, sum4xx(it)/sum_weight2
          it=max(it+1,nint(outtimefac*it))
        enddo 
        close(2) 

      if (WRITEDISTRIB) then
        write(disfile(4:11),'(I8.8)') NT 

        open(2,file=disfile,status='replace')

        write(2,*) 'Program ', VERSION
        write(2,*) 'Density distribution over entire NT'
        write(2,*) 'STEPSIG=', STEPSIG
        write(2,*) 'GAMMMA=', GAMMA 
        write(2,*) 'NT= ', NT
        write(2,*) 'NTSTART= ', NTSTART
        write(2,*) 'NTEND= ', NTEND
        write(2,*) 'L=', L
        write(2,*) 'IRINIT=',IRINIT
        write(2,*) 'NCONF=',totconf
        write (2,*)'=================================='
        write(2,*) '   ibin    x=(L/2*ibin)/NBIN  x/L   P(x)  P(x)*L  P(|x|)   L/2-x '
        do ibin_curr=-NBIN,NBIN 
          x= (L/2*ibin_curr)/NBIN
          PP= (sumdis(ibin_curr)*NBIN)/(L/2*totconf*(NTEND-NTSTART+1))
          PPsym= ( (0.5D0*sumdis(ibin_curr)+0.5D0*sumdis(-ibin_curr))*NBIN)/(L/2*totconf*(NTEND-NTSTART+1))
          Write(2,'(1X,I9,8(2X,E13.6))') ibin_curr, x, x/L, PP, PP*L, PPsym, L/2-x 
        enddo 
        close(2) 
      endif

#ifdef PARALLEL
      endif ! of if (myid==0)
#endif        
      
      
#ifdef PARALLEL
      call MPI_FINALIZE(ierr)
#endif 
      stop      
      END PROGRAM soft_fbm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine CORVEC(xr,Ns,M)
!
! generates a 1d array of Ns Gaussian random numbers xr
! correlated according to a (translationally invariant)
! user-supplied correlation function corfunc(is,Ns)
! 
! uses Fourier filtering method
!
! history
!      v0.9         Dec  7, 2013:        first version, uses Tao Pang FFT
!      v0.91        Oct 11, 2017:        uses much faster FFT by Ooura (in fftsg.f)        
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE corvec(xr,Ns,M,gam)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 

      integer(i4b)           :: Ns              ! number of sites, must be power of 2 
      integer(i4b)           :: M               ! Ns=2^M

      real(r8b)              :: xr(0:Ns-1)      ! random number array  
      real(r8b)              :: cr(0:Ns-1)      ! correlation function 
      integer(i4b)           :: is
      
      integer(i4b)           :: ip(0:int(2+sqrt(1.*Ns)))   ! workspace for FFT code (added int around second indx)
      real(r8b)              :: w(0:Ns/2-1)           ! workspace for FFT code 

      real(r8b)              :: gam                   ! FBM exponent, pass through to correlation function   
            
      real(r8b), external    :: gkiss05,erfcc
      external               :: rdft                   ! from Ooura's FFT package
      real(r8b),external     :: fbmcorfunc 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (Ns.ne.2**M) STOP 'Size indices do not match'
! Calculate correlation function 
      do is=0,Ns-1 
         cr(is)= fbmcorfunc(is,Ns,gam) 
      enddo
! Real FFT of correlation function
       ip(0)=0
       call rdft(Ns, 1, cr, ip, w)  
       
! Create array of independent Gaussian random numbers
      do is=0,Ns-1
         xr(is)= gkiss05()
      enddo
! Real FFT of input random numbers      
       call rdft(Ns, 1, xr, ip, w)  
! Filter the Fourier components of the random numbers
! as real-space correlations are symmmetric, FT of c is real
      do is=1,Ns/2-1
          xr(2*is)=xr(2*is)*sqrt(abs(cr(2*is)))*2.D0/Ns
          xr(2*is+1)=xr(2*is+1)*sqrt(abs(cr(2*is)))*2.D0/Ns
      enddo
      xr(0)=xr(0)*sqrt(abs(cr(0)))*2.D0/Ns
      xr(1)=xr(1)*sqrt(abs(cr(1)))*2.D0/Ns 
      
! FFT of filtrered random numbers (back to real space)
       call rdft(Ns, -1, xr, ip, w)  
       
! Transform from Gaussian distribution to flat distribution on (0,1)      
!      do is = 0,Ns-1
!        xr(is) = 1 -  0.5D0*erfcc(xr(is)/sqrt(2.0D0)) 
!      end do
      
      return
      END SUBROUTINE corvec 
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FBM correlation function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION fbmcorfunc(i,N,gam)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)              :: corfunc,fbmcorfunc
      real(r8b)              :: gam
      integer(i4b)           :: i,N
      integer(i4b)           :: dist

!      print *, 'gamma=',gam
      dist=min(i,N-i)
      if (dist.eq.0) then 
         corfunc=1.D0
      elseif (dist.lt.1000) then   
         corfunc = 0.5D0*( dble(dist+1)**(2.D0-gam) - 2.D0*(dble(dist)**(2.D0-gam)) + dble(dist-1)**(2.D0-gam) ) 
      else 
         corfunc = (2.D0-gam)*(1.D0-gam)*dble(dist)**(-gam) 
         corfunc = corfunc + (2.D0-gam)*(1.D0-gam)*(-gam)*(-1.D0-gam)*dble(dist)**(-gam-2.D0)/12.D0          
         corfunc = 0.5D0*corfunc
      endif   

      fbmcorfunc=corfunc
      
      return
      END FUNCTION fbmcorfunc 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gaussian random number generator gkiss05
!
! generates normally distributed independent random numbers
! (zero mean, variance 1) using Box-Muller method in polar form
!
! uniform random numbers provided by Marsaglia's kiss (2005 version)
!
! before using the RNG, call gkissinit(seed) to initialize
! the generator. Seed should be a positive integer.
!
!
! History:
!      v0.9     Dec  6, 2013:   first version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION gkiss05()
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      real(r8b)             :: gkiss05
      real(r8b), external   :: rkiss05

      real(r8b)             :: v1,v2,s,fac
      integer(i4b)          :: iset               ! switches between members of the Box-Muller pair
      real(r8b)             :: gset
      common /gausscom/gset,iset

      if (iset.ne.1) then
        do
          v1 = 2.D0 * rkiss05() - 1.D0
          v2 = 2.D0 * rkiss05() - 1.D0
          s = v1 * v1 + v2 * v2
          if ((s<1.D0) .and. (s>0.D0)) exit
        enddo
! Box-Muller transformation creates pairs of random numbers
        fac = sqrt(-2.D0 * log(s) / s)
        gset = v1 * fac
        iset = 1
        gkiss05 = v2 * fac
      else
        iset = 0
        gkiss05 = gset
      end if
      return
      END FUNCTION gkiss05


      SUBROUTINE gkissinit(iinit)
      implicit none
      integer,parameter     :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b)          :: iinit,iset
      real(r8b)             :: gset
      common /gausscom/gset,iset

      iset=0                         ! resets the switch between the members of the Box-Muller pair
      call kissinit(iinit)           ! initializes the rkiss05 RNG
      end subroutine gkissinit


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Random number generator KISS05 after a suggestion by George Marsaglia
! in "Random numbers for C: The END?" posted on sci.crypt.random-numbers
! in 1999
!
! version as in "double precision RNGs" in  sci.math.num-analysis
! http://sci.tech-archive.net/Archive/sci.math.num-analysis/2005-11/msg00352.html
!
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
! Overall period > 2^123
!
!
! A call to rkiss05() gives one random real in the interval [0,1),
! i.e., 0 <= rkiss05 < 1
!
! Before using rkiss05 call kissinit(seed) to initialize
! the generator by random integers produced by Park/Millers
! minimal standard LCG.
! Seed should be any positive integer.
!
! FORTRAN implementation by Thomas Vojta, vojta@mst.edu
! built on a module found at www.fortran.com
!
!
! History:
!        v0.9     Dec 11, 2010    first implementation
!        V0.91    Dec 11, 2010    inlined internal function for the SR component
!        v0.92    Dec 13, 2010    extra shuffle of seed in kissinit
!        v0.93    Aug 13, 2012    changed integer representation test to avoid data statements
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      FUNCTION rkiss05()
      implicit none

      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers
      real(r8b),parameter    :: am=4.656612873077392578d-10       ! multiplier 1/2^31

      real(r8b)             :: rkiss05
      integer(i4b)          :: kiss
      integer(i4b)          :: x,y,z,w              ! working variables for the four generators
      common /kisscom/x,y,z,w

      x = 69069 * x + 1327217885
      y= ieor (y, ishft (y, 13)); y= ieor (y, ishft (y, -17)); y= ieor (y, ishft (y, 5))
      z = 18000 * iand (z, 65535) + ishft (z, - 16)
      w = 30903 * iand (w, 65535) + ishft (w, - 16)
      kiss = ishft(x + y + ishft (z, 16) + w , -1)
      rkiss05=kiss*am
      END FUNCTION rkiss05


      SUBROUTINE kissinit(iinit)
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter     :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers

      integer(i4b) idum,ia,im,iq,ir,iinit
      integer(i4b) k,x,y,z,w,c1,c2,c3,c4
      real(r8b)    rkiss05,rdum
      parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
      common /kisscom/x,y,z,w

      !!! Test integer representation !!!
      c1=-8
      c1=ishftc(c1,-3)
!     print *,c1
      if (c1.ne.536870911) then
         print *,'Nonstandard integer representation. Stoped.'
         stop
      endif

      idum=iinit
      idum= abs(1099087573 * idum)               ! 32-bit LCG to shuffle seeds
      if (idum.eq.0) idum=1
      if (idum.ge.IM) idum=IM-1

      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         x=idum+1
      else
         x=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         y=idum+1
      else
         y=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         z=idum+1
      else
         z=idum
      endif
      k=(idum)/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum = idum + IM
      if (idum.lt.1) then
         w=idum+1
      else
         w=idum
      endif

      rdum=rkiss05()

      return
      end subroutine kissinit



 FUNCTION erfcc(x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates complementary error function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      real(r8b)   :: erfcc,x
      real(r8b)   :: t,z
      z=abs(x)
      t=1.D0/(1.D0+0.5D0*z)
      erfcc=t*exp(-z*z-1.26551223D0+t*(1.00002368D0+t*(.37409196D0+t*&
     &(.09678418D0+t*(-.18628806D0+t*(.27886807D0+t*(-1.13520398D0+t*&
     &(1.48851587D0+t*(-.82215223D0+t*.17087277D0)))))))))
      if (x.lt.0.D0) erfcc=2.D0-erfcc
      return
      END  