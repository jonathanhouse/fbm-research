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
#define VERSION 'Soft FBM Parallel (procs as walkers)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      implicit none
      integer,parameter      :: r8b= SELECTED_REAL_KIND(P=14,R=99)   ! 8-byte reals
      integer,parameter      :: i4b= SELECTED_INT_KIND(8)            ! 4-byte integers 
      integer,parameter      :: i8b= SELECTED_INT_KIND(18)            ! 8-byte integers 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      integer(i4b), parameter     :: M=26,NT=2**M               ! number of time steps (in which mu const) 
     
      real(r8b), parameter        :: GAMMA = 1.0D0              ! FBM correlation exponent 
      integer(i4b), parameter     :: NSETS = 10
      integer(i4b), parameter     :: NWALKS_IN_SET = 200 
      integer(i4b), parameter     :: NCONF=NSETS*NWALKS_IN_SET                    ! number of walkers



      real(r8b)                   :: force_weight = -0.25D0        ! multiplied against density gradient (negative as repelled by density)
      real(r8b), parameter        :: nonlin_factor = 1.D0         ! a*tanh(x/a) where a is nonlinear scale 
      character(6), parameter     :: FORCE_TYPE = 'LINEAR'         ! Nonlinear tanh = 'NONLIN', Linear force = 'LINEAR'
      character(4), parameter     :: GRAD_FORM = 'ASYM'           ! asymmetric gradient -> ASYM; symmetric gradient -> SYMM
      integer(i4b), parameter     :: GRAD_DX = 1                  ! step used in gradient formula : ex. GRAD_DX=1 w/ SYMM is two-point symmetric formula 
      integer(i4b), parameter     :: WINDOW = 3		             ! WINDOW*2 + 1 is width of window
      character(4), parameter     :: GRAD_TEST = 'NONE'
      character(4), parameter     :: FORCE_TEST = 'NONE'           ! random weight drawn from uniform dist -> RAND
      logical, parameter          :: WRITE_OUTPUT = .FALSE.


      real(r8b),parameter         :: L = 10000000.D0                 ! length of interval
      real(r8b),parameter         :: X0= 0.D0                  ! starting point

      real(r8b), parameter        :: STEPSIG=1.D0             ! sigma of individual step
      character(3)                :: STEPDIS='GAU'                 ! Gaussian = GAU, binary = BIN, box = BOX                   

      real(r8b),parameter         :: lambda = 0.2D0/STEPSIG
      real(r8b),parameter         :: wall_force = STEPSIG
      character(4)                :: WALL = 'HARD'

      logical,parameter           :: WRITEDISTRIB = .TRUE.        ! write final radial distribution    
      integer(i4b), parameter     :: NBIN =  5000000                   ! number of bins for density distribution
      integer(i4b), parameter     :: NTSTART=2**26-1000          ! begin and end of measuring distribution
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

      !real(r8b)              :: xx(0:NT)                        ! walker coordinates
      real(r8b)              :: xix(1:2*NT)                       ! increments

      real(r8b)              :: confxx(1:NT)                     ! average of xx after each time step  
      real(r8b)              :: conf2xx(1:NT)                    
      real(r8b)              :: sumxx(1:NT),sum2xx(1:NT)         ! sums over machines in MPI version
      !real(r8b)              :: auxxx(1:NT),aux2xx(1:NT) 

      real(r8b)		   :: window_mu, bin_mu, bin_mu2
      real(r8b)		   :: window_covar

      !real(r8b)              :: local_corr(1:NT)                  ! product of FBM step and gradient step at each time for a conf
      !real(r8b)              :: global_corr, xix_sum, grad_sum                      ! product of sum of FBM steps and sum of gradient steps for a conf
      !real(r8b)              :: sum_local(1:NT),aux_local(1:NT)
      !real(r8b)              :: sum_global, aux_global

      real(r8b)              :: grad                             ! density gradient 
      real(r8b)              :: force_step
      integer(i4b)           :: iconf, it, ibin, w, iset, iwalker, i                       ! configuration, and time counters   
      integer(i4b)           :: totconf,totwalks_in_set                         ! actual number of confs

      real(r8b)               :: set_history(-NBIN:NBIN)       ! denisty histogram used for gradient calculations 
      !real(r8b)               :: temp_history(-NBIN:NBIN)
      real(r8b), allocatable  :: temp_xx(:)
      real(r8b), allocatable  :: walkers_xix(:,:)

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

      !! actual number of walkers/set is an integer multiple of numprocs, likewise for totconf
      totwalks_in_set=(NWALKS_IN_SET/numprocs)*numprocs 
      totconf = totwalks_in_set*NSETS
      
      if (myid==0) then
         print *,'Program ',VERSION,' running on', numprocs, ' processes'
         print *,'--------------------------------------------------'
      endif ! of if (myid==0)
#else
      !totconf=NCONF
      print *,'Program ',VERSION,' running on single processor'
      print *,'--------------------------------------------------'
#endif 

      if(WRITEDISTRIB) then 
            allocate(xxdis(-NBIN:NBIN),sumdis(-NBIN:NBIN))
            !allocate(auxdis(-NBIN:NBIN))
            xxdis(:) = 0.D0
            sumdis(:) = 0.D0
            !auxdis(:) = 0.D0
      endif 

      confxx(:)=0.D0 
      conf2xx(:)=0.D0

      set_history(:) = 0.D0 
      
      allocate(temp_xx(1:int(totwalks_in_set/numprocs))) ! keep track of position for all walkers handled on a proc
      allocate(walkers_xix(1:int(totwalks_in_set/numprocs),1:2*NT)) ! keep track of walkers FBM steps in a given set

      !global_corr = 0.D0
      !local_corr(:) = 0.D0

! Loop over disorder configurations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL

      !!! PERFORM SETS OF WALKS WITH NWALKS_IN_SET walkers !!!
      set_loop: do iset=1,NSETS 

      !!! BEFORE EACH SET BEGINS IN EACH PROC WE MUST: 
      set_history(:) = 0.D0 ! 1. clear the set_history 
      set_history(0) = 0.D0 ! 2. init the starting distribution 
      temp_xx(:) = 0.D0 ! 3. init starting positions for all walkers in a proc 

      
      if (myid==0) then
            if (iset == 1) then
                  call system_clock(tnow,tcount)
                  write(*,'(A,I0,A,I0,A,I0,A,I0,A)') 'dis. set ', iset,'/', NSETS,&
                   ' with ', totwalks_in_set, '/', totconf, ' walkers'
            else
                  tlast=tnow
                  call system_clock(tnow)
                  write(*,'(A,I0,A,I0,A,I0,A,F0.3,A)') 'dis. conf. ', iset,' of ',NSETS,&
                ' (took ',(tnow-tlast)/(60*tcount),' minutes and ',mod(tnow-tlast,60*tcount)/(tcount*1.D0),' seconds)'
            end if 
      endif

#else
      !disorder_loop: do iwalker=1,totwalks_in_set
!         print *, 'dis. conf.', iconf
#endif 

      !!! FOR EVERY WALKER IN A PROC, CALCULATE AND SAVE ITS FBM NOISE !!! 
      i = 1
      init_xix_loop: do iwalker=myid+1,totwalks_in_set,numprocs

            iconf = iset*totwalks_in_set - (totwalks_in_set-1) + (iwalker-1) ! find unique iconf number for each walker 
            call gkissinit(IRINIT+iconf-1) 
            call corvec(xix,2*NT,M+1,GAMMA)  ! create x-increments (correlated Gaussian random numbers)
            xix(:) = xix(:)*STEPSIG
            walkers_xix(i,:) = xix(:)
            i = i + 1 

      end do init_xix_loop

      !!! TIME LOOP !!!  
      time_loop: do it=1, NT
 
            ! for each proc at this time step, we want to find it's walkers contribution to our set history dist 
            !temp_history(:) = 0.D0 

            walks_on_proc: do iwalker=1,int(totwalks_in_set/numprocs)

                  !! 1. calculate current bin for iwalker !!
                  ibin=nint( temp_xx(iwalker)*NBIN/LBY2 ) 

                  grad = 0.D0 ! reset grad in case we went past the wall in soft wall case
		      window_mu = 0.D0
		      bin_mu = 0.D0
		      bin_mu2 = 0.D0
		      window_covar = 0.D0

                  xix(it) = walkers_xix(iwalker,it) ! store current walkers FBM step here 

                  !! 2. calculate gradient for iwalker from full-time set_history !!
                  if( (GRAD_FORM .eq. 'ASYM') .and. (GRAD_TEST .eq. 'NONE') ) then

                        if( abs(ibin).le.(NBIN-GRAD_DX) ) then ! if within exclusive (-NBIN,NBIN), calculate gradient with current bin and bin in the direction particle wants to move 
                                    
                              if ( xix(it) .gt. 0.D0 ) then ! if FBM noise points rightward, calculalte gradient with current and immediate right bin
                                    grad = ( set_history(ibin+GRAD_DX) - set_history(ibin) ) / (GRAD_DX*(LBY2/NBIN))
            
                              else if ( xix(it) .lt. 0.D0 ) then ! if FBM points leftward, calculate gradient with current and immediate left bin
                                    grad = ( set_history(ibin) - set_history(ibin-GRAD_DX) ) / (GRAD_DX*(LBY2/NBIN))
                                          
                              else ! else, xix=0, so we flip a coin 
                                    grad = ( set_history(ibin+GRAD_DX) - set_history(ibin) ) / (GRAD_DX*(LBY2/NBIN))
                                    if (rkiss05() < 0.5D0) then
                                          grad = ( set_history(ibin) - set_history(ibin-GRAD_DX) ) / (GRAD_DX*(LBY2/NBIN))
                                    end if 
                              end if 


                        else ! we're in region where calculating gradient towards wall walks off distribution
                              if (ibin .gt. 0) then ! we're at right wall
                                    grad = ( set_history(ibin) - set_history(ibin-GRAD_DX) ) / (GRAD_DX*(LBY2/NBIN)) ! take gradient to the left 
                              else ! we're at left wall 
                                    grad = ( set_history(ibin+GRAD_DX) - set_history(ibin) ) / (GRAD_DX*(LBY2/NBIN)) ! take gradient to the right 
                              end if 
                        end if 

                  end if 


                  if(GRAD_FORM .eq. 'SYMM') then

                        if( abs(ibin).le.(NBIN-GRAD_DX) ) then ! if within exclusive (-NBIN,NBIN), calculate gradient with current bin and bin in the direction particle wants to move 
                              grad = ( set_history(ibin+GRAD_DX) - set_history(ibin-GRAD_DX) ) / (GRAD_DX*2.D0*LBY2/NBIN)
                              
                        else ! we're in region where calculating gradient towards wall walks off distribution
                              if (ibin .gt. 0) then ! we're at right wall
                                    grad = ( set_history(ibin) - set_history(ibin-GRAD_DX) ) / (GRAD_DX*(LBY2/NBIN)) ! take gradient to the left 
                              else ! we're at left wall 
                                    grad = ( set_history(ibin+GRAD_DX) - set_history(ibin) ) / (GRAD_DX*(LBY2/NBIN)) ! take gradient to the right 
                              end if 
                        end if 

                  end if 


                  !! 3. calculate the step contribution from gradient term !! 
                  if (FORCE_TYPE .eq. 'NONLIN') then
                        force_step = nonlin_factor*STEPSIG*tanh(force_weight*grad/nonlin_factor) 

                  else ! FORCE_TYPE = 'linear'
                        force_step = force_weight*grad 
                  endif 

                  !! 4. find and save walker's new position !! 
                  temp_xx(iwalker) = temp_xx(iwalker) + xix(it) + force_step

                  if (WALL .eq. 'SOFT') then
                        !temp_xx(iwalker) = temp_xx(iwalker) + wall_force*exp(-lambda*(xx(it-1)+LBY2)) - wall_force*exp(lambda*(xx(it-1)-LBY2)) 

                  else ! WALL .eq. 'HARD'
                        if ( abs(temp_xx(iwalker)).gt.LBY2 ) then ! stopping boundaries
                              temp_xx(iwalker) = temp_xx(iwalker) - xix(it) - force_step  ! undo the step just taken 
                              print *, 'Walker at t=', it, ' reached the wall.'
                        endif 
                  end if

                  confxx(it)=confxx(it) + temp_xx(iwalker)
                  conf2xx(it)=conf2xx(it) + temp_xx(iwalker)*temp_xx(iwalker)
                  

                  if (WRITE_OUTPUT) then
                        if (myid==0) then
 
                   !write(*,'(I0.1, A, F0.7,A,F0.10,A,F0.3,A,F0.3,A,I0.1,A,I0.1,A,I0.1,A)')  it, ' : ',&
                   !xx(it), ' = ',&
                   !xx(it-1), ' + ',xix(it), ' + ', force_step, " [",&
                   !int(set_history(ibin-1)), ",", int(set_history(ibin)), &
                   !",", int(set_history(ibin+1)), "]"
 
                        endif
                   endif


                  !! 5. find walker's new bin from calculated position & increment necessary distributions !! 
                  !ibin=nint( temp_xx(iwalker)*NBIN/LBY2 )

                  !temp_history(ibin)=temp_history(ibin)+1.0D0 

                  if ( (ibin.ge.-NBIN) .and. (ibin.le.NBIN)) then
                        if( (it.ge.NTSTART) .and. (it.le.NTEND) .and. WRITEDISTRIB) then
                              xxdis(ibin) = xxdis(ibin) + 1.0D0 ! build steady-state distribution 
                        end if 
                  end if
                  
            end do walks_on_proc


            !!! AFTER TIME LOOP, WE MUST UPDATE FULL-TIME SET_HISTORY DISTRIBUTION !!! 

            !! parent receives temp histories from each process, and sums with current set_history !! 
            if (myid .eq. 0) then 

                  do iwalker=1,int(totwalks_in_set/numprocs) 
                        ibin=nint( temp_xx(iwalker)*NBIN/LBY2 )
                        set_history(ibin) = set_history(ibin) + 1.D0
                  end do

                  !set_history(:) = set_history(:) + temp_history(:)
                  do id=1,numprocs-1
                        CALL MPI_RECV(temp_xx,int(totwalks_in_set/numprocs), MPI_DOUBLE_PRECISION, id, &
                                                                        1, MPI_COMM_WORLD,status, ierr)
                        do iwalker=1,int(totwalks_in_set/numprocs) 
                              ibin=nint( temp_xx(iwalker)*NBIN/LBY2 )
                              set_history(ibin) = set_history(ibin) + 1.D0
                        end do
                        !set_history(:) = set_history(:) + temp_history(:)
                  end do 

            !! children send temp histories to parent !! 
            else if (myid .ne. 0) then
                  call MPI_SEND(temp_xx,int(totwalks_in_set/numprocs),MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
            end if 

            !! wait until all machines have sent their temp histories to parent, and parent has parsed them all !!
            call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

            !! parent sends back new set_history to each of the procs !! 
            if (myid .eq. 0) then
                  do id=1, numprocs-1
                        call MPI_SEND(set_history,2*NBIN+1,MPI_DOUBLE_PRECISION,id,2,MPI_COMM_WORLD,ierr)
                  end do 
            !! children receive new set_histroy from parent !! 
            else if (myid .ne. 0) then
                  call MPI_RECV(set_history,2*NBIN+1,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,status,ierr)
            end if 

      end do time_loop

end do set_loop  

! Now collect and analyze data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PARALLEL

      !! for children, send data arrays to parent !! 
     if (myid.ne.0) then                                                 
         call MPI_SEND(confxx,NT,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
         call MPI_SEND(conf2xx,NT,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)

         if(WRITEDISTRIB) then
            call MPI_SEND(xxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
         end if

      !! for parent: 
      else
            !! add it's data to the sum !! 
            sumxx(:)=confxx(:)
            sum2xx(:)=conf2xx(:)
            sumdis(:)=xxdis(:)

         !! receive data from each child, and add it to the data sum !! 
            do id=1,numprocs-1                                                  
                  call MPI_RECV(confxx,NT,MPI_DOUBLE_PRECISION,id,3,MPI_COMM_WORLD,status,ierr)
                  call MPI_RECV(conf2xx,NT,MPI_DOUBLE_PRECISION,id,4,MPI_COMM_WORLD,status,ierr)
                  sumxx(:)=sumxx(:)+confxx(:) 
                  sum2xx(:)=sum2xx(:)+conf2xx(:) 

                  if(WRITEDISTRIB) then
                        call MPI_RECV(xxdis,2*NBIN+1,MPI_DOUBLE_PRECISION,id,5,MPI_COMM_WORLD,status,ierr)
                        sumdis(:)=sumdis(:)+xxdis(:)
                  end if 
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
        write(2,*) 'FORCE_TYPE', FORCE_TYPE
        write(2,*) 'FORCE_TEST', FORCE_TEST
        write(2,*) 'GRAD_FORM', GRAD_FORM
        write(2,*) 'GRAD_DX', GRAD_DX
        write(2,*) 'GRAD_TEST', GRAD_TEST
       ! write(2,*) '<sum(grad)*sum(xix)>', sum_global/totconf
	  write(2,*) 'WINDOW_WIDTH', 2*WINDOW + 1
        write (2,*)'=================================='
        write(2,*) '   time         <r>         <r^2>'      
        it=1
        do while(it.le.NT)  
          Write(2,'(1X,I8,6(2X,E13.6))')  it, sumxx(it)/totconf, sum2xx(it)/totconf
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
        do ibin=-NBIN,NBIN 
          x= (L/2*ibin)/NBIN
          PP= (sumdis(ibin)*NBIN)/(L/2*totconf*(NTEND-NTSTART+1))
          PPsym= ( (0.5D0*sumdis(ibin)+0.5D0*sumdis(-ibin))*NBIN)/(L/2*totconf*(NTEND-NTSTART+1))
          Write(2,'(1X,I9,8(2X,E13.6))') ibin, x, x/L, PP, PP*L, PPsym, L/2-x 
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