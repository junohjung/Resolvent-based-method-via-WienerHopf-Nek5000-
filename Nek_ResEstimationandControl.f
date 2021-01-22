c-----------------------------------------------------------------------
C
C  USER SPECIFIED ROUTINES:
C
C     - boundary conditions
C     - initial conditions
C     - variable properties
C     - local acceleration for fluid (a)
C     - forcing function for passive scalar (q)
C     - general purpose routine for checking errors etc.


c ============================================================================
c           Modules used for obtainind and using Resolvent-based control 


      include 'SavePerturbations.f' ! Module for saving perturbations to disk
      include 'Load_FLD_To.f'       ! Module for loading flow snapshots to an array memory
      include 'FLD2Force.f'         ! Uses the abobe module to use flow snapshots as an external force.
      include 'Shapes.f'            ! Defines domain "Shapes", used as support for sensors and actuators
      include 'Control.f'           ! Control module based on a control law based on the convolution of a control Kernel and sensor readings
c ============================================================================

c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg

      udiff =0.
      utrans=0.
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg)
      use FLD2Force
      use Control
      use Shapes
      include 'SIZE'
      include 'TOTAL'
      
      include 'NEKUSE'
      real fffx,fffy,fffz,tmp_rand(2)
      logical,save :: firstCall = .true.
      real, save :: rms = 0
      real tmp
      integer , save :: lastISTEP = -1 , icount = 0
      integer ix,iy,iz, eg,e,ijke ,ntot , nloc
      integer currIter,sensRange(2),actRange(2)
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      COMMON /ResMode/ currIter,sensRange,actRange
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      e = gllel(eg)
      
      if (firstCall) then
        if (ix==1 .and. iy==1 .and. iz==1 .and. eg==1 ) then
          write(*,*) 'AAAA Setting Seeds', firstCall
        endif
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        seed=nidd
        CALL RANDOM_SEED(PUT = seed)
c        call random_init(.true., .true.)
        DEALLOCATE(seed)
        firstCall=.false.
      endif
      
      ! Check for current run type, and setups external forces accordling
      if ( currIter  ==1 ) then
          ! Direct run, as in "Resolvent-based optimal estimation of transitional and turbulent flows"
          ! Previous run from an adjoint run are used as external forcing using the "Load_FLD_To" module 
          call FLD2Force_GetF(ix,iy,iz,eg,fffx,fffy,fffz)
          ffx = fffx
          ffy = fffy
          ffz = fffz
      elseif (currIter == 3 .or. currIter == 4  ) then
          ! Generates a pseudo random external force for random and control run
          ! Needs to be replaces by a gaussian random force at some point... 
           tmp = 6.14*cos(x*485. + time*3.78)
           tmp = 6.14*sin(y*357 + tmp*352.75)
           tmp = sin(100.45*x+250.48*y+335.54*time + tmp*352.23)        
           ffx = tmp
           tmp = 6.14*cos(x*985.54 + time)
           tmp = 6.14*sin(y*457.12 + tmp*252)
           tmp = sin(137.58*x+218.9*y+335.75*time + tmp*652.45)
           ffy = tmp
           ffz = 0
c        tmp = bm1(ix,iy,iz,e)
c        call random_number(tmp_rand)
c        ffx= 2*(tmp_rand(1)-0.5)*3/sqrt(dt*bm1(ix,iy,iz,e))
c        ffx= sqrt(-2*log(tmp_rand(1)))*cos(2*pi*tmp_rand(2))/sqrt(dt)
c        call random_number(tmp_rand)
c        ffy= sqrt(-2*log(tmp_rand(1)))*cos(2*pi*tmp_rand(2))/sqrt(dt)
c        ffy= 2*(tmp_rand(1)-0.5)*3/sqrt(dt*bm1(ix,iy,iz,e))
c        call random_number(tmp_rand)
c        ffz= sqrt(-2*log(tmp_rand(1)))*cos(2*pi*tmp_rand(2))/sqrt(dt)
c        ffz= 2*(tmp_rand(1)-0.5)*3/sqrt(dt*bm1(ix,iy,iz,e))
c        if (ix==1 .and. iy==1 .and. iz==1 .and. eg==1 ) then
c          rms = rms+ffz**2.0
c          write(*,*) 'AAAA', dt, ffx,ffy,ffz,' rms ',
c     $                      (rms*dt)/(dt*istep)
c          if (istep==600) call exitt
c        endif
      else
          !For currrun==0, the adjoint run, no external forces are used. Intead sensor shape is used as an initial condition. 
          ffx = 0
          ffy = 0
          ffz = 0        
      endif
      
      if (x>-5  ) then
          !Filters external forces for x > -5, as Pff(x,x') = h(x-5)delta(x-x') is assumed.
          ffx = 0 
          ffy = 0 
          ffz = 0 
      endif

      
      if (currIter == 4) then
          !If on control run, uses the "Control" module to compute the actuation of each actuator.
          do i=1,(actRange(2)-actRange(1)+1) ! loop from 1 to the number of actuators
            call getShape(actRange(1)-1+i,x,y,z,fffx,fffy,fffz) ! get the shape of the i-th actuator
            
            ! Add actuation to the external forces
            ffx = ffx + Control_Actuation(i)*fffx
            ffy = ffy + Control_Actuation(i)*fffy
            ffz = ffz + Control_Actuation(i)*fffz
          enddo
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
      qvol   = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      use FLD2Force
      use Control
      use Shapes
      implicit none
      include 'SIZE'
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! PARAM
      include 'SOLN'            ! V[XYZ]
      include 'GEOM'            ! xm1,ym1
      include 'ADJOINT'
      include 'MASS'

      real,external::glsum
      integer ntot,i,j 
      real pertNorm(lpert) , maxnorm 
      character(len=100) FLDReference 
      character(len=100) fileAdd1 
      real readings(30)
      real,save :: absTol = 0
      real,save ::  relativeTol = 1e-4
      real , save :: exitTime = 1e99
      logical, save :: runCompleted = .false. 
      logical, save :: firstSave = .true. 
      integer currIter,sensRange(2),actRange(2)
      COMMON /ResMode/ currIter,sensRange,actRange
      ntot = nx1*ny1*nz1*nelv
      
      ! Save Mass Matrix for posprocessing
      ! Mass matrix can be usefull for some posprocessing. It is not nedded for 
      ! control, but can be helpful to compute perturbation norms, etc 
      if (NIO==0) write(*,*) "Saving Mass Matrix on istep==0"
      if (istep==0) then
            call exportShapes()
            ifxyo = .true.
            call outpost2(bm1,bm1,
     $                    bm1,bm2,
     $                    bm1,0,'bm1')
            if (currIter==-3) call exitt

      endif
      
      
      ! If on a Control run (currIter=4), update the Control module. 
      if (currIter == 4) then
        ! Initialize control on the first time step.
        if (istep ==0) then
           if (NIO==0) write(*,*) 'Init Control Module...  '
           fileAdd1 = 'ControlLawParams.txt' ! Control parameter file. Includes number of sensors, actuators, delta t and control laws.
           call ControlInit(fileAdd1)
        endif 
        
        write(*,*) 'Update Control Module...  '      
        ! Use the shapes module to obtain current sensor readings. Print them for debug purposes.
        readings=ProjectOnShapes(1)
        if (NIO==0) write(*,*) 'Current Readings ',
     $       time,readings(sensRange(1):sensRange(2))
        call ControlUpdateReadings(
     $       time,readings(sensRange(1):sensRange(2)))
       
       ! Update sensor history on the "Control" module 
        call Control_CalcActuations()
       ! Print the "Control" module  sensor history  for debug.
        call ControlPrintReadingsHist()
       ! Print the "Control" module  current actuations.
        call ControlPrintActuation()
       ! Save the "Control"  module history (actuation and readings). Useful for debug.
        call SaveControlHist()
      endif     
      
      ! If on a Direct run (currIter == 1), use a previous adjoint runs as  external force.       
      if (NIO==0) write(*,*) 'Updating Forces...'
      if ( currIter  == 1) then
            if (istep==0) then ! initialize  the "FLD2ForceInit" module 
                  fileAdd1 = 'ForceFLDLists.txt' ! file containing the adjoint run snapshots (in reverse order)
                  call FLD2ForceInit(fileAdd1)                  
            endif
            ! update the module (reads new files if nedded)
            call  FLD2ForceUpdate()
      endif

      ! Turn on the adjoint solver if on an Adjoint run (curriter==0)
      if (istep ==0) then       
          ifadj = (currIter  == 0)
      endif
      
      ! Save current sensor readings (used latter to obtain control laws).
      call SaveProjUp2Shapes(1,1)



      ! Calculate current perturbation norm (used as a stopping criteria for 
      ! the direct/adjoint/actuator runs. Also useful to evaluate control performance.
      if (NIO==0) write(*,*) 'Computing run Norm...'
      maxnorm = 0
      do j=1,npert
            pertNorm(j) = 0
            do i=1,ntot
              if (xm1(i,1,1,1)>0) then
               pertNorm(j) = pertNorm(j) + 
     c         ( vxp(i,j)**2 + vyp(i,j)**2 )*bm1(i,1,1,1)
     	      endif
            enddo
            pertNorm(j) = glsum(pertNorm(j),1)
            pertNorm(j) = sqrt(pertNorm(j))
            
            maxnorm = max(maxnorm,pertNorm(j))  
      end do
      ! sets run tolerance: run stops when the current norm is 10^-5 (relative norm) of the maximum norm.
      absTol = max(absTol,maxnorm*relativeTol)  
      if (NIO==0) write(*,*) 'AbsTol : ', absTol,
     $   ' rel', relativeTol  ,' pertNorm ', maxnorm
      
      ! Writes current perturbation norm to file
      if (NIO==0) then
        if (istep == 0) then
          open (unit=98, file='runNorm.txt',
     $            status='replace',action='write')
        else
          open (unit=98, file='runNorm.txt',
     $        position="append" ,status='old',action='write')
        endif
        write(98,*) time , maxnorm , exitTime !,runCompleted
        close(98)
      endif
      
      
      ! Uses "SavePerturbation" module to save flow snapshots. The prefix is a function of the current run
      if (NIO==0) write(*,*) 'Save Flow...'
      if ( maxnorm > absTol ) then
        ifxyo = firstSave
        if (currIter==0) call SavePerturbations('a') ! adjoint run
        if (currIter==1) call SavePerturbations('d') ! direct  run
        if (currIter==2) call SavePerturbations('f') ! actuator run (f for force)
        if (currIter==3) call SavePerturbations('r') ! random run
        if (currIter==4) call SavePerturbations('c') ! control run
        firstSave = .false.
      endif
      
      if (NIO==0) write(*,*) 'Checking Exit condtions...'
      ! setup exit flag if perturbation norm is smaller than tolerance. 
      ! If on a direct run also requires that all adjoint files to have been read
      if (currIter==1) then
            runCompleted = .not.  FLD2Force_StillActive() ! set end of run after force is no longer active         
      else
            runCompleted =.true.
      endif
           
      if ( maxnorm < absTol .and. runCompleted ) then
      if (NIO==0) write(*,*) 'Creating runCompled file...'
        open (unit=98, file='runCompleted',
     $            status='replace',action='write')
        write(98,*) 'Done' 
        close(98)
        
        call exitt
      end if
      
      return
      end

c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
c     NOTE ::: This subroutine MAY NOT be called by every process
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      if (JP>0) then
            ux=0.0
            uy=0.0
            uz=0.0
            temp=0.0
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      use Shapes
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer currIter,sensRange(2),actRange(2)
      COMMON /ResMode/ currIter,sensRange,actRange
      ! Uses sensor/actuator shapes as initial conditions for the 
      ! adjoit/actuator run, and null initial conditions otherwise
      if (JP==0) then
            ux = 0.0
            uy = 0.0
            uz = 0.0
      else
        if ( currIter  == 0 .or. currIter  == 2) then
          call getShape(sensRange(1)-1+JP,x,y,z,ux,uy,uz)
        else
            ux = 0.0
            uy = 0.0
            uz = 0.0
        endif
      endif
      
      temp=0
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      use Shapes
      include 'SIZE'
      include 'TOTAL'
      integer :: a
      character(len=100) fileAdd1 
c      integer nio
      integer currIter,sensRange(2),actRange(2)
      COMMON /ResMode/ currIter,sensRange,actRange
      
      ! Reads the run type, as well as sensor and actuator ranges 
      ! i.e. range 1-3 indicates that shapes 1,2 and 3 are used.  
      open (unit=98, file='currRun.txt',
     $                status='old',action='read')
      read(98,*) currIter
      close(98)
      open (unit=98, file='sensRange.txt',
     $                status='old',action='read')
      read(98,*) sensRange
      close(98)
      open (unit=98, file='actRange.txt',
     $                status='old',action='read')
      read(98,*) actRange
      close(98)      
      
      ! Initi "Shapes" module. 'Shapes.txt' provides a parametric list of shapes.
      fileAdd1='shapes.txt'
      call ReadShapesFile(fileAdd1)
      
c     call platform_timer(0) ! not too verbose
c     call platform_timer(1) ! mxm, ping-pong, and all_reduce timer
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'      

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
c      include 'MASS'
      integer currIter,sensRange(2),actRange(2)
      COMMON /ResMode/ currIter,sensRange,actRange
      
      !Save Undeformed Mass Matrix (before loading external fld) for posprocessing.
      if (istep==0) then
            ifxyo = .true.
            call outpost2(bm1,bm1,
     $                    bm1,bm2,
     $                    bm1,0,'bm0')
      endif
      if (currIter==-3) call exitt
      
      return
      end



c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
