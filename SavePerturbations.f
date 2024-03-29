
c========================================================
c     AdjEst_Output, called from AdjEst_Init
c
c     Outposts perturbations for adjoint, direct, actuator
c       and control runs.
c========================================================
      subroutine SavePerturbations(prefix)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer i
      real  ,save   :: lastSavedTime=-1e99 
      real IOTIME
c      integer IOSTEP
      character prefix
      logical exportFlag 
      
      IOTIME = param(14)
      exportFlag = .false.
      if (IOSTEP .ne. 0) exportFlag =
     $      exportFlag .or. mod(ISTEP,abs(IOSTEP))==0

      if (IOTIME .ne. 0) exportFlag =
     $      exportFlag .or. time-lastSavedTime>=abs(IOTIME)
c     ------ Writes Output Files --------------------
      if (exportFlag) then
        lastSavedTime = time
        call SavePerturbationsNow(prefix)
      end if
      return
      end subroutine


      subroutine SavePerturbationsSetIPS(i)
      integer :: i
      integer :: ipscalout 
      COMMON /savePerts/ ipscalout
      ipscalout = i
      end subroutine
      

      subroutine SavePerturbationsNow(prefix)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer ,save :: currSaveId=0
      integer ntot,i , j
      CHARACTER(LEN=3)  :: currPert
      character prefix
      logical previfxyo 
      logical,save :: firstSave =.true.  
      real pertNorm(lpert)
      integer :: ipscalout = 0
      COMMON /savePerts/ ipscalout
      
      ntot = nx1*ny1*nz1*nelv
      previfxyo = ifxyo
      ifxyo = currSaveId==0
      currSaveId = currSaveId+1

c     ------ Writes Output Files --------------------

      if (NIO==0)  then
        write(*,*) "AdjEst : Outposting ",
     $              npert, "perturbations"
      end if
      if (ipscalout==0) then
        j = 1+jp
      else
        j = ipscalout
      endif

c     Output perturbation states
      do i = 1,npert
         pertNorm(i) = 0
         do j=1,ntot 
           pertNorm(i) = pertNorm(i) + 
     $         (vxp(j,I)**2. + vyp(j,I)**2.)*bm1(j,1,1,1)
           if (ldim==3)  pertNorm(i) = pertNorm(i) + 
     $         (vzp(j,I)**2.)*bm1(j,1,1,1)
         enddo
         
         pertNorm(i) = glsum(pertNorm(i),1)
         pertNorm(i) = sqrt(pertNorm(i))

         write(currPert,10) prefix, i
   10    FORMAT(A1,I0.2)
         call outpost2(VXP(1,i),VYP(1,i),
     $                 VZP(1,i),PRP(1,i),
     $                 tp,ipscalout,currPert)
        
      end do
      if (NIO==0) then
        if (firstSave) then
          open (unit=98, file='pertNorm.txt',
     $            status='replace',action='write')
          firstSave = .false.
        else
          open (unit=98, file='pertNorm.txt',
     $         position="append" ,status='old',action='write')
        endif
        write(98,*) time , pertNorm(:)
        write(*,*) "Perturbation Norm " ,time , pertNorm(:)
        close(98)
      endif
      ifxyo = previfxyo
      return
      end subroutine
