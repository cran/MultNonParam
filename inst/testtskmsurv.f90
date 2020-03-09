      program testkm
implicit none
      integer ngrp,ntot,ii,efg,jj,nstats
      double precision,allocatable,dimension(:)::time,cens
      integer,allocatable,dimension(:)::group,delta,rt
      double precision stats(4),pvs(4)
      character (len=12) :: names
      ngrp=2
      write(6,*) "ntot"
      read(5,*) ntot
!     write(6,*) "In program tester about to allocate"
      allocate(time(ntot),cens(ntot),group(ntot),delta(ntot),rt(ntot))
      do jj=1,2
         call loadup(ntot,time,cens,group)
         call fixdata(ntot,time,cens,delta)
!        write(6,*) "In program tester finished loop, time",(time(ii),ii=1,ntot)
         call rankem(time,rt,ntot)
!        write(6,*) "Before tskmsurv rt",(rt(ii),ii=1,ntot)
         call tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats)
!        write(6,*) "Returned from tskmsurv, nstats",nstats
!        write(6,*) "stats",(stats(ii),ii=1,nstats)
         call tskmsurvpv(ntot,rt,delta,ngrp,group,nstats,pvs,efg) 
         write(6,*) "pvs",(pvs(ii),ii=1,nstats)
      end do
      deallocate(time,cens,group,delta,rt)
      end

      subroutine loadup(ntot,time,cens,group)
      integer ntot
      double precision time(ntot),cens(ntot)
      integer group(ntot)
      integer ii
      do ii=1,ntot
!        write(6,*) "In program tester about to call two random uniforms"
         time(ii)=-log(rand(0))
         cens(ii)=-2*log(rand(0))
!        write(6,*) "In program tester finishd call to two unforms"
         group(ii)=(2*ii-1)/ntot+1
!        write(6,*) "In program tester set group "
      end do
      return 
      end
      subroutine fixdata(ntot,time,cens,delta)
      integer ntot
      double precision time(ntot),cens(ntot)
      integer delta(ntot)
      integer ii
      do ii=1,ntot
         delta(ii)=0
!        write(6,*) "In program tester set delta"
         if(time(ii)<=cens(ii)) delta(ii)=1
!        write(6,*) "In program tester set delta"
         time(ii)=min(time(ii),cens(ii))
!        write(6,*) "In program tester set time"
      end do
      return 
      end
