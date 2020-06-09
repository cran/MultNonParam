      program testkm
implicit none
      integer ngrp,ntot,ii,efg,jj,nstats,cnt
      double precision,allocatable,dimension(:)::time,cens,stats,statsmat
      integer,allocatable,dimension(:)::group,delta,rt
      double precision pvs(4)
!     character (len=12) :: names
      ngrp=2
      write(6,*) "ntot"
      read(5,*) ntot
!     write(6,*) "In program tester about to allocate"
      jj=max(ntot,1)
      allocate(time(jj),cens(jj),group(jj),delta(jj),rt(jj))
      if(ntot.lt.0) then
         call readme(ntot,time,delta,group)
         deallocate(time,cens,group,delta,rt)
         allocate(time(ntot),cens(ntot),group(ntot),delta(ntot),rt(ntot))
         call readme(ntot,time,delta,group)
      else
         call loadup(ntot,time,cens,group)
         call fixdata(ntot,time,cens,delta)
      end if
      write(6,*) "Data generated"
!     write(6,*) "In program tester finished loop, time",(time(ii),ii=1,ntot)
      call rankem(time,rt,ntot)
!     write(6,*) "Before tskmsurv rt",(rt(ii),ii=1,ntot)
      allocate(stats(4))
      call tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats)
      deallocate(stats)
      write(6,*) "Returned from tskmsurv, nstats",nstats
      allocate(stats(nstats))
      do ii=1,nstats
         stats(ii)=0.0d0
      end do
      call tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats)
      write(6,*) "stats",(stats(ii),ii=1,nstats)
      cnt=0
      allocate(statsmat(1))
      write(6,*) "cnt before first call to tskmsurvpv",cnt
      call tskmsurvpv(ntot,rt,delta,ngrp,group,nstats,pvs,cnt,statsmat,efg) 
      write(6,*) "cnt after first call to tskmsurvpv",cnt
      deallocate(statsmat)
      allocate(statsmat(nstats*cnt))
      write(6,*) "cnt before second call to tskmsurvpv",cnt
      call tskmsurvpv(ntot,rt,delta,ngrp,group,nstats,pvs,cnt,statsmat,efg)
      write(6,*) "cnt after second call to tskmsurvpv",cnt
      open(53,file="temp")
      do jj=1,cnt
         write(53,*) (statsmat((jj-1)*nstats+ii),ii=1,nstats)
      end do!jj
      close(53)
      deallocate(statsmat)
      write(6,*) "pvs",(pvs(ii),ii=1,nstats)
      deallocate(time,cens,group,delta,rt)
      end

      subroutine loadup(ntot,time,cens,group)
      integer ntot
      double precision time(ntot),cens(ntot)
      integer group(ntot)
      integer ii
      do ii=1,ntot
!        write(6,*) "In program tester about to call two random uniforms"
         group(ii)=(2*ii-1)/ntot+1
         time(ii)=-log(rand(0))+group(ii)
         cens(ii)=-2*log(rand(0))
!        write(6,*) "In program tester finishd call to two unforms"
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
      subroutine readme(ntot,time,delta,grp)
      integer ntot,ii
      double precision time(max(ntot,1))
      integer delta(max(ntot,1)),grp(max(ntot,1))
      open(52,file="testinput",status="old",position="rewind")
      if(ntot<0) then
         ii=0
         ntot=0
         do while(ii.eq.0)
            ntot=ntot+1
            read(52,*,iostat=ii) time(1),delta(1),grp(1)
!           write(6,*) "ii",ii
         end do
         ntot=ntot-1
!        write(6,*) "ntot",ntot
      else
         write(6,*) "Reading existing file",ntot
         do ii=1,ntot
!           write(6,*) "ii",ii
            read(52,*) time(ii),delta(ii),grp(ii)
!           write(6,*) "ii",ii
         end do
      end if
      close(52)
      return
      end
