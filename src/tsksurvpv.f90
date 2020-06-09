      subroutine tskmsurvpv(nobs,rt,delta,ngrp,group,npv,pvs,cnt,statsmat,efg) 
implicit none
      integer nobs,rt(nobs),delta(nobs),ngrp,group(nobs),newtot,efg,ii,npv,cnt,space
      integer,allocatable,dimension(:)::perm,count
      double precision,allocatable,dimension(:)::stats,stats0
      double precision statsmat(max(npv*cnt,1))
      double precision pvs(max(npv,1))
!     character (*) :: names
      character(20) :: fmt
!     write(6,*) "cnt",cnt,"npv",npv,"max(npv*cnt,1)",max(npv*cnt,1)
      allocate(perm(nobs))
      efg=0
!     write(6,*) "After initperm perm=",(perm(ii),ii=1,nobs)
      if(npv.eq.0) then
         allocate(stats0(1))
!        call tskmsurv(nobs,rt,delta,ngrp,group,npv,stats0,names)
         call tskmsurv(nobs,rt,delta,ngrp,group,npv,stats0)
         cnt=0
         deallocate(stats0)
      else
         space=cnt
         allocate(stats(npv),stats0(npv),count(npv))
         do ii=1,npv
            count(ii)=0
         end do 
!        call tskmsurv(nobs,rt,delta,ngrp,group,npv,stats0,names)
         call tskmsurv(nobs,rt,delta,ngrp,group,npv,stats0)
!        write(6,*) "Observed statistics",(stats0(ii),ii=1,npv)
         call initperm(nobs,ngrp,group,perm,efg) 
         if(efg.eq.0) then
            newtot=nobs
            cnt=0
!           write(6,*) "Reset cnt to zero"
!           fmt='(a4,  (1x,i2))'
!           write(fmt(5:6),'(i2)') nobs
            do while(newtot.gt.0)
!              write(6,fmt) "perm",(perm(ii),ii=1,nobs)
!              call tskmsurv(nobs,rt,delta,ngrp,perm,npv,stats,names)
               call tskmsurv(nobs,rt,delta,ngrp,perm,npv,stats)
               cnt=cnt+1
               do ii=1,npv
                  if(stats(ii).ge.stats0(ii)) count(ii)=count(ii)+1
!                 if(cnt.le.space) write(6,*) "cnt",cnt,"space",space
!                 if(cnt.le.space) statsmat((cnt-1)*npv+ii)=1.0d0
                  if(cnt.le.space) statsmat((cnt-1)*npv+ii)=stats(ii)
!                 write(6,*) (cnt-1)*npv+ii
               end do
!              write(6,*) "New statistics",(stats(ii),ii=1,npv),"Compare with", &
!                 (stats0(ii),ii=1,npv),"count",(count(ii),ii=1,npv),"Overall count",cnt
               call nextp(perm,newtot,1)
            end do
            do ii=1,npv
               pvs(ii)=dble(count(ii))/dble(cnt)
            end do
         end if
         deallocate(stats,stats0,count)
      end if
!     write(6,*) "cnt",cnt,"npv",npv,"space",space
!     write(6,*) "statsmat",(statsmat(ii),ii=1,npv*space)
      deallocate(perm)
      return
      end
