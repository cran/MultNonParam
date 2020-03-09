!     subroutine tskmsurvpv(ntot,rt,delta,ngrp,group,npv,pvs,names,efg) 
      subroutine tskmsurvpv(ntot,rt,delta,ngrp,group,npv,pvs,efg) 
implicit none
      integer ntot,rt(ntot),delta(ntot),ngrp,group(ntot),newtot,efg,ii,npv,cnt
      integer,allocatable,dimension(:)::perm,count
      double precision,allocatable,dimension(:)::stats,stats0
      double precision pvs(max(npv,1))
!     character (*) :: names
      allocate(perm(ntot))
      efg=0
      cnt=0
      call initperm(ntot,ngrp,group,perm,efg) 
!     write(6,*) "After initperm perm=",(perm(ii),ii=1,ntot)
      if(npv.eq.0) then
         allocate(stats0(1))
!        call tskmsurv(ntot,rt,delta,ngrp,group,npv,stats0,names)
         call tskmsurv(ntot,rt,delta,ngrp,group,npv,stats0)
         deallocate(stats0)
      else
         allocate(stats(npv),stats0(npv),count(npv))
         do ii=1,npv
            count(ii)=0
         end do 
!        call tskmsurv(ntot,rt,delta,ngrp,group,npv,stats0,names)
         call tskmsurv(ntot,rt,delta,ngrp,group,npv,stats0)
!        write(6,*) "Observed statistics",(stats0(ii),ii=1,npv)
         if(efg.eq.0) then
            newtot=ntot
            do while(newtot.gt.0)
!              call tskmsurv(ntot,rt,delta,ngrp,perm,npv,stats,names)
               call tskmsurv(ntot,rt,delta,ngrp,perm,npv,stats)
               cnt=cnt+1
               do ii=1,npv
                  if(stats(ii).ge.stats0(ii)) count(ii)=count(ii)+1
               end do
!              write(6,*) "New statistics",(stats(ii),ii=1,npv),"Compare with", (stats0(ii),ii=1,npv),"count",(count(ii),ii=1,npv)
               call nextp(perm,newtot,1)
            end do
            do ii=1,npv
               pvs(ii)=dble(count(ii))/dble(cnt)
            end do
         end if
         deallocate(stats,stats0)
      end if
      return
      end
