      subroutine tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats,names)
implicit none
! Arguments
! ntot   integer   1      number of observations
! rt     integer   ntot   ranks of event times
! delta  integer   ntot   1 for events, 0 for censored
! ngrp   integer   1      number of groups
! group  integer   ntot   group indicator, should be in 1:ngrp
! nstats integer   1      number of statistics to calculate
! stats  double    nstats statistics
! names  character nstats names
      integer ntot,rt(ntot),delta(ntot),ngrp,group(ntot)
      double precision,allocatable,dimension(:):: s
      integer,allocatable,dimension(:):: y,dv,cv,grpsz
      integer ii,jj,dd,this,mxr,mir,next,yt,dt,ct,nstats
      double precision pooleds,lastpooleds,stats(max(nstats,1)),sa
      character*3 names(max(nstats,1))
!     character*58 fmt
!     fmt="(a4,1x,i2,1x,2(a6,1x,  (i2,1x)),a8,  (1x,f5.3),a6,1x,f5.3)"
!     write(fmt(22:23),'(i2)') ngrp
!     write(fmt(36:37),'(i2)') ngrp
!     write(6,*) "|",fmt,"|"
!     write(6,*) "Entered tskmsurv",nstats
!     write(6,*) "group",(group(ii),ii=1,ntot)
      if(nstats.eq.0) then
         nstats=4
      else
         names(1)="KS ";names(2)="AD ";names(3)="CV1";names(4)="CV2"
         allocate(s(ngrp),y(ngrp),dv(ngrp),cv(ngrp),grpsz(ngrp))
         do jj=1,ngrp
            y(jj)=0;s(jj)=1.0d0
         end do
         yt=0; pooleds=1.0d0
         mxr=rt(1); mir=rt(1)
         do ii=1,ntot
            y(group(ii))=y(group(ii))+1
            if(rt(ii).lt.mir) mir=rt(ii)
            if(rt(ii).gt.mxr) mxr=rt(ii)
         end do
         do jj=1,ngrp
            grpsz(jj)=y(jj)
            yt=yt+y(jj)
         end do
         this=mir
         do jj=1,nstats
            stats(jj)=0.0d0
         end do
!        write(6,*) "In  tskmsurv this",this,"mxr",mxr
         do while(this.le.mxr)
!           write(6,*) "Top of do loop y(1)",y(1),"y(2)",y(2),"yt",yt
            do jj=1,ngrp
               dv(jj)=0; cv(jj)=0
            end do
            dt=0; ct=0
            next=mxr+1
            do ii=1,ntot
               if((rt(ii).gt.this).and.(rt(ii).lt.next)) next=rt(ii)
               if(rt(ii).eq.this) then
                  dv(group(ii))=dv(group(ii))+delta(ii)
                  cv(group(ii))=cv(group(ii))+1-delta(ii)
                  dt=dt+delta(ii)
                  ct=ct+1-delta(ii)
               end if
            end do
            do jj=1,ngrp
               if(dv(jj).gt.0) s(jj)=s(jj)*(1.0d0-dv(jj)/dble(y(jj)))
               y(jj)=y(jj)-dv(jj)-cv(jj)
            end do
            lastpooleds=pooleds
            if(dt.gt.0) pooleds=pooleds*(1.0d0-dt/dble(yt))
            yt=yt-dt-ct
            stats(1)=max(stats(1),abs(s(1)-s(2)))
!           write(6,*) "stats(1)",stats(1)
!           if((pooleds.eq.0.0d0).and.(s(1).ne.s(2)).and.(lastpooleds.gt.pooleds)) write(6,*) "Warning a pooled",pooleds,lastpooleds
!           if(pooleds.eq.1.0d0) write(6,*) "Warning b"
            sa=(s(1)*grpsz(1)+s(2)*grpsz(2))/dble(grpsz(1)+grpsz(2))
            if((s(1).ne.s(2)).and.(lastpooleds.gt.pooleds)) stats(2)=&
               stats(2)+(lastpooleds-pooleds)*abs(s(1)-s(2))/(sa*(1.0d0-sa))
            stats(3)=stats(3)+(lastpooleds-pooleds)*abs(s(1)-s(2))
            stats(4)=stats(4)+(lastpooleds-pooleds)*abs(s(1)-s(2))**2
!           write(6,fmt) "this",this, "events",(dv(jj),jj=1,ngrp), "censor",(cv(jj),jj=1,ngrp), &
!             "survival",(s(jj),jj=1,ngrp),"pooled",pooleds
!           write(6,*) "y(1)",y(1),"y(2)",y(2),"yt",yt
!           write(6,*) "Intermediate stats",(stats(ii),ii=1,nstats)
            this=next
         end do!while
!        write(6,*) "Finished loop"
         deallocate(s,y,dv,cv,grpsz)
      end if
!     write(6,*) "Leaving tskmsurv"
      return
      end
