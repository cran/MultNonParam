!     subroutine tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats,names)
      subroutine tskmsurv(ntot,rt,delta,ngrp,group,nstats,stats)
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
      integer ii,jj,this,mxr,mir,next,yt,dt,ct,nstats
      double precision pooleds,lastpooleds,stats(max(nstats,1)),sa
!     character (*) :: names
#ifdef DEBUGME
      character*58 fmt
      fmt="(a4,1x,i2,1x,2(a6,1x,  (i3,1x)),a8,  (1x,f5.3),a6,1x,f5.3)"
      write(fmt(22:23),'(i2)') ngrp
      write(fmt(36:37),'(i2)') ngrp
      write(6,*) "|",fmt,"|"
      write(6,*) "Entered tskmsurv",nstats," ngrp",ngrp
      write(6,*) "group",(group(ii),ii=1,ntot)
      write(6,*) "delta",(delta(ii),ii=1,ntot)
#endif
      if(nstats.eq.0) then
         nstats=4
      else
!        names="KS AD CV1CV2"
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
#ifdef DEBUGME
         write(6,*) "In  tskmsurv this",this,"mxr",mxr,"s",s
#endif
         do while(this.le.mxr)
#ifdef DEBUGME
            write(6,*) "Top of do loop y(1)",y(1),"y(2)",y(2),"yt",yt
#endif
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
#ifdef DEBUGME
            write(6,*) "dv",dv,"dt",dt
#endif
            do jj=1,ngrp
               if(dv(jj).gt.0) s(jj)=s(jj)*(1.0d0-dv(jj)/dble(y(jj)))
               y(jj)=y(jj)-dv(jj)-cv(jj)
            end do
            lastpooleds=pooleds
            if(dt.gt.0) then
               pooleds=pooleds*(1.0d0-dt/dble(yt))
#ifdef DEBUGME
               write(6,*) "Reducing pooled survival estimate"
#endif
            end if
            yt=yt-dt-ct
            stats(1)=max(stats(1),abs(s(1)-s(2)))
#ifdef DEBUGME
            write(6,*) "stats(1)",stats(1),"s(1)",s(1),"s(2)",s(2)
            if((pooleds.eq.0.0d0).and.(s(1).ne.s(2)).and.(lastpooleds.gt.pooleds)) write(6,*) "Warning a pooled",pooleds,lastpooleds
            if(pooleds.eq.1.0d0) write(6,*) "Warning b"
#endif
            sa=(s(1)*grpsz(1)+s(2)*grpsz(2))/dble(grpsz(1)+grpsz(2))
            if((s(1).ne.s(2)).and.(lastpooleds.gt.pooleds)) stats(2)=&
               stats(2)+(lastpooleds-pooleds)*abs(s(1)-s(2))/(sa*(1.0d0-sa))
            stats(3)=stats(3)+(lastpooleds-pooleds)*abs(s(1)-s(2))
            stats(4)=stats(4)+(lastpooleds-pooleds)*abs(s(1)-s(2))**2
#ifdef DEBUGME
            write(6,fmt) "this",this, "events",(dv(jj),jj=1,ngrp), "censor",(cv(jj),jj=1,ngrp), &
              "survival",(s(jj),jj=1,ngrp),"pooled",pooleds
            write(6,*) "y(1)",y(1),"y(2)",y(2),"yt",yt
            write(6,*) "Intermediate stats",(stats(ii),ii=1,nstats)
#endif
            this=next
         end do!while
#ifdef DEBUGME
         write(6,*) "Finished loop"
#endif
         deallocate(s,y,dv,cv,grpsz)
      end if
#ifdef DEBUGME
      write(6,*) "Leaving tskmsurv"
#endif
      return
      end
