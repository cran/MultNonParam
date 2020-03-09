      subroutine rankem(time,rt,nobs)
implicit none
      integer nobs,ii,count
      integer rt(nobs)
      double precision time(nobs),mint,maxt,next
      logical done
!     write(6,*) "On input to rankem nobs",nobs
!     write(6,*) "On input to rankem time",(time(ii),ii=1,nobs)
      mint=time(1)
      maxt=time(1)
!     write(6,*) "Set maxt to ",maxt
!     write(6,*) "In ramkem before first loop"
      do ii=1,nobs
!        write(6,*) "ii=",ii
         mint=min(time(ii),mint)
         maxt=max(time(ii),maxt)
      end do
!     write(6,*) "In ramkem finished first loop, mint",mint," maxt ",maxt
      done=.false.
!     write(6,*) "Just set next to ",next
      count=0
      do while(.not.done)
         next=max(maxt,0.0d0)+2.0d0
         count=count+1
         do ii=1,nobs
            if(time(ii).gt.mint) next=min(time(ii),next)
         end do
         do ii=1,nobs
            if(time(ii).eq.mint) then 
               rt(ii)=count
!              write(6,*) "Setting rank ",count," in slot ",ii
            end if
         end do
         mint=next
!        write(6,*) "next",next
         done=(next.gt.maxt)
      end do
      return
      end
