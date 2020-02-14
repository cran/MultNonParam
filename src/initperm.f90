      subroutine initperm(ntot,ngrp,group,perm,efg) 
! efg=0 ok
!    =1 too many groups
      integer ntot,ngrp,group(ntot),perm(ntot),efg
      integer ii,jj,kk
      integer,allocatable,dimension(:)::grplabs,count
      logical found
      allocate(grplabs(ngrp),count(ngrp))
      jj=1
      grplabs(1)=group(1)
      count(1)=1
      ii=1
!     write(6,*) "ntot",ntot
      do while((jj.le.ngrp).and.(ii.lt.ntot))
!        write(6,*) "ntot in loop",ntot," ii ",ii," check ",(ii.le.ntot)
         ii=ii+1
         found=.false.
         do kk=1,jj
            if(group(ii).eq.grplabs(kk)) then
               count(jj)=count(jj)+1
               found=.true.
!              write(6,*) "ii=",ii," indicator ",group(ii)," adding to group ",kk," label ",grplabs(kk)
            end if
         end do
         if(.not.found) then
!           write(6,*) "ii=",ii," indicator ",group(ii)," forming group ",jj+1
            jj=jj+1
            if(jj.le.ngrp) then
               grplabs(jj)=group(ii)
               count(jj)=1
            end if
         end if
      end do
      if(jj.gt.ngrp) then
         efg=1
      else
         if(jj.ne.2) then
            efg=2
         else
!           write(6,*) "count",count(1),count(2)
            if(grplabs(1).gt.grplabs(2)) then
               ii=grplabs(2);grplabs(2)=grplabs(1);grplabs(1)=ii
               ii=count(2);count(2)=count(1);count(1)=ii
            end if
            do ii=1,count(1)
               perm(ii)=grplabs(1)
            end do
            do ii=1,count(2)
               perm(ii+count(1))=grplabs(2)
            end do
         end if
      end if
      return
      end
