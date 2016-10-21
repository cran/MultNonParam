     subroutine cntperms(perm,n,nb,be,tot)
implicit none
     integer nb
     integer be(nb)
     integer n,perm(n)
     integer i,j,start,ngrp
     integer,dimension(:),allocatable::cntgrp
     double precision tot
     ngrp=1
     do i=1,n
!        if(perm(i).le.0) write(6,*) "Group number too low."
!        if(perm(i).ge.100) write(6,*) "Group number too high?"
        ngrp=max(ngrp,perm(i))
     end do
     allocate(cntgrp(ngrp))
     tot=1.0d0
!    write(6,*) "About to cycle through blocks"
     do i=1,nb
        if(i.eq.1) start=1
        if(i.gt.1) start=be(i-1)+1
        tot=tot*gamma(dble(be(i)-start+2))
!       write(6,*) "tot",tot
        do j=1,ngrp
           cntgrp(j)=0
        end do
        do j=start,be(i)
           cntgrp(perm(j))=cntgrp(perm(j))+1
        end do
        do j=1,ngrp
            tot=tot/gamma(dble(cntgrp(j)+1))
        end do
     end do
!     write(6,*) tot
     return
     end






   subroutine labelblock(be,nb,blk,nx)
   integer nb
   integer nx,blk(nx),be(nb),jj,kk,start
   do jj=1,nb
      if(jj.eq.1) start=1
      if(jj.gt.1) start=be(jj-1)+1
      do kk=start,be(jj)
         blk(kk)=jj
      end do
   end do
   return
   end

   
   
   subroutine grpmeans(means,ngrp,gm,nx,x,used,grp)
   integer nx,ngrp,used(ngrp),grp(nx)
   double precision x(nx),means(ngrp),gm
   integer jj
   do jj=1,ngrp
          means(jj)=0.0d0
          used(jj)=0
   end do
   do jj=1,nx
          means(grp(jj))=means(grp(jj))+x(jj)
          used(grp(jj))=used(grp(jj))+1
   end do
!  write(6,*) "Mark a, ngrp=",ngrp
   gm=0.0d0
   do jj=1,ngrp
      gm=gm+means(jj)
      if(used(jj).gt.0) then
         means(jj)=means(jj)/used(jj)
      else
!        write(6,*) "Empty group"
      end if
   end do
!  write(6,*) "Mark a"
   gm=gm/nx
   return
   end



!    subroutine writeanova(means,blkmn,gm,ngrp,bss,ess,nx,aov,nb)
! implicit none
!    integer nx,ngrp,jj,nb,df(3)
!    double precision means(ngrp),gm,bss,ess,aov,blkmn(nb)
!    character*20 fmt
!    fmt='(a11,1x,  (f9.4,1x))'
!   write(fmt(9:10),'(i2)') ngrp
!   write(6,fmt) "group means",(means(jj),jj=1,ngrp)
!   if(nb.gt.1) then
!      write(fmt(9:10),'(i2)') nb
!      write(6,fmt) "block means",(blkmn(jj),jj=1,nb)
!    end if
!    df(1)=ngrp-1
!    df(2)=nx-1-(ngrp-1)-(nb-1)
!      df(3)=nx-1
!   write(6,'(a10,1x,f9.4,1x,a1,1x,f9.4)') "grand mean",gm,"f",aov
!   write(6,'(a47)') "source  df sum of squares   mean squares  f"
!      write(6,'(a5,1x,i5,3(1x,f12.4))') "model",df(1),bss,bss/df(1),aov
!      write(6,'(a5,1x,i5,3(1x,f12.4))') "error",df(2),ess,ess/df(2)
!      write(6,'(a5,1x,i5,3(1x,f12.4))') "total",df(3),bss+ess
!   return
!   end

   
   subroutine chkgrps(ngrp,first,nx,grp)
   integer ngrp,nx,grp(nx)
   logical first
   if(ngrp.le.0) then
      first=.true.
!     write(6,*) "Ping"
      ngrp=0
      do jj=1,nx
!        if(grp(jj).le.0) write(6,*) "Bad group"
         ngrp=max(ngrp,grp(jj))
      end do
!     write(6,*) "Ping"
   else
      first=.false.
   end if
   return
   end
   
   
   double precision function aov2(x,grp,nx,ngrp,nb,be)
implicit none
   integer nx,nb
   integer grp(nx),be(nb)
   double precision x(nx),tss,ess,gm,temp,bss
   integer ngrp,jj,df(3)
   integer,dimension(:),allocatable::used,blk,usedb
   double precision,dimension(:),allocatable::means,meansb
   logical first
!  write(6,*) "Entered aov2, grp",(grp(jj),jj=1,nx)
!  write(6,*) "ngrp",ngrp
   call chkgrps(ngrp,first,nx,grp)
!  write(6,*) "ngrp",ngrp
   allocate(meansb(nb),blk(nx),usedb(nb))!Do this smarter
   allocate(means(ngrp),used(ngrp))
   call grpmeans(means,ngrp,gm,nx,x,used,grp)
   call labelblock(be,nb,blk,nx)!This needs ony be done once
   call grpmeans(meansb,nb,temp,nx,x,usedb,blk)!This needs ony be done once
!  write(6,*) "Done calculating means in aov2"
   tss=0.0d0
   ess=0.0d0
   bss=0.0d0
   df(1)=ngrp-1
   df(2)=nx-1-(ngrp-1)-(nb-1)
   df(3)=nx-1
   do jj=1,nx
      bss=bss+(means(grp(jj))-gm)**2
      ess=ess+(x(jj)-means(grp(jj))-meansb(blk(jj))+gm)**2
!     if(first) write(6,*) "Residuals",x(jj)-means(grp(jj))-meansb(blk(jj))+gm
      tss=tss+(x(jj)-gm)**2
   end do
!  write(6,*) "Done calculating sums of squares in aov2"
   aov2=(bss/df(1))/(ess/df(2))
!  write(6,*) (grp(jj),jj=1,nx)
!   if(first) call writeanova(means,meansb,gm,ngrp,bss,ess,nx,aov2,nb)
!  write(6,*) "About to dealocate in aov2, ngrp",ngrp,"nb",nb
!  write(6,*) (means(jj),jj=1,ngrp)
!  write(6,*) (used(jj),jj=1,ngrp)
   deallocate(means,used)
   deallocate(meansb,blk,usedb)!Do this smarter.
!  write(6,*) "Leaving aov2"
   return
   end

   
   
   
   double precision function aov(x,grp,nx,ngrp)
implicit none
   integer nx
   integer grp(nx)
   double precision x(nx),tss,ess,gm,bss
!  double precision temp(1)
   integer ngrp,jj
   integer,dimension(:),allocatable,save::used
   double precision,dimension(:),allocatable::means
   logical first
   call chkgrps(ngrp,first,nx,grp)
   allocate(means(ngrp),used(ngrp))
   call grpmeans(means,ngrp,gm,nx,x,used,grp)
   tss=0.0d0
   ess=0.0d0
   do jj=1,nx
      ess=ess+(x(jj)-means(grp(jj)))**2
      tss=tss+(x(jj)-gm)**2
   end do
!  write(6,*) "Mark a"
   bss=tss-ess
   aov=(bss/(ngrp-1))/(ess/(nx-ngrp))
!  if(first) call writeanova(means,temp,gm,ngrp,bss,ess,nx,aov,1)
   deallocate(means,used)
   return
   end


   
        subroutine next(perm,n,b)
implicit none
     integer n,perm(n),b
     integer i,j,k,l,mxperm
     logical go
!    write(6,*) "in next n=",n,"b=",b
!    write(6,*) "in next perm=",(perm(i),i=1,n)
     go=.TRUE.
     mxperm=perm(b)
     do i=b,n
        mxperm=max(mxperm,perm(i))!I know that I am looking at the first element twice.
     end do
     i=n-1
     do while(go)
!       write(6,*) "in next i=",i
        if(i.ge.b) then
           if(perm(i)<perm(i+1)) go=.FALSE.
        end if
        if(i.eq.0) go=.FALSE.
        if(go) i=i-1
     end do
! i is now the index of last entry with the entry above it in ascending order
     if(i>0) then
! I think that the next line can be removed
!       j=n
        k=mxperm+1
        l=n+1
        do j=n,i+1,-1
           if((perm(j)>perm(i)).and.(perm(j)<k)) then
              k=perm(j)
              l=j
           end if
        end do
! l is now the index of the smallest entry after i and greater than perm(i)
        j=perm(l)
        perm(l)=perm(i)
        perm(i)=j
!       write(6,*) "Swapping entries",i,"and",l
!       write(6,*) "Intermediate",(perm(j),j=1,n)
        do j=1,(n-i)/2
           l=perm(i+j)
           perm(i+j)=perm(n+1-j)
           perm(n+1-j)=l
        end do
     else
        n=-n
     end if
!    write(6,*) "Before return in next n=",n
     return
     end
     
     subroutine nextb(perm,n,nb,be,first)
implicit none
     integer nb
     integer be(nb)
     integer n,perm(n)
     integer i,j,k,l
     integer,dimension(:),allocatable,save::operm
     logical first
!    write(6,*) "nextb, n=",n,"be",(be(i),i=1,nb)
     if(.not.allocated(operm)) then
        allocate(operm(n))
     end if
     if(first) then
        do i=1,n
           operm(i)=perm(i)
        end do
        first=.false.
     end if
     l=1
     do while(l.gt.0)
        if(l.eq.1) j=1
        if(l.gt.1) j=be(l-1)+1
        k=be(l)
!       write(6,*) "Before next l=",l," k=",k,"j=",j,"perm",(perm(i),i=1,n)
        call next(perm,k,j)
!       write(6,*) "After next k=",k,"l=",l,"perm",(perm(i),i=1,n)
        if(k.lt.0) then
           if(l.lt.nb) then
!             write(6,*) "In nextb branch a"
              k=-k
              l=l+1
              do i=j,k
                 perm(i)=operm(i)
              end do
           else
!             write(6,*) "In nextb branch b"
              nb=-nb
              n=-n
              l=-1!Stop looping through blocks, since we've already completed permuting all blocks
           end if
        else
!          write(6,*) "In nextb flipping l"
           l=-1!Stop looping through blocks, since we have found the next permutation without running of the end.
        end if
     end do
!    write(6,*) "At end of nextb nb",nb
     return
     end

     subroutine aovp(n,grpi,nb,be,x,tot,pv)
implicit none
! Perform permutation analysis of variance.  
! Permute group memberships, perhaps with blocking.
! Report counts of permutations, and those with statistic exceeding
! p-value.
     integer nb
     integer n,ngrp,grpi(n),be(nb),npprt
     integer i
     integer(kind=8) ms(2),rt,cnt(2)
     double precision x(n),aovn,aov,aov2,aovo,tot,remt,pv
!    character*20 fmt
     logical first
! Next line stops uninitialized warning.
     aovn=0.0d0
!    fmt='(a2,i3,a5,  i3,f9.4)'
!    write(fmt(11:12),'(i2)') n
     ngrp=0
!    write(6,*) "ngrp=",ngrp
     aovo=-1.0d0
     cnt(1)=0
     cnt(2)=0
!    open(8,file="temp")
     first=.true.
     i=0
     npprt=1000000
     call cntperms(grpi,n,nb,be,tot)
!     write(6,*) "tot",tot
     call system_clock(ms(1),rt)
!    write(6,*) "Mark 0 ngrp=",ngrp
     do while(n.gt.0)
        i=i+1
!       write(6,*) "ngrp=",ngrp
        if(nb.eq.1) aovn=aov(x,grpi,n,ngrp)
        if(nb.gt.1) aovn=aov2(x,grpi,n,ngrp,nb,be)
!       write(8,fmt) "n=",n,"grpi=",(grpi(i),i=1,abs(n)),aovn
        if(aovo.lt.0.0d0) aovo=aovn
        cnt(1)=cnt(1)+1
        if(aovn.ge.aovo) cnt(2)=cnt(2)+1
!       write(6,*) "About to generate next permutation"
        if(nb.eq.1) call next(grpi,n,nb)
        if(nb.gt.1) call nextb(grpi,n,nb,be,first)
!       write(6,*) "Generated next permutation"
        if(i.eq.((i/npprt)*npprt)) then
!           write(6,'(a30,f20.17)') "Proportion completed",i/tot
           call system_clock(ms(2),rt)
           remt=((ms(2)-ms(1))*tot/(i))/rt
!           write(6,*) rt
!           write(6,*) "Projected completion in days:",remt/(24*60*60)
        end if
     end do
     pv=dble(cnt(2))/dble(cnt(1))
     return
     end

     subroutine betatest(n,x,y,pval)
implicit none
     integer (kind=8) count(2),npprt
     integer,allocatable,dimension(:)::perm,be
!    character*15 fmt
     integer n,i,b
     double precision x(n),y(n),pval,mx,cp,cp0,tot
     allocate(perm(n),be(1))
     npprt=1
     do i=1,3
        npprt=1000*npprt
     end do
     be(1)=n
     b=1
     mx=0.0d0
     do i=1,n
        mx=mx+x(i)
        perm(i)=i
     end do
     call cntperms(perm,n,1,be,tot)
!     write(6,*) "tot=",tot
     mx=mx/n
     do i=1,n
        x(i)=x(i)-mx
     end do 
     cp0=cp(n,x,y,perm)
!    fmt='(a2,i3,a5,  i3)'
!    write(fmt(11:12),'(i2)') n
     count(1)=1
     count(2)=1
     do while(n.gt.0)
!        call projtime(count(2),tot,npprt)
!       write(6,fmt) "n=",n,"perm=",(perm(i),i=1,abs(n))
        call next(perm,n,b)
        if(cp(n,x,y,perm).ge.cp0) count(1)=count(1)+1
        count(2)=count(2)+1
     end do
     pval=dble(count(1))/dble(count(2))
     return
     end
     
     double precision function cp(n,x,y,perm)
implicit none
     integer n,perm(n),i
     double precision x(n),y(n)
     cp=0.0
     do i=1,n
        cp=cp+x(i)*y(perm(i))
     end do
     return
     end
