      program testprobest
implicit none
      integer r,M,N,nstr,ngrp
      integer,dimension(:),allocatable::grpv,gn,strv,ustr
      double precision,dimension(:,:),allocatable::ym,xm,Vf11
      double precision,dimension(:),allocatable::delta,b
      logical,dimension(:,:),allocatable::zmat
      integer ii,jj
      r=2
      N=8
      M=0
      nstr=1
      ngrp=2
      allocate(ym(N,r),grpv(N),strv(N),zmat(N,r),delta(r),b(r),Vf11(r,r),xm(N,1))
      open(30,file="temp.csv")
      read(30,*)
      do ii=1,N
         read(30,*) (ym(ii,jj),jj=1,r),grpv(ii),strv(ii)
         do jj=1,r
            zmat(ii,jj)=.FALSE.
         end do
         write(6,*) "ym(ii,jj),jj=1,r),grpv(ii),strv(ii)",(ym(ii,jj),jj=1,r),grpv(ii),strv(ii)
      end do
      do jj=1,r
         delta(jj)=0.0d0
         b(jj)=0.0d0
         do ii=1,r
            Vf11(ii,jj)=0.0d0
         end do
      end do 
      allocate(ustr(nstr),gn(ngrp))
      ustr(1)=1
      gn(1)=2;gn(2)=1
      call probest(r,m,n,grpv,ngrp,gn,strv,ustr,nstr,ym,xm,zmat,delta,b,vF11)
      write(6,*) "b=",(b(jj),jj=1,r)
      end
