   subroutine probest(r,M,N,grpv,ngrp,gn,strv,ustr,nstr,ym,xm,zmat,delta,b,Vf11)
implicit none
   integer r,M,N,grpv(N),strv(N),nstr,ustr(nstr),ngrp,gn(ngrp)
   double precision ym(N,r),xm(N,M),delta(r),b(r),Vf11(r,r)
   double precision,allocatable,dimension(:,:,:):: su
   double precision,allocatable,dimension(:,:):: sut1,G
   double precision,allocatable,dimension(:):: sut2
   logical zmat(N,r)
   double precision,allocatable,dimension(:,:,:)::u1,ut1,u2
   double precision,allocatable,dimension(:,:)::ut2
   integer,dimension(:),allocatable::vv,tt
   double precision yd
   integer ii,kk,hh,sstr,jj,jjp,i1,i2,i3,td,mm,nj,njp,njk,njpk,countn
   allocate(u1(N,N,r),ut1(N,N,M),ut2(N,N),u2(N,N,r),tt(N),vv(N))
   allocate(su(2,N,r),sut1(N,M),sut2(N),G(N,2*r+M+1))
   do ii=1,N
      if(grpv(ii).eq.gn(2)) then
         tt(ii)=1
      else 
         tt(ii)=-1
      end if
      do jj=1,r
         su(1,ii,jj)=0.0d0
         su(2,ii,jj)=0.0d0
      end do!jj
      do mm=1,M
         sut1(ii,mm)=0.0d0
      end do!mm
      sut2(ii)=0.0d0
   end do!ii=1,N
   do hh=1,nstr
      do kk=1,r
         call matchv(sstr,vv,ustr(hh),zmat,strv,N,r,kk)
         do jj=1,sstr
            njk=countn(N,vv(jj),strv,ustr(hh),grpv,zmat,r,kk)
            do jjp=1,sstr
               njpk=countn(N,vv(jjp),strv,ustr(hh),grpv,zmat,r,kk)
               yd=ym(vv(jj),kk)-ym(vv(jjp),kk)
               td=tt(vv(jj))-tt(vv(jjp))
               if(grpv(vv(jj)).eq.gn(1)) yd=yd+delta(kk)
               if(grpv(vv(jjp)).eq.gn(1)) yd=yd-delta(kk)
               i1=0
               i2=0
               i3=0
               if((td.gt.0).and.(yd.gt.0.0d0)) i1=1
               if((td.lt.0).and.(yd.lt.0.0d0)) i1=1
               if(td.ne.0) i2=1
               if(yd.eq.0.0d0) i3=1
               u1(vv(jj),vv(jjp),kk)=(i1+0.5d0*i2*i3)/(njk+njpk+1.0d0)
               u2(vv(jj),vv(jjp),kk)=i2/(njk+njpk+1.0d0)
               if(vv(jj).ne.vv(jjp)) then
!                 write(6,*) "jj=",vv(jj),"jjp=",vv(jjp),"u1=",u1(vv(jj),vv(jjp),kk),"i1=",i1
                  su(1,vv(jj),kk)=su(1,vv(jj),kk)+u1(vv(jj),vv(jjp),kk)
                  su(2,vv(jj),kk)=su(2,vv(jj),kk)+u2(vv(jj),vv(jjp),kk)
               end if
            end do!jjp
         end do!jj
      end do!kk=1,r
      call matchv(sstr,vv,ustr(hh),zmat,strv,N,r,0)
      do jj=1,sstr
         nj=countn(N,vv(jj),strv,ustr(hh),grpv,zmat,r,0)
         do jjp=1,sstr
            td=tt(vv(jj))-tt(vv(jjp))
            njp=countn(N,vv(jjp),strv,ustr(hh),grpv,zmat,r,0)
            ut2(vv(jj),vv(jjp))=0.5d0*abs(td)/float(nj+njp)
            do mm=1,M
               ut1(vv(jj),vv(jjp),mm)=0.5d0*td*(xm(vv(jj),mm)-xm(vv(jjp),mm))/float(nj+njp)
               if(vv(jj).ne.vv(jjp)) sut1(vv(jj),mm)=sut1(vv(jj),mm)+ut1(vv(jj),vv(jjp),mm)
            end do!mm
            sut2(vv(jj))=sut2(vv(jj))+ut2(vv(jj),vv(jjp))
         end do!jjp
      end do!jj
   end do!hh=1,nstr
   do ii=1,N
      sut2(ii)=sut2(ii)/float(N-1)
      do mm=1,M
         sut1(ii,mm)=sut1(ii,mm)/float(N-1)
         G(ii,r+mm)=sut1(ii,mm)
      end do!mm
      do kk=1,r
         su(1,ii,kk)=su(1,ii,kk)/float(N-1)
         su(2,ii,kk)=su(2,ii,kk)/float(N-1)
         G(ii,kk)=su(1,ii,kk)
         G(ii,kk+r+M)=su(2,ii,kk)
      end do!kk
      G(ii,2*r+M+1)=sut2(ii)
   end do!ii
!  write(6,*) "sut2",(sut2(ii),ii=1,N)
   call finish(M,r,G,N,b,Vf11)
   deallocate(u1,ut1,ut2,u2,tt,vv,su,sut1,sut2,G)
   return
   end
   subroutine matchv(sstr,vv,hh,zmat,strv,N,r,kk)
implicit none
   integer r,N,sstr,strv(N),vv(N),hh,kk
   logical zmat(N,r)
   integer ii
   sstr=0
   do ii=1,N
      if(kk.gt.0) then
         if(zmat(ii,kk).and.(strv(ii).eq.hh)) then
            sstr=sstr+1
            vv(sstr)=ii
         end if
      else
         if(strv(ii).eq.hh) then
            sstr=sstr+1
            vv(sstr)=ii
         end if
      end if
   end do
   return
   end
   integer function countn(N,ij,strv,hh,grpv,zmat,r,kk)
   integer N,strv(N),hh,grpv(N),r,kk,ij
   logical zmat(N,r)
   integer ii
   countn=0
   do ii=1,N
      if(kk>0) then
         if((strv(ii).eq.hh).and.(grpv(ii).eq.grpv(ij)).and.zmat(ii,kk)) countn=countn+1
      else
         if((strv(ii).eq.hh).and.(grpv(ii).eq.grpv(ij))) countn=countn+1
      end if
   end do!ii
   return
   end
!##############################################
     subroutine finish(M,r,G,N,b,Vf11)
implicit none
     integer M,r,N
     double precision G(N,2*r+M+1),b(r),Vf11(r,r)
     double precision,allocatable,dimension(:)::Gbar,ff
     double precision,allocatable,dimension(:,:)::VGbar,H,Vf22,Vf21,Vf12
     integer,allocatable,dimension(:)::ipiv
     integer ii,jj,kk,ll
     allocate(Gbar(2*r+M+1),VGbar(2*r+M+1,2*r+M+1),ff(r+M),H(r+M,2*r+M+1))
     do ii=1,2*r+M+1
        Gbar(ii)=0.0d0
        do kk=1,N
           Gbar(ii)=Gbar(ii)+G(kk,ii)
        end do!kk
        Gbar(ii)=Gbar(ii)/N
     end do!ii
     do ii=1,2*r+M+1
        do jj=ii,2*r+M+1
           VGbar(ii,jj)=0.0d0
           do kk=1,N
              VGbar(ii,jj)=VGbar(ii,jj)+ (G(kk,ii)-Gbar(ii))*(G(kk,jj)-Gbar(jj))
           end do!kk
           VGbar(ii,jj)=(VGbar(ii,jj)/(N-1))*(4.0d0/N)
           if(jj.ne.ii) VGbar(jj,ii)=VGbar(ii,jj)
        end do!jj
!       write(6,*) (VGbar(ii,jj),jj=1,2*r+M+1)
     end do!ii
     do ii=1,r
        ff(ii)=Gbar(ii)/Gbar(r+M+ii)
     end do!ii
     if(M.gt.0) then
       do ii=1,M
          ff(ii+r)=Gbar(r+ii)/Gbar(2*r+M+1)
       end do!ii
     end if
     do ii=1,M+r
        do jj=1,2*r+M+1
           H(ii,jj)=0.0d0
        end do!jj
        H(ii,ii)=ff(ii)/Gbar(ii)
     end do!ii
     do ii=1,r
        H(ii,r+M+ii)=-ff(ii)/Gbar(r+M+ii)
     end do!ii
     if(M.gt.0) then
        do ii=1,M
           H(r+ii,2*r+M+1)=-ff(r+ii)/Gbar(2*r+M+1)
        end do 
     end if
     do ii=1,r
        b(ii)=ff(ii)
        do jj=1,r
           Vf11(ii,jj)=0.0d0
           do kk=1,2*r+M+1
              do ll=1,2*r+M+1
                 Vf11(ii,jj)=Vf11(ii,jj)+H(ii,kk)*H(jj,ll)*VGbar(kk,ll)
              end do
           end do
        end do
!       write(6,*) "Vf11",(Vf11(ii,jj),jj=1,r)
     end do
!    write(6,*) "b",(b(ii),ii=1,r)
     if(M.gt.0) then
        allocate(ipiv(M), Vf21(M,r),Vf22(M,M),Vf12(r,M))
        do ii=1,r
           do jj=r+1,r+M
              Vf21(jj-r,ii)=0.0d0
              do kk=1,2*r+M+1
                 do ll=1,2*r+M+1
                    Vf21(jj-r,ii)=Vf21(jj-r,ii)+H(ii,kk)*H(jj,ll)*VGbar(kk,ll)
                 end do
              end do
              Vf12(ii,jj-r)=Vf21(jj-r,ii)
           end do
!          write(6,*) "Vf21",(Vf21(jj,ii),jj=1,M)
        end do
        do ii=r+1,r+M
           ff(ii-r)=ff(ii)
           do jj=r+1,r+M
              Vf22(ii-r,jj-r)=0.0d0
              do kk=1,2*r+M+1
                 do ll=1,2*r+M+1
                    Vf22(ii-r,jj-r)=Vf22(ii-r,jj-r)+H(ii,kk)*H(jj,ll)*VGbar(kk,ll)
                 end do
              end do
           end do
!          write(6,*) "Vf22",(Vf22(ii-r,jj),jj=1,M)
        end do
! Lapack routine to calculate LU decomposition.  Feeds matrix inversion routine.
        call dgetrf(M,M,Vf22,M,ipiv,ii)
        call dgetrs('N',M,r,Vf22,M,ipiv,Vf21,M,ii)
!       write(6,*) "Before ff",(ff(ii),ii=1,M)
        call dgetrs('N',M,1,Vf22,M,ipiv,ff,M,ii)
!       write(6,*) "After ff",(ff(ii),ii=1,M)
        do ii=1,r
           do kk=1,M
              b(ii)=b(ii)-Vf12(ii,kk)*ff(kk)
              do jj=1,r
                 Vf11(ii,jj)=Vf11(ii,jj)-Vf12(ii,kk)*Vf21(kk,jj)
              end do
           end do
        end do
        deallocate(ipiv,Vf21,Vf22,Vf12)
     end if
     deallocate(Gbar,VGbar,ff,H)
     return
     end
