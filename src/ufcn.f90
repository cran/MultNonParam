     module uucache
     integer, parameter:: j8=selected_int_kind(15)
     integer,dimension(:),allocatable,save::nnvec
     integer (kind=j8) ,dimension(:),allocatable,save ::uuvec
     integer,save:: initnn
     integer(kind=j8),save,dimension(:),allocatable::nfac
     logical,save::ccdf
     end module
     subroutine initmod(nn,cdf)
use uucache
implicit none
     integer nn,ii
     logical cdf
     if(allocated(nnvec)) then
        deallocate(nnvec)
        deallocate(uuvec)
        deallocate(nfac)
     end if
     allocate(nnvec(nn),nfac(nn))
!    write(6,*) "Initializing module"
     nfac(1)=1
     do ii=1,nn
        if(ii.eq.1) nnvec(ii)=1
        if(ii.gt.1) nnvec(ii)=nnvec(ii-1)+ii*(ii-1)/2+1
        if(ii.gt.1) nfac(ii)=nfac(ii-1)*ii
     end do
     allocate(uuvec(nnvec(nn)+1+nn*(nn-1)/2))
     do ii=1,nnvec(nn)+nn*(nn-1)/2+1
        uuvec(ii)=-1
     end do
     initnn=nn
     ccdf=cdf
     return
     end
     recursive function uu(nn,ss,cdf) result(newu)
use uucache
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     logical cdf
     integer nn,ss,ii
     integer(kind=i8) newu,pc,pullcache
     logical,save::init=.false.
     pc=-2
     if(init.and.(ccdf.neqv.cdf)) init=.false.
!    write(6,*) "init",init
     if(init) then
        if(nn.gt.initnn) then 
           init=.false.
!           write(6,*) "Triggering reinitialization nn=",nn,"initnn",initnn
        end if
     end if
     if(.not.init) then
        init=.true.
        call initmod(nn,cdf)
     end if
     if((nn.gt.1).and.(ss.ge.0)) then
        pc=pullcache(nn,ss,cdf)
        if(pc.ge.0)  then
!          write(6,*) "good"
           newu=pc
        else
           newu=0
           do ii=0,nn-1
              newu=newu+uu(nn-1,ss-ii,cdf)
           end do
        end if
     else
        if(ss.lt.0) newu=0
        if(ss.gt.(nn*(nn-1)/2)) then
           if(cdf) then
!             write(6,*) "cdf branch"
              newu=1
           else
              newu=0
           end if
!          write(6,*) "newu,cdf",newu,cdf
        end if
        if((ss.ge.0).and.(ss.le.(nn*(nn-1)/2))) newu=1
     end if
     if(newu.gt.0) call putcache(nn,ss,newu)
!    if(pc.ge.0) write(6,*) "pc,newu",pc,newu,"ss",ss,"nn",nn
     return
     end
     subroutine putcache(nn,ss,val)
use uucache
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     integer (kind=i8) val
     integer nn,ss
!    write(6,*) "Putting ",val," into slot ",nnvec(nn)+ss
     uuvec(nnvec(nn)+ss)=val
     return
     end
     function pullcache(nn,ss,cdf)
use uucache
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     integer nn,ss,ii
     integer(kind=i8) pullcache
     logical cdf
     ii=0
!    write(6,*) "nnvec",(nnvec(ii),ii=1,initnn),"nn=",nn,"ss=",ss
     if((ss.ge.0).and.(ss.le.(nn*(nn-1)/2))) then 
        pullcache=uuvec(nnvec(nn)+ss)
     else
        if(ss.lt.0) pullcache=0
        if(ss.gt.(nn*(nn-1)/2)) then
           if(cdf) then
!             write(6,*) "CDF BRANCH"
              pullcache=nfac(nn)
           else 
              pullcache=0
           end if
        end if
     end if
!    write(6,*) "pullcache",pullcache,"ss=",ss,"nn=",nn,"position",nnvec(nn)+ss
     return
     end
     subroutine dconcordant(ss,nn,dc)
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     integer ss,nn
     double precision dc
     integer(kind=i8) uu,dd
     integer ii
     dd=1
     do ii=2,nn
        dd=dd*ii
     end do
     dc=uu(nn,ss,.false.)/dble(dd)
     return
     end
     subroutine pconcordant(ss,nn,dc)
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     integer ss,nn
     double precision dc
     integer(kind=i8) uu,dd
     integer ii
     dd=1
     do ii=2,nn
        dd=dd*ii
     end do
     dc=uu(nn,ss,.true.)/dble(dd)
     return
     end
     subroutine qconcordant(qq,nn,ss)
! ss is output
implicit none
     integer, parameter:: i8=selected_int_kind(15)
     integer ss,nn,ii
     double precision qq
     integer(kind=i8) cc,nfac,pp,uu
     nfac=1
     do ii=2,nn
        nfac=nfac*ii
     end do
     cc=ceiling(qq*nfac,8)
     pp=0
     ss=0
!    write(6,*) "cc=",cc,"pp=",pp,"ss=",ss,"qq*nfac",ceiling(qq*nfac,8)
     do while(pp.lt.cc)
        ss=ss+1
        pp=uu(nn,ss,.true.)
!       write(6,*) "ss,pp",ss,pp
     end do
     return
     end
