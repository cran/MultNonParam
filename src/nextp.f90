     subroutine nextp(perm,n,b)
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
