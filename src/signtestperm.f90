       subroutine signtestperm(y,s,n,out,verbose)
implicit none
       integer n
       double precision y(n),s(n),signtestone
       logical verbose
       character star(1)
       integer ii,kk,out,marka,markb
       logical,dimension(:),allocatable::index
       logical done
       double precision testobs,testnew
       allocate(index(n))
       if(verbose) ii=0
       do ii=1,n
          index(ii)=(y(ii).gt.0.0d0)
!         write(6,*) "index(ii)",index(ii),"y(ii)",y(ii)
       end do!ii
       testobs=signtestone(index,s,n,marka)
       do ii=1,n
          index(ii)=.False.
       end do!ii
       done=.False.
       out=0
       do while(.not.done)
          kk=-1
          do while(kk.le.0)
!            write(6,*) "kk",kk,(index(ii),ii=1,n)
             if(kk.ge.(-n)) then
                if(index(-kk)) then 
                   kk=kk-1
                else
                   kk=-kk
                end if
             else
                kk=n+1
             end if
          end do 
          if(kk.gt.n) then
             done=.True.
          else
             do ii=1,kk-1
                index(ii)=.False.
             end do!ii
             index(kk)=.True.
          end if
          if(.not.done) then
             testnew=signtestone(index,s,n,markb)
             if(testnew.ge.testobs) out=out+1
             star=" "
             if(marka.eq.markb) star="*"
!            if(verbose) write(6,*) (index(ii),ii=1,n),testnew,testobs,out,star
          end if
       end do
       deallocate(index)
       return
       end
       double precision function signtestone(index,s,n,mark)
implicit none
       integer n,ii,mark
       double precision s(n)
       logical index(n)
       signtestone=0.0d0
       mark=0
       do ii=1,n
          mark=2*mark
          if(index(ii)) then
             mark=mark+1
             signtestone=signtestone+s(ii)
          end if
       end do!ii
       return
       end
