       subroutine signtestperm(y,s,n,out)
implicit none
       integer n
       double precision y(n),s(n),signtestone
       integer ii,kk,out
       logical index(n)
       logical done
       double precision testobs,testnew
       do ii=1,n
          index(ii)=(y(ii).gt.0.0d0)
!         write(6,*) "index(ii)",index(ii),"y(ii)",y(ii)
       end do!ii
       testobs=signtestone(index,s,n)
       do ii=1,n
          index(ii)=.False.
       end do!ii
       done=.False.
       out=0
       do while(.not.done)
          kk=-1
          do while(kk.le.0)
             if(index(-kk)) then 
                kk=kk-1
             else
                kk=-kk
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
             testnew=signtestone(index,s,n)
             if(testnew.ge.testobs) out=out+1
!            write(6,*) (index(ii),ii=1,n),testnew,testobs,out
          end if
       end do
       return
       end
       double precision function signtestone(index,s,n)
       integer n,ii
       double precision s(n)
       logical index(n)
       signtestone=0.0d0
       do ii=1,n
          if(index(ii)) signtestone=signtestone+s(ii)
       end do!ii
       return
       end
