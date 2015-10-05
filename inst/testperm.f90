       program testperm
       integer n,out
       double precision,allocatable,dimension(:):: y,s
       n=5
       allocate(y(n),s(n))
       y(1)=1.0
       y(2)=-2.0
       y(3)=3.0
       y(4)=-4.0
       y(5)=5.0
       do out=1,n
          s(out)=out
       end do
       call signtestperm(y,s,n,out)
       end
