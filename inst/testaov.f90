       program testaov
implicit none
       integer n,out,ii,nb
       double precision pv
       double precision,allocatable,dimension(:):: y
       integer,allocatable,dimension(:):: grp,be
       write(6,*) "Number of data points"
       read(5,*) n
       open(30,file="lastinput")
       write(30,*) n," #Number of data poitns"
       allocate(y(n),grp(n))
       write(6,*) "Input data points (y, group)"
       do ii=1,n
          read(5,*) y(ii),grp(ii)
          write(30,*) y(ii),grp(ii)," #Data"
       end do
       write(6,*) "How many blocks?"
       read(5,*) nb
       allocate(be(nb))
       write(30,*) nb, " #Number of blocks"
       write(6,*) "Input block endpoints"
       read(5,*) (be(ii),ii=1,nb)
       write(30,*) (be(ii),ii=1,nb), " #Block ends"
       call aovp(n,grp,nb,be,y,ii,pv)
       end
