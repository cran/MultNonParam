       program testperms
       integer nobj,ii
       integer,allocatable,dimension(:):: perm
       write(6,*) "nobj?"
       read(5,*) nobj
       allocate(perm(nobj))
       do ii=1,nobj
          perm(ii)=0
          if(ii.gt.(nobj/2)) perm(ii)=1
       end do
       do while(nobj.gt.0)
          write(6,*) (perm(ii),ii=1,nobj)
          call nextp(perm,nobj,1)
       end do
       end
