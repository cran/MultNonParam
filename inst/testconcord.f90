     program test
implicit none
     integer(kind=8) uu,tot,myuu
     integer nn,ss
     logical cdf
     double precision dc,dtot,dcc
     cdf=.false.
     write(6,*) "nn"
     read(5,*) nn
     if(nn.gt.0) then
        tot=0
        dtot=0.0d0
        do ss=0,nn*(nn-1)/2
           myuu=uu(nn,ss,cdf)
           call dconcordant(ss,nn,dcc)
           tot=tot+myuu
           dtot=dtot+dcc
           write(6,*) myuu,tot,dcc,dtot
        end do
        call pconcordant(nn*(nn-1)/2-1,nn,dc)
        write(6,*) dc
        call qconcordant(0.025d0,nn,ss)
        write(6,*) "ss=",ss
!       tot=0
        do ss=0,nn*(nn-1)/2
           myuu=uu(nn,ss,.true.)
!          tot=tot+myuu
           write(6,*) myuu,tot
        end do
     else
        nn=-nn
        call qconcordant(0.05d0,nn,ss)
        write(6,*) "ss=",ss
     end if
     end
