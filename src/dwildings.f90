
# 1 "wilding.F90"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "wilding.F90"



   double precision recursive function wilding(u1,u2,m1,n1,m2,n2) result(out)

implicit none
   integer u1,u2,m1,n1,m2,n2



   double precision zero, one

   logical cc(4)
! write(6,'(a15,6(1x,a3,i3))') "Calling wilding with", "u1=",u1,"u2=",u2, "m1=",m1,"m2=",m2, "n1=",n1,"n2=",n2
   zero=0
   one=1
   cc(1)=(m1.lt.0).or.(m2.lt.0).or.(n1.lt.0).or.(n2.lt.0)
   cc(2)=(u1.lt.0).or.(u1.gt.(m1*n1)).or.(u2.lt.0).or.(u2.gt.((m1+m2)*(n1+n2)))
   cc(3)=((n1.eq.0).or.(m1.eq.0)).and.(u1.ne.0)
   cc(4)=((n2.eq.0).or.(m2.eq.0)).and.(u2.ne.u1)
   if(cc(1).or.cc(2).or.cc(3).or.cc(4)) then
! write(6,*) "Branch A"
      out=zero
   else
      if((m1.eq.0).and.(m2.eq.0).and.(n1.eq.0).and.(n2.eq.0)) then
! write(6,*) "Branch B"
         out=one
      else
! write(6,*) "Branch C"
         out=(m1*wilding(u1,u2,m1-1,n1,m2,n2)+&
             n1*wilding(u1-m1,u2-m1-m2,m1,n1-1,m2,n2)+&
             m2*wilding(u1,u2,m1,n1,m2-1,n2)+&
             n2*wilding(u1,u2-m1-m2,m1,n1,m2,n2-1))/(m1+m2+n1+n2)
      end if
   end if
   return
   end











     subroutine wildings(u1,u2,m1,n1,m2,n2,out)
implicit none
     integer u1,u2,m1,n1,m2,n2



     double precision out,wilding

     out=wilding(u1,u2,m1,n1,m2,n2)
     return
     end
