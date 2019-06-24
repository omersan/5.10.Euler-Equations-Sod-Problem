!-----------------------------------------------------------------------------!
!WENO Solver for Euler equations (one-dimensional)
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: July 2015
!-----------------------------------------------------------------------------!

program euler
implicit none
integer::nx,ns,nt
real*8 ::dt,tm,dx,ds,x,t
integer::i,k
real*8,allocatable::u(:,:,:)


!reading input file
open(7,file='input_sol.txt')
read(7,*)nx 	!number of intervals in x
read(7,*)ns 	!number of record in t
read(7,*)dt     !time step 
read(7,*)tm	    !maximum time 
close(7)

!spatial grid size
dx = 1.0d0/dfloat(nx)

!time interval for storing data
ds = tm/dfloat(ns)

!max number of time for numerical solution
nt = nint(tm/dt)

!numerical solution with RK3 + CRWENO schemes
allocate(u(1:nx,1:3,0:ns))
call numerical(nx,ns,nt,dx,dt,u)

!write solutions in Tecplot format

!write time variable solutions
open(3, file="euler_sol_time.plt")
write(3,*) 'variables ="x","rho"'
  
do k=0,ns
  	t = dfloat(k)*ds
	write(3,'(a16,i8,a10,f10.4,a3)')'zone f=point i=',nx,',t="time',t,'"'
		
	do i=1,nx
		x = -0.5d0*dx + dfloat(i)*dx
        write(3,*) x, u(i,1,k)
    end do

end do
close(3)

!write 3D view solutions
open(4, file="euler_sol_3Dtx.plt")
write(4,*) 'variables ="t","x","rho"'
write(4,*)'zone f=point i=',ns+1,',j=',nx

do i=1,nx
  
 		x = -0.5d0*dx + dfloat(i)*dx
        
		do k=0,ns
  		t = dfloat(k)*ds		

        write(4,*) t, x, u(i,1,k)
        end do

end do
close(4)

open(5, file="euler_sol_3Dxt.plt")
write(5,*) 'variables ="x","t","rho"'
write(5,*)'zone f=point i=',nx,',j=',ns+1
     
do k=0,ns
  	t = dfloat(k)*ds
    
		do i=1,nx
 		x = -0.5d0*dx + dfloat(i)*dx
    
        write(5,*) x, t, u(i,1,k)
        end do

end do
close(5)


open(7, file="euler_sol_final.plt")
write(7,*) 'variables ="x","rho"'
do i=1,nx
x = -0.5d0*dx + dfloat(i)*dx
write(7,*) x,u(i,1,ns)
end do
        
end


!-----------------------------------------------------------------------------!
!compute numerical solutions
!	* 3rd-order Runge-Kutta for temporal 
!	* 5th-order WENO Scheme for spatial
!-----------------------------------------------------------------------------!
subroutine numerical(nx,ns,nt,dx,dt,u)
implicit none
integer::nx,ns,nt,i,j,k,freq,m
real*8 ::dx,dt,gamma
real*8 ::un(1:nx,1:3),ut(1:nx,1:3),r(1:nx,1:3)
real*8 ::u(1:nx,1:3,0:ns)

!u : stored solution 
!un: numerical solution at time n


!initial condition: Sod problem
call ic(nx,dx,gamma,un)

    
!initial time recording
do m=1,3 
do i=1,nx
u(i,m,0)=un(i,m)
end do
end do

!time integration (with RK3)
k=0 !record index
freq = int(nt/ns) !record frequency
    
do j=1,nt

print*, j

	call rhs(nx,dx,gamma,un,r)
    
	do m=1,3 
	do i=1,nx
    ut(i,m) = un(i,m) + dt*r(i,m)
    end do
    end do

	call rhs(nx,dx,gamma,ut,r)

	do m=1,3 
	do i=1,nx
    ut(i,m) = 0.75d0*un(i,m) +0.25d0*ut(i,m) + 0.25d0*dt*r(i,m)
    end do
    end do

	call rhs(nx,dx,gamma,ut,r)

	do m=1,3 
	do i=1,nx
    un(i,m) = 1.0d0/3.0d0*un(i,m) +2.0d0/3.0d0*ut(i,m) + 2.0d0/3.0d0*dt*r(i,m)
    end do
    end do

	!record data
    if(mod(j,freq).eq.0) then
    k=k+1
    
    do m=1,3 
    do i=1,nx
    u(i,m,k)=un(i,m)
    end do
    end do
    
    end if
        
end do


return
end


!-----------------------------------------------------------------------------------!
!Initial conditions and problem definition
!-----------------------------------------------------------------------------------!
subroutine ic(nx,dx,gamma,q)
implicit none
integer::nx,i
real*8 ::dx,x,r,u,p,e
real*8 ::rhoL,rhoR,uL,uR,pL,pR,x0,gamma
real*8 ::q(1:nx,1:3)

!Sod problem

    rhoL=1.0d0
    uL=0.0d0
    pL=1.0d0	

    rhoR=0.125d0
    uR=0.0d0
    pR=0.1d0

    x0    = 0.5d0
    gamma = 1.4d0
   
!construction initial conditions for conserved variables 


do i=1,nx
  
    x =-0.5d0*dx + dfloat(i)*dx 
  	if(x.gt.x0) then
        r=rhoR
		u=uR
    	p=pR
    else
		r=rhoL
		u=uL
    	p=pL
    end if
    	
	e=p/(r*(gamma-1.0d0))+0.5d0*u*u
        
    !conservative variables 
	q(i,1)=r
	q(i,2)=r*u
	q(i,3)=r*e
end do



return
end


!-----------------------------------------------------------------------------!
!Compute rhs for numerical solution of inviscid Burger equation
!  r = -u*u' 
!
!we are using conservative flux form: r = -(u*u/2)'
!-----------------------------------------------------------------------------!
subroutine rhs(nx,dx,gamma,u,r)
implicit none
integer::nx,i,m
real*8 ::dx,gamma
real*8 ::u(1:nx,1:3),r(1:nx,1:3)
real*8,allocatable ::uL(:,:),uR(:,:)
real*8,allocatable ::fL(:,:),fR(:,:)
real*8,allocatable ::f(:,:)

allocate(uL(1:nx+1,1:3))
allocate(uR(1:nx+1,1:3))
allocate(fL(1:nx+1,1:3))
allocate(fR(1:nx+1,1:3))
allocate(f(1:nx+1,1:3))


!compute upwind reconstruction for conserved variable u
call weno5L(nx,u,uL)

!compute downwind reconstruction for conserved variable u
call weno5R(nx,u,uR)

!compute fluxes
call fluxes(nx,gamma,uR,fR)
call fluxes(nx,gamma,uL,fL)

!compute Riemann solver
!call Rusanov(nx,gamma,uR,uL,fR,fL,f)
call HLLC(nx,gamma,uR,uL,fR,fL,f)

!compute RHS
do m=1,3
do i=1,nx	
	r(i,m) = -(f(i+1,m)-f(i,m))/dx
end do
end do

return
end


!---------------------------------------------------------------------------!
!Riemann solver: Rusanov
!---------------------------------------------------------------------------!
subroutine Rusanov(nx,gamma,uR,uL,fR,fL,f)
implicit none
integer::nx,i,m
real*8::gamma
real*8,dimension(1:nx+1,1:3)::uR,uL,fR,fL,f
real*8,dimension(1:nx+1)::cc

!compute wavespeed (Jacobian = df/du)
!call wavespeed(nx,gamma,u,cc)
call wavespeed2(nx,gamma,uR,uL,cc)

!Riemann solver: Rusanov
do m=1,3
do i=1,nx+1
f(i,m) = 0.5d0*(fR(i,m)+fL(i,m)) - 0.5d0*cc(i)*(uR(i,m)-uL(i,m))
end do
end do


return
end

!---------------------------------------------------------------------------!
!Riemann solver: HLLC
!---------------------------------------------------------------------------!
subroutine HLLC(nx,gamma,uR,uL,fR,fL,f)
implicit none
integer::nx,i,m
real*8::gamma,gm
real*8,dimension(1:nx+1,1:3)::uR,uL,fR,fL,f
real*8::rhLL,uuLL,eeLL,ppLL,rhRR,uuRR,eeRR,ppRR,aaLL,aaRR
real*8::SL,SR,Sp,PLR,Ds(3)

gm=gamma-1.0d0

!Compute D
Ds(1) = 0.0d0
Ds(2) = 1.0d0
       
       
do i=1,nx+1
  	
	!Left state:
	rhLL = uL(i,1)
	uuLL = uL(i,2)/rhLL	
	eeLL = uL(i,3)/rhLL
    ppLL = gm*(eeLL*rhLL - 0.5d0*rhLL*(uuLL*uuLL))
    aaLL = dsqrt(dabs(gamma*ppLL/rhLL))
 	
  
    !Right state:
	rhRR = uR(i,1)
	uuRR = uR(i,2)/rhRR
	eeRR = uR(i,3)/rhRR
    ppRR = gm*(eeRR*rhRR - 0.5d0*rhRR*(uuRR*uuRR))
    aaRR = dsqrt(dabs(gamma*ppRR/rhRR))
 

	!Compute SL and SR
		SL = min(uuLL,uuRR) - max(aaLL,aaRR)
		SR = max(uuLL,uuRR) + max(aaLL,aaRR)

    !Compute compound speed
    	Sp = (ppRR - ppLL + rhLL*uuLL*(SL-uuLL) - rhRR*uuRR*(SR-uuRR)) &
            /(rhLL*(SL-uuLL) - rhRR*(SR-uuRR)) !never get zero

    !Compute compound pressure
    	PLR= 0.5d0*(ppLL + ppRR + rhLL*(SL-uuLL)*(Sp-uuLL) &
                                + rhRR*(SR-uuRR)*(Sp-uuRR))
    !Compute D
    	Ds(3) = Sp
        
	!compute HLLC flux in x-direction
	if(SL.ge.0.0d0) then
		do m=1,3  
    	f(i,m)=fL(i,m)     
		end do	
	else if (SR.le.0.0d0) then
		do m=1,3  
    	f(i,m)=fR(i,m)     
		end do	
	else if (Sp.ge.0.0d0 .and. SL.le.0.0d0) then 
		do m=1,3  
        f(i,m)=(Sp*(SL*uL(i,m)-fL(i,m))+SL*PLR*Ds(m))/(SL-Sp)   	   
		end do	
    else if (Sp.le.0.0d0 .and. SR.ge.0.0d0) then
      	do m=1,3  
    	f(i,m)=(Sp*(SR*uR(i,m)-fR(i,m))+SR*PLR*Ds(m))/(SR-Sp)    
		end do	      
	end if
                     
end do	



return
end


!---------------------------------------------------------------------------!
!compute wave speed
!---------------------------------------------------------------------------!
subroutine wavespeed(nx,gamma,q,ps)
implicit none
integer::nx,i
real*8::gamma,l1,l2,l3
real*8::q(1:nx,3),ps(1:nx+1),rad(1:nx+1)

!Spectral radius of Jacobian
do i=1,nx
l1=dabs(q(i,2)/q(i,1))
l2=dabs(q(i,2)/q(i,1) + dsqrt(gamma*((gamma-1.0d0)*(q(i,3)-0.5d0*q(i,2)*q(i,2)/q(i,1)))/q(i,1)))
l3=dabs(q(i,2)/q(i,1) - dsqrt(gamma*((gamma-1.0d0)*(q(i,3)-0.5d0*q(i,2)*q(i,2)/q(i,1)))/q(i,1)))
rad(i) = max(l1,l2,l3)
end do


!Propagation speed
do i=2,nx
ps(i) = max(rad(i),rad(i-1))
end do
ps(1) = ps(2)
ps(nx+1) = ps(nx)

return
end


!---------------------------------------------------------------------------!
!compute wavespeed based on Roe avarage
!---------------------------------------------------------------------------!
subroutine wavespeed2(nx,gamma,uR,uL,ps)
implicit none
integer::nx,i
real*8::gamma,gm,uu,hh,aa,alpha
real*8::uL(1:nx+1,1:3),uR(1:nx+1,1:3),ps(1:nx+1)
real*8::rhLL,uuLL,eeLL,ppLL,hhLL,rhRR,uuRR,eeRR,ppRR,hhRR

gm=gamma-1.0d0

do i=1,nx+1

	!Left and right states:
	rhLL = uL(i,1)
	uuLL = uL(i,2)/rhLL	
	eeLL = uL(i,3)/rhLL
    ppLL = gm*(eeLL*rhLL - 0.5d0*rhLL*(uuLL*uuLL))
    hhLL = eeLL + ppLL/rhLL


	rhRR = uR(i,1)
	uuRR = uR(i,2)/rhRR
	eeRR = uR(i,3)/rhRR
    ppRR = gm*(eeRR*rhRR - 0.5d0*rhRR*(uuRR*uuRR))
    hhRR = eeRR + ppRR/rhRR

	alpha = 1.0d0/(dsqrt(rhLL) + dsqrt(rhRR))
    
	!Roe averages
	uu = (dsqrt(rhLL)*uuLL + dsqrt(rhRR)*uuRR)*alpha
	hh = (dsqrt(rhLL)*hhLL + dsqrt(rhRR)*hhRR)*alpha
	aa = dsqrt(dabs(gm*(hh - 0.5d0*(uu*uu))))
   
	!characteristic speed for Rusanov flux
	ps(i) = aa + dabs(uu)
  
end do


return
end


!---------------------------------------------------------------------------!
!interface flux reconstruction formula
!computing fluxes from conserved quantities
!---------------------------------------------------------------------------!
subroutine fluxes(nx,gamma,q,f)
implicit none
integer::nx,i
real*8::gamma
real*8::q(1:nx+1,3),f(1:nx+1,3)

do i=1,nx+1
f(i,1) = q(i,2)
f(i,2) = q(i,2)*q(i,2)/q(i,1) + (gamma-1.0d0)*(q(i,3)-0.5d0*q(i,2)*q(i,2)/q(i,1))
f(i,3) = q(i,2)*q(i,3)/q(i,1) + (gamma-1.0d0)*q(i,2)/q(i,1)*(q(i,3)-0.5d0*q(i,2)*q(i,2)/q(i,1))
end do

return
end

!-----------------------------------------------------------------------------!
!WENO5 reconstruction for upwind direction (positive & left to right)
!u(i): solution value at nodes i; i=1,2,...,N
!f(j): recontructed value at nodes j=i-1/2; j=1,...,N+1  
!periodic boundary condition
!-----------------------------------------------------------------------------!
subroutine weno5L(n,u,f)
implicit none
integer::n
real*8 ::u(1:n,1:3),f(1:n+1,1:3)
integer::i,m
real*8 ::a,b,c,d,e,w

do m=1,3

i=0
  	a = u(i+3,m)
  	b = u(i+2,m)
  	c = u(i+1,m)
  	d = u(i+1,m)
  	e = u(i+2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w

i=1
  	a = u(i+1,m)
  	b = u(i,m)
  	c = u(i,m)
  	d = u(i+1,m)
  	e = u(i+2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w

i=2
  	a = u(i-1,m)
  	b = u(i-1,m)
  	c = u(i,m)
  	d = u(i+1,m)
  	e = u(i+2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w
       
do i=3,n-2
  	a = u(i-2,m)
  	b = u(i-1,m)
  	c = u(i,m)
  	d = u(i+1,m)
  	e = u(i+2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w
end do

i=n-1
  	a = u(i-2,m)
  	b = u(i-1,m)
  	c = u(i,m)
  	d = u(i+1,m)
  	e = u(i+1,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w
    
i=n
  	a = u(i-2,m)
  	b = u(i-1,m)
  	c = u(i,m)
  	d = u(i,m)
  	e = u(i-1,m)
  	call weno5(a,b,c,d,e,w)
  	f(i+1,m) = w

end do

return
end

!-----------------------------------------------------------------------------!
!WENO5 reconstruction for downwind direction (negative & right to left)
!u(i): solution value at nodes i; i=1,2,...,N
!f(j): recontructed value at nodes j=i-1/2; j=1,2,...,N+1 
!periodic boundary condition 
!-----------------------------------------------------------------------------!
subroutine weno5R(n,u,f)
implicit none
integer::n
real*8 ::u(1:n,1:3),f(1:n+1,1:3)
integer::i,m
real*8 ::a,b,c,d,e,w


do m=1,3
  
i=1
  	a = u(i+2,m)
  	b = u(i+1,m)
  	c = u(i,m)
  	d = u(i,m)
  	e = u(i+1,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w

i=2
  	a = u(i+2,m)
  	b = u(i+1,m)
  	c = u(i,m)
  	d = u(i-1,m)
  	e = u(i-1,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w
        
do i=3,n-2
  	a = u(i+2,m)
  	b = u(i+1,m)
  	c = u(i,m)
  	d = u(i-1,m)
  	e = u(i-2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w
end do

i=n-1
  	a = u(i+1,m)
  	b = u(i+1,m)
  	c = u(i,m)
  	d = u(i-1,m)
  	e = u(i-2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w

i=n
  	a = u(i-1,m)
  	b = u(i,m)
  	c = u(i,m)
  	d = u(i-1,m)
  	e = u(i-2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w

i=n+1
  	a = u(i-3,m)
  	b = u(i-2,m)
  	c = u(i-1,m)
  	d = u(i-1,m)
  	e = u(i-2,m)
  	call weno5(a,b,c,d,e,w)
  	f(i,m) = w

end do
  
return
end

!----------------------------------------------------------------------------------!
!WENO5 
!----------------------------------------------------------------------------------!
subroutine weno5(a,b,c,d,e,f)
implicit none
real*8 ::a,b,c,d,e,f
real*8 ::q1,q2,q3
real*8 ::s1,s2,s3
real*8 ::a1,a2,a3
real*8 ::eps

q1 = a/3.0d0 - 7.0d0/6.0d0*b + 11.0d0/6.0d0*c
q2 =-b/6.0d0 + 5.0d0/6.0d0*c + d/3.0d0
q3 = c/3.0d0 + 5.0d0/6.0d0*d - e/6.0d0

s1 = 13.0d0/12.0d0*(a-2.0d0*b+c)**2 + 0.25d0*(a-4.0d0*b+3.0d0*c)**2
s2 = 13.0d0/12.0d0*(b-2.0d0*c+d)**2 + 0.25d0*(d-b)**2
s3 = 13.0d0/12.0d0*(c-2.0d0*d+e)**2 + 0.25d0*(3.0d0*c-4.0d0*d+e)**2

!Jiang-Shu estimator
eps = 1.0d-6
a1 = 1.0d-1/(eps+s1)**2
a2 = 6.0d-1/(eps+s2)**2
a3 = 3.0d-1/(eps+s3)**2

!Shen-Zha estimator
!eps = 1.0d-20
!a1 = 1.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s1))**2)
!a2 = 6.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s2))**2)
!a3 = 3.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s3))**2)

f = (a1*q1 + a2*q2 + a3*q3)/(a1 + a2 + a3)

return
end subroutine



