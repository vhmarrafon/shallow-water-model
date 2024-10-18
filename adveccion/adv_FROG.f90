program version_zero
implicit none
real::ci(100),cf(100),cp(100)
integer::i,n
real,parameter::dx=0.01,dt=0.009,u=1
integer,parameter::xt=100,nt=90

!!!!! NUNCA ADICIONAR TEMPO NA MATRIZ
! CFC cdt/dx<1, dt<dx/c
!write(*,*)1/dx !n de pontos em x = 100
!write(*,*)0.81/dt !n de pontos em t = 90

!PRIMEIRO TEMPO FRON
 cp=0
 cf=0
 ci=0

do i=2,xt-1
	cp(i)=exp(-200*((i*dx-0.25)**2))
enddo

!primeiro tempo
 cp(1)=0
 cp(xt)=0
do i=2,xt-1
	cf(i)=-1*u*(dt/dx)*(cp(i)-cp(i-1))+cp(i)
enddo
	ci=cp
	cp=cf
do n=2,nt-1

write(*,*)cp(xt-1),n
	do i=2,xt-1
		cf(i)=-1*u*((2*dt)/(2*dx))*(cp(i+1)-cp(i-1))+ci(i)
		!write(*,*)cf(i),ci(i),cp(i)

	enddo
!!! sem (i) salva para todos os valores de i
!!! write(10) -> salva em binario
	ci=cp	
	cp=cf
	cp(1)=0
	cp(xt)=0
	cf(i)=0
	cf(xt)=0
	write(11)cf


enddo
end

























