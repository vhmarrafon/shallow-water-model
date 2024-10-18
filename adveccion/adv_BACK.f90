program version_zero
implicit none
real::ci(100),cf(100)
integer::i,n
real,parameter::dx=0.01,dt=0.009,u=1
integer,parameter::xt=100,nt=90

!!!!! NUNCA ADICIONAR TEMPO NA MATRIZ
! CFC cdt/dx<1, dt<dx/c
!write(*,*)1/dx !n de pontos em x = 100
!write(*,*)0.81/dt !n de pontos em t = 90

!FRONT


do i=1,xt
	ci(i)=exp(-200*(i*dx-0.25)**2)
enddo

do n=1,nt

	do i=1,xt
		if (i==xt .or. i==1)then
			cf(i)=0	
		else
		cf(i)=-1*u*(dt/dx)*(ci(i)-ci(i-1))+ci(i)
		endif
		!write(*,*)cf(n,i)

	enddo
!!! sem (i) salva para todos os valores de i
!!! write(10) -> salva em binario
	write(10)cf
	ci=cf
enddo
end

























