program version_1

implicit none
real,dimension(140,95)::np,nf,up,uf,vp,vf,uout,vout
real::c1,c2,c3,c4,c5
integer::i,n,j
real,parameter::dx=100000,dy=50000,dt=4000,H=1,f=0.000,r=0.0000001,g=9.81
integer,parameter::xt=140,yt=95,nt=1000
!coriolis f=0.00010284

!!!!! NUNCA ADICIONAR TEMPO NA MATRIZ
! CFC cdt/dx<1, dt<dx/c
!write(*,*)1/dx !n de pontos em x = 100
!write(*,*)0.81/dt !n de pontos em t = 90
 c1=(dt*f)/4
 c2=(dt*g)/dx
 c3=(dt*g)/dy
 c4=(H*dt)/dx
 c5=(H*dt)/dy
up=0
vp=0

!!ORDEM do DO, tempo,x,y
do i=1,xt
	do j=1,yt
		np(i,j)=exp((-1.*((i-0.5*xt)/5.)**2)-1.*(((j-0.5*yt)/5.)**2))
		!write(*,*)np(i,j)
	enddo
enddo



!escrevendo o primeiro timestep
write(20)np
!write(20)up
!write(20)vp


do n=1,nt

	do i=1,xt
		do j=1,yt
!! CONDICAO DE CONTORNO PARA Y
			if(j==1 .or. j==yt)then
				vp(i,j)=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! CONDICOES DE CONTORNO PARA i==1 e i==xt

			elseif(i==xt .and. j/=1 .and. j/=yt )then
			!!!! i+1 = 1
				vf(i,j)=vp(i,j)-c1*(up(i,j)+up(i,j-1)+up(1,j-1)+up(1,j))-c3*(np(i,j)-np(i,j-1))-r*dt*vp(i,j)
				nf(i,j)=np(i,j)-c4*(up(1,j)-up(i,j))-c5*(vp(i,j+1)-vp(i,j))
			elseif(i==1 .and. j/=1 .and. j/=yt )then
			!!!! i-1 = xt
				uf(i,j)=up(i,j)+c1*(vp(i,j)+vp(i,j+1)+vp(xt,j+1)+vp(xt,j))-c2*(np(i,j)-np(xt,j))-r*dt*up(i,j)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else
				uf(i,j)=up(i,j)+c1*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-c2*(np(i,j)-np(i-1,j))-r*dt*up(i,j)
				vf(i,j)=vp(i,j)-c1*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-c3*(np(i,j)-np(i,j-1))-r*dt*vp(i,j)
				nf(i,j)=np(i,j)-c4*(up(i+1,j)-up(i,j))-c5*(vp(i,j+1)-vp(i,j))
				vout(i,j)=(vp(i,j+1)+vp(i,j))/2
				uout(i,j)=(up(i+1,j)+up(i,j))/2
			endif
		enddo
	enddo

	write(20)nf
	!write(20)uf
	!write(20)vf
	np=nf
	up=uf
	vp=vf

enddo
end

























