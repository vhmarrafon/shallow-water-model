program version_1

implicit none
real,dimension(0:141,0:96)::up,uf,vp,vf,ui,vi
real,dimension(1:140,1:95)::np,nf,ni,uout,vout
real::c1,c2,c3,c4,c5
integer::i,n,j,leap
real,parameter::dx=100000,dy=50000,dt=4000,H=1,f=0.000,r=0.0000001,g=9.81
integer,parameter::xt=140,yt=95,nt=1000
!coriolis f=0.00010284
!leap=0
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
ui=0
vi=0

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
if(n<=4)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #########################################################################################
	do i=1,xt
		do j=1,yt
			nf(i,j)=np(i,j)-c4*(up(i+1,j)-up(i,j))-c5*(vp(i,j+1)-vp(i,j))
		enddo
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     #####################################################################################
	do i=0,xt+1
		do j=0,yt+1
!! CONDICAO DE CONTORNO PARA Y
			if(j==0 .or. j==yt+1 .or. i==0 .or. i==xt+1)then
				vf(i,j)=0
				uf(i,j)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else
				uf(i,j)=up(i,j)+c1*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-c2*(np(i,j)-np(i-1,j))
				vf(i,j)=vp(i,j)-c1*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-c3*(np(i,j)-np(i,j-1))
			endif
		enddo
	enddo

else
	
! #############################       LEAP FROG     #########################################################################
! ###########################################################################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #########################################################################################
!if(leap==1)then	
	do i=1,xt
		do j=1,yt
			nf(i,j)=ni(i,j)-2*c4*(up(i+1,j)-up(i,j))-2*c5*(vp(i,j+1)-vp(i,j))
		enddo
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     #####################################################################################
	do i=0,xt+1
		do j=0,yt+1
!! CONDICAO DE CONTORNO PARA Y
			if(j==0 .or. j==yt+1 .or. i==0 .or. i==xt+1)then
				vf(i,j)=0
				uf(i,j)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else
				uf(i,j)=ui(i,j)+2*c1*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-2*c2*(np(i,j)-np(i-1,j))
				vf(i,j)=vi(i,j)-2*c1*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-2*c3*(np(i,j)-np(i,j-1))
			endif
		enddo
	enddo

endif

	write(20)nf
	!write(20)uf
	!write(20)vf
	ni=np
	ui=up
	vi=vp
	np=nf
	up=uf
	vp=vf
enddo
end
