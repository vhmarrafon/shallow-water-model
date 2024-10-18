program VERSION_4
implicit none
real,dimension(1:140,1:95)::up,uf,vp,vf,ui,vi
real,dimension(1:139,1:94)::np,nf,ni,uout,vout
real::c1,c2,c3,c4,c5,fb
integer::i,n,j
real,parameter::dx=100000,dy=50000,dt=4000,H=1,f=0.000,r=0.0000001,g=9.81,beta=2.42
integer,parameter::xt=140,yt=95,nt=4000
!coriolis f=0.00010284
!!!!! NUNCA ADICIONAR TEMPO NA MATRIZ
! CFC cdt/dx<1, dt<dx/c
!write(*,*)1/dx !n de pontos em x = 100
!write(*,*)0.81/dt !n de pontos em t = 90
 c1=(dt)/4
 c2=(dt*g)/dx
 c3=(dt*g)/dy
 c4=(H*dt)/dx
 c5=(H*dt)/dy
up=0
vp=0
ui=0
vi=0

!!ORDEM do DO, tempo,x,y
do i=1,xt-1
	do j=1,yt-1
		np(i,j)=2*exp((-1.*((i-0.5*xt)/5.)**2)-1.*(((j-0.5*yt)/5.)**2))
		!write(*,*)np(i,j)
	enddo
enddo



!escrevendo o primeiro timestep
!write(20)np
!write(20)up
!write(20)vp


do n=1,nt
if(n<=4)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #########################################################################################
	do i=1,xt-1
		do j=1,yt-1
			nf(i,j)=np(i,j)-c4*(up(i+1,j)-up(i,j))-c5*(vp(i,j+1)-vp(i,j))
		enddo
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     #####################################################################################
	do i=1,xt
		do j=1,yt
!! CONDICAO DE CONTORNO PARA Y
				fb=(j-0.5*yt)*beta*dy*(0.00000000001)
!! CONDICAO DE CONTORNO
			if(j==1 .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)+c1*fb*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-c2*(np(i,j)-np(i-1,j))-r*up(i,j)*dt
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)+c1*fb*(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))-c2*(np(i,1)-np(i-1,1))-r*up(i,j)*dt


			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-c1*fb*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-c3*(np(i,j)-np(i,j-1))-r*vp(i,j)*dt
				uf(i,j)=ui(i,j)+c1*fb*(vp(i,j)+vp(i,1)+vp(xt,1)+vp(xt,j))-c2*(np(i,j)-np(xt-1,j))-r*up(i,j)*dt
			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-c1*fb*(up(i,j)+up(i,j-1)+up(1,j-1)+up(i+1,j))-c3*(np(i,j)-np(i,j-1))-r*vp(i,j)*dt
				uf(i,j)=ui(i,j)+c1*fb*(vp(i,j)+vp(i,1)+vp(xt,1)+vp(xt,j))-c2*(np(1,j)-np(i-1,j))-r*up(i,j)*dt
			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else
				!write(*,*)fb
				uf(i,j)=up(i,j)+c1*fb*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-c2*(np(i,j)-np(i-1,j))-r*up(i,j)*dt
				vf(i,j)=vp(i,j)-c1*fb*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-c3*(np(i,j)-np(i,j-1))-r*vp(i,j)*dt
			endif
		enddo
	enddo

else
	
! #############################       LEAP FROG     #########################################################################
! ###########################################################################################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #########################################################################################
!if(leap==1)then	
	do i=1,xt-1
		do j=1,yt-1
			nf(i,j)=ni(i,j)-2*c4*(up(i+1,j)-up(i,j))-2*c5*(vp(i,j+1)-vp(i,j))
		enddo
	enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     #####################################################################################
	do i=1,xt
		do j=1,yt
				fb=(j-0.5*yt)*beta*dy*(0.00000000001)
!! CONDICAO DE CONTORNO
			if(j==1 .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)+2*c1*fb*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-2*c2*(np(i,j)-np(i-1,j))-r*up(i,j)*2*dt
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)+2*c1*fb*(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))-2*c2*(np(i,1)-np(i-1,1))-r*up(i,j)*2*dt

			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-2*c1*fb*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-2*c3*(np(i,j)-np(i,j-1))-r*vp(i,j)*2*dt
				uf(i,j)=ui(i,j)+2*c1*fb*(vp(i,j)+vp(i,j+1)+vp(xt,j+1)+vp(xt,j))-2*c2*(np(i,j)-np(xt-1,j))-r*up(i,j)*2*dt
			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-2*c1*fb*(up(i,j)+up(i,j-1)+up(1,j-1)+up(1,j))-2*c3*(np(1,j)-np(1,j-1))-r*vp(i,j)*2*dt
				uf(i,j)=ui(i,j)+2*c1*fb*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-2*c2*(np(1,j)-np(i-1,j))-r*up(i,j)*2*dt

			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			else	

				uf(i,j)=ui(i,j)+2*c1*fb*(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))-2*c2*(np(i,j)-np(i-1,j))-r*up(i,j)*2*dt
				vf(i,j)=vi(i,j)-2*c1*fb*(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))-2*c3*(np(i,j)-np(i,j-1))-r*vp(i,j)*2*dt
			endif
		enddo
	enddo

endif
	do i=1,xt-1
		do j=1,yt-1
			uout(i,j)=(uf(i,j)+uf(i+1,j))/2
			vout(i,j)=(vf(i,j)+vf(i,j+1))/2
		enddo
	enddo
	write(20)nf
	write(20)uout
	write(20)vout
	ni=np
	ui=up
	vi=vp
	np=nf
	up=uf
	vp=vf
write(*,*)'TEMPO =',n
enddo
end
