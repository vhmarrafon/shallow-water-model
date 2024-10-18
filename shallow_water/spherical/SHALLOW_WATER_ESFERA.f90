program VERSION_2
implicit none
real,dimension(1:240,1:91)::up,uf,vp,vf,ui,vi,lon,lat
real,dimension(1:239,1:90)::np,nf,ni,uout,vout
real,dimension(1:91)::M,f_cor
real::Nc,f,dhdx,dhdy,npdw,umed,vmed
integer::i,n,j
real,parameter::dt=1200,H=200,pi=3.1416,r=0.0000001,g=9.81,beta=2.42,rad=3.1416/180,r_earth=6371000
real,parameter::dx=1.5,dy=1.5
real,parameter::lat0=-10,lon0=300
integer,parameter::xt=240,yt=91,nt=2000
!coriolis f=0.00010284
!!!!! NUNCA ADICIONAR TEMPO NA MATRIZ
! CFC cdt/dx<1, dt<dx/c
!write(*,*)1/dx !n de pontos em x = 100
!write(*,*)0.81/dt !n de pontos em t = 90
Nc=1/r_earth
up=0
vp=0
ui=0
vi=0

!!ORDEM do DO, tempo,x,y
do i=1,xt
	do j=1,yt

		lon(i,j)=(i-1)*dx
		lat(i,j)=(j-(yt+1)/2.)*dy
		f_cor(j)=2*(2*3.1416/(3600*24))*sin(lat(i,j)*rad)
		M(j)=1/(r_earth*cos(lat(i,j)*rad))
		if(i==1 .or. i==xt)write(*,*)lat(i,j)
	enddo
enddo

do i=1,xt-1
	do j=1,yt-1

		np(i,j)=exp((-1.*((lon(i,j)-lon0)/10.)**2)-1.*(((lat(i,j)-lat0)/10.)**2))
		!write(*,*)np(i,j)
	enddo
enddo

do n=1,nt
if(n<=4)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #########################################################################################
	do i=1,xt-1
		do j=1,yt-1
			npdw=(up(i+1,j)-up(i,j))/dx+(vp(i,j+1)*cos(lat(i,j+1)*rad)-vp(i,j)*cos((lat(i,j)+lat(i,j+1))*rad/2))/dy
			nf(i,j)=np(i,j)-H*M(j)*dt*npdw
		enddo
	enddo
	do i=1,xt
		do j=1,yt
			if(i>=xt-1 .and. j<yt-1 .and. j/=1)then
				umed=(up(i,j)+up(i,j-1)+up(1,j-1)+up(1,j))/4
				dhdy=Nc*g*dt*(np(1,j+1)-np(1,j))/dy
				uf(i,j)=0
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-dt*Nc*(umed**2)*tan(lat(i,j)*rad)


			elseif(i==1 .and. j<yt-1 .and. j/=1)then
				umed=(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))/4
				dhdy=Nc*g*dt*(np(i+1,j+1)-np(i,j))/dy
				uf(i,j)=0
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-dt*Nc*(umed**2)*tan(lat(i,j)*rad)


			elseif(j>=yt-1 .and. i<xt-1 .and. i/=1)then
				vmed=(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))/4
				dhdx=M(j)*g*dt*(np(i+1,1)-np(i,1))/dx
				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+Nc*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,1))*rad/2)
				vf(i,j)=0


			elseif(j==1 .and. i<xt-1 .and. i/=1)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				dhdx=M(j)*g*dt*(np(i+1,j)-np(i,j))/dx
				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+Nc*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,j+1))*rad/2)
				vf(i,j)=0

			elseif(j==1 .or. j>=yt-1 .and. i>=xt-1)then
				uf(i,j)=0
				vf(i,j)=0
			elseif(j==1 .or. j>=yt-1 .and. i==1)then
				uf(i,j)=0
				vf(i,j)=0
			else
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))/4
				dhdx=M(j)*g*dt*(np(i+1,j)-np(i,j))/dx
				dhdy=Nc*g*dt*(np(i+1,j+1)-np(i,j))/dy
				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+Nc*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,j+1))*rad/2)
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-dt*Nc*(umed**2)*tan(lat(i,j)*rad)
			endif
		enddo		
	enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  ###########################################################################################
else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  LEAP FROG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  LEAP FROG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  LEAP FROG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  LEAP FROG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,xt-1
		do j=1,yt-1
			!! variacao do vento utilizado para calculo de nf
			npdw=(up(i+1,j)-up(i,j))/dx+(vp(i,j+1)*cos(lat(i,j+1)*rad)-vp(i,j)*cos((lat(i,j)+lat(i,j+1))*rad/2))/dy
			nf(i,j)=ni(i,j)-H*M(j)*2*dt*npdw
		enddo
	enddo
	do i=1,xt
		do j=1,yt

			if(i>=xt-1 .and. j/=yt .and. j/=1)then
				umed=(up(i,j)+up(i,j-1)+up(1,j-1)+up(1,j))/4
				dhdy=Nc*g*dt*2*(np(1,j+1)-np(1,j))/dy
				uf(i,j)=0
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed-2*dt*Nc*(umed**2)*tan(lat(i,j)*rad)


			elseif(i==1 .and. j/=yt .and. j/=1)then
				umed=(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))/4
				dhdy=Nc*g*dt*2*(np(i,j+1)-np(i,j))/dy
				uf(i,j)=0
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed-2*dt*Nc*(umed**2)*tan(lat(i,j)*rad)


			elseif(j>=yt-1 .and. i/=xt .and. i/=1)then
				vmed=(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))/4
				dhdx=M(j)*g*2*dt*(np(i+1,1)-np(i,1))/dx
				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+Nc*2*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,1))*rad/2)
				vf(i,j)=0

			elseif(j==1 .and. i/=xt .and. i/=1)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				dhdx=M(j)*g*2*dt*(np(i+1,j)-np(i,j))/dx
				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+Nc*2*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,j+1))*rad/2)
				vf(i,j)=0

			elseif(j==1 .or. j==yt .and. i==xt)then
				uf(i,j)=0
				vf(i,j)=0
			elseif(j==1 .or. j==yt .and. i==1)then
				uf(i,j)=0
				vf(i,j)=0
			else
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(i,j)+up(i,j-1)+up(i+1,j-1)+up(i+1,j))/4
				dhdx=M(j)*g*2*dt*(np(i+1,j)-np(i,j))/dx
				dhdy=Nc*g*dt*2*(np(i,j+1)-np(i,j))/dy
				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+Nc*2*dt*up(i,j)*vmed*tan((lat(i,j)+lat(i,j+1))*rad/2)
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed-2*dt*Nc*(umed**2)*tan(lat(i,j)*rad)
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
ui=up
vi=vp
ni=np
np=nf
vp=vf
up=uf
enddo

!write(20)np
end
