program shallow_water_v1
implicit none
real,dimension(1:240,1:91)::ui,up,uf,vi,vp,vf
real,dimension(1:239,1:90)::ni,np,nf,uout,vout
real,dimension(1:91)::lat,M,dx
real,dimension(1:240)::lon
real::dy,Nc,vmed,umed,umed0,vmed0,dhdx,dhdy
real,parameter::H=200,dlon=1.5,dlat=1.5,dt=600,rad=3.1416/180,g=9.81,lon0=300,lat0=-10,r_earth=6371000
integer::i,j,n
integer,parameter::xt=240,yt=91,nt=2000
! CFL udt/dx<1, dt<dx/c
!!!!!!!!!!!!!!!!!!!!!!!!!!! ################################## PRE PROCESSAMENTO
!!!!! CALCULO DE PARAMETROS 
do i=1,xt
	lon(i)=i*dlon
	!write(*,*)lon(i)
enddo
! dx e dy em coordenadas esfericas	
do j=1,yt	
	lat(j)=(j-0.5*(yt+1))*dlat
	dx(j)=dlon*r_earth*cos(lat(j)*rad)*rad
	M(j)=1/(r_earth*cos(lat(j)*rad))

	
enddo
up=0
vp=0
dy=r_earth*dlat*rad
!write(*,*)1/M
!write(*,*)dy


!!!!! CALCULO DA ICBC
do i=1,xt-1
	do j=1,yt-1
		np(i,j)=exp(-1*((lon(i)-lon0)/10)**2-1*((lat(j)-lat0)/10)**2)
	enddo
enddo
!write(10)np
!!!!!!!!!! #################################################################### CALCULO DAS EQUACOES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     forward (n<=4)
do n=1,nt
	if(n<=4)then
!!!!!!!!!!! 	conservacao do momento
	do i=1,xt
		do j=1,yt
!! condicoes de contorno

			if(j==1 .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=up(i,j)-g*dt*(np(i,j)-np(i-1,j))/dx(j)
	
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=up(i,j)-g*dt*(np(i,1)-np(i-1,1))/dx(1)

			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vf(i,j)=vp(i,j)-g*dt*(np(i,j)-np(i,j-1))/dy
				uf(i,j)=up(i,j)-dt*g*(np(i,j)-np(xt-1,j))/dx(yt-1)
			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vf(i,j)=vp(i,j)-g*dt*(np(1,j)-np(1,j-1))/dy
				uf(i,j)=up(i,j)-dt*g*(np(1,j)-np(i-1,j))/dx(1)
			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0
!!!!!!!!!!!!!!
			else

				!vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				!umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				uf(i,j)=up(i,j)-g*dt*(np(i,j)-np(i-1,j))/dx(j)
				vf(i,j)=vp(i,j)-g*dt*(np(i,j)-np(i,j-1))/dy
			endif
		enddo		
	enddo



!!!!!!!!!!!      conservacao de massa
	do i=1,xt-1
		do j=1,yt-1
			dhdx=(up(i+1,j)-up(i,j))/dx(j)
			dhdy=(vp(i,j+1)*cos(lat(j+1)*rad)-vp(i,j)*cos(lat(j)*rad))/(dlat*rad)
			nf(i,j)=np(i,j)-H*dt*(dhdx+M(j)*dhdy)
		enddo
	enddo

	else
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!! 	conservacao do momento LEAP FROG
	do i=1,xt
		do j=1,yt
!! teste de estabilidade
! CFL udt/dx<1, dt<dx/u
			if(up(i,j)*dt/dx(j)>=1)then
				write(*,*)"SIMULACAO INSTAVEL, TEMPO =",n
				!write(*,*)up(i,j)*dt/dx
				GO TO 10
			
			endif
!! condicoes de contorno               LEAP FROG
			if(j==1 .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)-g*2*dt*(np(i,j)-np(i-1,j))/dx(j)
	
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vf(i,j)=0
				uf(i,j)=ui(i,j)-g*2*dt*(np(i,1)-np(i-1,1))/dx(1)

			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-g*2*dt*(np(i,j)-np(i,j-1))/dy
				uf(i,j)=ui(i,j)-2*dt*g*(np(i,j)-np(xt-1,j))/dx(yt-1)
			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vf(i,j)=vi(i,j)-g*2*dt*(np(1,j)-np(1,j-1))/dy
				uf(i,j)=ui(i,j)-2*dt*g*(np(1,j)-np(i-1,j))/dx(1)
			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0
!!!!!!!!!!!!!!
			else
				!vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				!umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				uf(i,j)=ui(i,j)-g*2*dt*(np(i,j)-np(i-1,j))/dx(j)
				vf(i,j)=vi(i,j)-g*2*dt*(np(i,j)-np(i,j-1))/dy
			endif
		enddo		
	enddo



!!!!!!!!!!!      conservacao de massa        LEAP FROG
	do i=1,xt-1
		do j=1,yt-1
			dhdx=(up(i+1,j)-up(i,j))/dx(j)
			dhdy=(vp(i,j+1)*cos(lat(j+1)*rad)-vp(i,j)*cos(lat(j)*rad))/(dlat*rad)
			nf(i,j)=ni(i,j)-H*dt*2*(dhdx+M(j)*dhdy)
		enddo
	enddo
endif
!!!!!!!!!! #################################################################### POS - PROCESSAMENTO
	ni=np
	vi=vp
	ui=up
	do i=1,xt-1
		do j=1,yt-1
			uout(i,j)=(uf(i,j)+uf(i+1,j))/2
			vout(i,j)=(vf(i,j)+vf(i,j+1))/2
		enddo
	enddo
	write(20)nf
	write(20)uout
	write(20)vout
	np=nf
	vp=vf
	up=uf
enddo
10 CONTINUE
end
