program shallow_water_v3
implicit none
integer,parameter::xt=240,yt=106,nt=4000
real,dimension(1:xt,1:yt)::ui,up,uf,vi,vp,vf
real,dimension(1:xt-1,1:yt-1)::ni,np,nf,uout,vout
real,dimension(1:yt)::lat,M,dx,f_cor
real,dimension(1:xt)::lon
real::dy,Nc,vmed,umed,umed0,vmed0,dhdx,dhdy,dudx,dudy,dvdx,dvdy,curvx,curvy,disx,disy,dnudx,dnvdy,dudx2,dudy2,dvdx2,dvdy2
real,parameter::H=200,dlon=1.5,dlat=1.5,dt=200,rad=3.1416/180,g=9.81,lon0=300,lat0=-10,r_earth=6371000,r=0.0000001,nic=0.0015
integer::i,j,n

! CFL udt/dx<1, dt<dx/c
!!!!!!!!!!!!!!!!!!!!!!!!!!! ################################## PRE PROCESSAMENTO
!!!!! CALCULO DE PARAMETROS 
do i=1,xt
	lon(i)=i*dlon
	!write(*,*)lon(i)
enddo
! dx e dy em coordenadas esfericas	
do j=1,yt	
	lat(j)=(j-0.5*yt)*dlat
	dx(j)=dlon*r_earth*cos(lat(j)*rad)*rad
	M(j)=1/(r_earth*cos(lat(j)*rad))
	f_cor(j)=2*((2*3.1416)/(3600*24))*sin(lat(j)*rad)
	write(*,*)lat(j)
	
enddo

vp=0
up=0
Nc=1/r_earth
dy=r_earth*dlat*rad
!write(*,*)1/M
!write(*,*)dy


!!!!! CALCULO DA ICBC
do i=1,xt-1
	do j=1,yt-1
		!if(lat(j)<0)np(i,j)=2*exp(-1*((lon(i)-lon0)/10)**2-1*((lat(j)-lat0)/10)**2)
		!if(lat(j)>0)np(i,j)=-2*exp(-1*((lon(i)-lon0)/10)**2-1*((lat(j)+lat0)/10)**2)
		np(i,j)=2*exp(-1*((lon(i)-lon0)/10)**2-1*((lat(j)-lat0)/10)**2)
		!if(np(i,j)>0.1)then
		!	np(i,j)=1
		!else
		!	np(i,j)=0
		!endif
	enddo
enddo
write(10)np
!!!!!!!!!! #################################################################### CALCULO DAS EQUACOES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     forward (n<=4)
do n=1,nt
	if(n<=4)then
!!!!!!!!!!! 	conservacao do momento
	do i=1,xt
		do j=1,yt
			disx=dt*r*up(i,j)
			disy=dt*r*vp(i,j)
!! condicoes de contorno

			if(j==1 .and. i/=1 .and. i/=xt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				dhdx=g*dt*(np(i,j)-np(i-1,j))/dx(j)
				curvx=dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*dt*(vp(i,j)-vp(i,yt))/dy
				dudy=vp(i,j)*dt*(up(i,j)-up(i,yt))/dy
				dvdx=up(i,j)*dt*(vp(i,j)-vp(i-1,j))/dx(j)
				dudx2=nic*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,yt))

				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=0
	
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vmed=(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))/4
				dhdx=g*dt*(np(i,1)-np(i-1,1))/dx(1)
				curvx=dt*Nc*up(i,j)*vmed*tan(((lat(1)+lat(j))/2)*rad)
				curvy=dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*dt*(vp(i,j)-vp(i-1,j))/dx(j)
				dudx2=nic*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*(dt/(dy**2))*(up(i,1)-2*up(i,j)+up(i,j-1))

				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=0

			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(xt,j+1)+vp(xt,j))/4
				umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				dhdx=dt*g*(np(i,j)-np(xt-1,j))/dx(yt-1)
				dhdy=g*dt*(np(i,j)-np(i,j-1))/dy
				curvx=dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*dt*(up(i,j)-up(xt,j))/dx(j)
				dvdy=vp(i,j)*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*dt*(vp(i,j)-vp(xt,j))/dx(j)

				dudx2=nic*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(xt,j))
				dudy2=nic*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*(dt/(dx(j)**2))*(vp(i+1,j)-2*vp(i,j)+vp(xt,j))
				dvdy2=nic*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))

				uf(i,j)=up(i,j)-dhdy+dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-curvy-disy-dvdy-dvdx+dudx2+dudy2


			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(1,j-1)+up(1,j)+up(i,j)+up(i,j-1))/4
				dhdx=dt*g*(np(1,j)-np(i-1,j))/dx(1)
				dhdy=g*dt*(np(1,j)-np(1,j-1))/dy
				curvx=dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*dt*(vp(i,j)-vp(i-1,j))/dx(j)

				dudx2=nic*(dt/(dx(j)**2))*(up(1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*(dt/(dx(j)**2))*(vp(1,j)-2*vp(i,j)+vp(i-1,j))
				dvdy2=nic*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))

				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-curvy-disy-dvdy-dvdx+dudx2+dudy2

			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0
!!!!!!!!!!!!!!
			else

				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				dhdx=g*dt*(np(i,j)-np(i-1,j))/dx(j)
				dhdy=g*dt*(np(i,j)-np(i,j-1))/dy
				curvx=dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*dt*(vp(i,j)-vp(i-1,j))/dx(j)

				dudx2=nic*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*(dt/(dx(j)**2))*(vp(i+1,j)-2*vp(i,j)+vp(i-1,j))
				dvdy2=nic*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))

				uf(i,j)=up(i,j)-dhdx+dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vp(i,j)-dhdy-dt*f_cor(j)*umed-curvy-disy-dvdy-dvdx+dvdx2+dvdy2
			endif
		enddo		
	enddo



!!!!!!!!!!!      conservacao de massa
	do i=1,xt-1
		do j=1,yt-1
			dudx=(up(i+1,j)-up(i,j))/dx(j)
			dvdy=(vp(i,j+1)*cos(lat(j+1)*rad)-vp(i,j)*cos(lat(j)*rad))/(dlat*rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   DNUDX
			if(i==xt-1)then
			dnudx=(dt/dx(j))*((np(1,j)+np(i,j))*up(1,j)/2-(np(i,j)+np(i-1,j))*up(i,j)/2)

			elseif(i==1)then
			dnudx=(dt/dx(j))*((np(i+1,j)+np(i,j))*up(i+1,j)/2-(np(i,j)+np(xt-1,j))*up(i,j)/2)

			else
			dnudx=(dt/dx(j))*((np(i+1,j)+np(i,j))*up(i+1,j)/2-(np(i,j)+np(i-1,j))*up(i,j)/2)

			endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   DNVDY
			if(j==yt-1)then
			dnvdy=(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			elseif(j==1)then
			dnvdy=(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			else
			dnvdy=(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			endif
			nf(i,j)=np(i,j)-H*dt*(dudx+M(j)*dvdy)-dnudx-dnvdy
		enddo
	enddo

	else
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!##################################################################### LEAP FROG
!!!!!!!!!!! 	conservacao do momento LEAP FROG
	do i=1,xt
		do j=1,yt
			disx=2*dt*r*up(i,j)
			disy=2*dt*r*vp(i,j)
!! teste de estabilidade
! CFL udt/dx<1, dt<dx/u
			if(up(i,j)*dt/dx(j)>=1)then
				write(*,*)"SIMULACAO INSTAVEL, TEMPO =",n
				!write(*,*)up(i,j)*dt/dx
				GO TO 10
			
			endif
!! condicoes de contorno               LEAP FROG
			if(j==1 .and. i/=1 .and. i/=xt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				dhdx=g*2*dt*(np(i,j)-np(i-1,j))/dx(j)
				curvx=2*dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				dudx=up(i,j)*2*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*2*dt*(vp(i,j)-vp(i,yt))/dy
				dudy=vp(i,j)*2*dt*(up(i,j)-up(i,yt))/dy
				dvdx=up(i,j)*2*dt*(vp(i,j)-vp(i-1,j))/dx(j)
				dudx2=nic*2*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*2*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,yt))


				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=0
	
			elseif(j==yt .and. i/=1 .and. i/=xt)then
				vmed=(vp(i,j)+vp(i,1)+vp(i-1,1)+vp(i-1,j))/4
				dhdx=g*2*dt*(np(i,1)-np(i-1,1))/dx(1)
				curvx=2*dt*Nc*up(i,j)*vmed*tan(((lat(1)+lat(j))/2)*rad)
				dudx=up(i,j)*2*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*2*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*2*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*2*dt*(vp(i,j)-vp(i-1,j))/dx(j)
				dudx2=nic*2*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*2*(dt/(dy**2))*(up(i,1)-2*up(i,j)+up(i,j-1))


				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=0


			elseif(i==1 .and. j/=1 .and. j/=yt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(xt,j+1)+vp(xt,j))/4
				umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				dhdx=2*dt*g*(np(i,j)-np(xt-1,j))/dx(yt-1)
				dhdy=g*2*dt*(np(i,j)-np(i,j-1))/dy
				curvx=2*dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=2*dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*2*dt*(up(i,j)-up(xt,j))/dx(j)
				dvdy=vp(i,j)*2*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*2*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*2*dt*(vp(i,j)-vp(xt,j))/dx(j)
				dudx2=nic*2*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(xt,j))
				dudy2=nic*2*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*2*(dt/(dx(j)**2))*(vp(i+1,j)-2*vp(i,j)+vp(xt,j))
				dvdy2=nic*2*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))

				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed-curvy-disy-dvdy-dvdx+dvdx2+dvdy2



			elseif(i==xt .and. j/=1 .and. j/=yt)then
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(1,j-1)+up(1,j)+up(i,j)+up(i,j-1))/4
				dhdx=2*dt*g*(np(1,j)-np(i-1,j))/dx(1)
				dhdy=g*2*dt*(np(1,j)-np(1,j-1))/dy
				curvx=2*dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=2*dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*2*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*2*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*2*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*2*dt*(vp(i,j)-vp(i-1,j))/dx(j)

				dudx2=nic*2*(dt/(dx(j)**2))*(up(1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*2*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*2*(dt/(dx(j)**2))*(vp(1,j)-2*vp(i,j)+vp(i-1,j))
				dvdy2=nic*2*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))


				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed-curvy-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed+curvx-disy-dvdy-dvdx+dvdx2+dvdy2



			elseif(i==1 .or. i==xt .and. j==1)then
				vf(i,j)=0
				uf(i,j)=0
			elseif(i==1 .or. i==xt .and. j==yt)then
				vf(i,j)=0
				uf(i,j)=0
!!!!!!!!!!!!!!
			else
				vmed=(vp(i,j)+vp(i,j+1)+vp(i-1,j+1)+vp(i-1,j))/4
				umed=(up(i+1,j-1)+up(i+1,j)+up(i,j)+up(i,j-1))/4
				dhdx=g*2*dt*(np(i,j)-np(i-1,j))/dx(j)
				dhdy=g*2*dt*(np(i,j)-np(i,j-1))/dy
				curvx=2*dt*Nc*up(i,j)*vmed*tan(((lat(j+1)+lat(j))/2)*rad)
				curvy=2*dt*Nc*(umed**2)*tan(lat(j)*rad)
				dudx=up(i,j)*2*dt*(up(i,j)-up(i-1,j))/dx(j)
				dvdy=vp(i,j)*2*dt*(vp(i,j)-vp(i,j-1))/dy
				dudy=vp(i,j)*2*dt*(up(i,j)-up(i,j-1))/dy
				dvdx=up(i,j)*2*dt*(vp(i,j)-vp(i-1,j))/dx(j)

				dudx2=nic*2*(dt/(dx(j)**2))*(up(i+1,j)-2*up(i,j)+up(i-1,j))
				dudy2=nic*2*(dt/(dy**2))*(up(i,j+1)-2*up(i,j)+up(i,j-1))

				dvdx2=nic*2*(dt/(dx(j)**2))*(vp(i+1,j)-2*vp(i,j)+vp(i-1,j))
				dvdy2=nic*2*(dt/(dy**2))*(vp(i,j+1)-2*vp(i,j)+vp(i,j-1))

				uf(i,j)=ui(i,j)-dhdx+2*dt*f_cor(j)*vmed+curvx-disx-dudx-dudy+dudx2+dudy2
				vf(i,j)=vi(i,j)-dhdy-2*dt*f_cor(j)*umed-curvy-disy-dvdy-dvdx+dvdx2+dvdy2
			endif
		enddo		
	enddo



!!!!!!!!!!!      conservacao de massa        LEAP FROG
	do i=1,xt-1
		do j=1,yt-1
			dudx=(up(i+1,j)-up(i,j))/dx(j)
			dvdy=(vp(i,j+1)*cos(lat(j+1)*rad)-vp(i,j)*cos(lat(j)*rad))/(dlat*rad)
			if(i==xt-1)then
			dnudx=2*(dt/dx(j))*((np(1,j)+np(i,j))*up(1,j)/2-(np(i,j)+np(i-1,j))*up(i,j)/2)

			elseif(i==1)then
			dnudx=2*(dt/dx(j))*((np(i+1,j)+np(i,j))*up(i+1,j)/2-(np(i,j)+np(xt-1,j))*up(i,j)/2)

			else
			dnudx=2*(dt/dx(j))*((np(i+1,j)+np(i,j))*up(i+1,j)/2-(np(i,j)+np(i-1,j))*up(i,j)/2)

			endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   DNVDY
			if(j==yt-1)then
			dnvdy=2*(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			elseif(j==1)then
			dnvdy=2*(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			else
			dnvdy=2*(dt/dy)*((np(i,j+1)+np(i,j))*vp(i,j+1)/2-(np(i,j)+np(i,j-1))*vp(i,j)/2)

			endif
			nf(i,j)=ni(i,j)-H*dt*2*(dudx+M(j)*dvdy)-dnudx-dnvdy
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
