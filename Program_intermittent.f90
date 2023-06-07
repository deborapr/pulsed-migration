!Model two islands sea level oscilation 
program lake
	implicit none 
	integer n, num, nger, m, idum, i, j, nfilhos, NG, nrep, contnrep,t_prox
	integer medesp,ngrupo, ngrupo1,ngrupo2
	real*8 c,q, tmut, conect, c_prox  
	parameter (q=0.05d0,m=2000, num=400) 
	integer contl1,contl2,ne,ns
	integer pai, mae, dif, nachou, ntent, intq,m1, m2, nl,nmeio, mesp, acada, l, ll
	integer  contmae, contpai, ncontesp,ncontesp1,ncontesp2, npais1,npais2,npais12 , pai0,novo , lagoa, cont !SA
	real*8   mut,rand
	Integer, ALLOCATABLE :: P(:,:), Pn(:,:),Pmeio(:,:),auxiliar(:)
	intq=nint(m*q) 
	nl=m+2

OPEN(UNIT=62,FILE='input_sea.txt',STATUS='old')
    read(62,*)  nger,c,nrep,tmut, acada
OPEN(UNIT=60,FILE='Temporal.txt',STATUS='UNKNOWN')
	write(60,*) 'pop ', 'rep ','mig ','tem ', 'esp ', 'abund1 ',  'abund2 '
	close(60)
	
		medesp=0
		ALLOCATE( P(num, m+2), Pn(num, m+2), Pmeio(num/2,m+2),auxiliar(m+2)) 
			do contnrep=1,nrep
				close(62)
				OPEN(UNIT=62,FILE='input_sea.txt',STATUS='old')
				read(62,*) nger, c 
				P=0
				Pn=0
			    ! Initial distribution	
				P(1:num/2,m+2)=1 
				P(num/2+1:num,m+2)=2 
				P(:,m+1)=1
				mesp=1
				ns=1
				read(62,*)  t_prox,c_prox
				!Dynamics
				do i=1, nger
					if(i.eq.t_prox.and.t_prox.lt.nger) then
						c=c_prox
						read(62,*)t_prox,c_prox
					endif
				! REPRODUCTION
					Pn=0
					contmae=0
					nfilhos=0 
					npais1=count(P(:,m+2).eq.1)
					npais2=count(P(:,m+2).eq.2)
					
					do while (nfilhos<num)
						if (nfilhos<num/2) then
							npais12=npais1
							pai0=1	
						else 
							npais12=npais2
							pai0=npais1+1
						endif
						call random_number(rand)
						mae=pai0+Int(rand*(npais12)) 
						nachou=0
						ntent=npais12
						do contpai=1, ntent
								call random_number(rand)	
								pai=pai0+Int(rand*(npais12)) 
								if (mae.ne.pai) then
									dif=0
									do j=1, m  
										dif=dif+abs(P(mae,j)-P(pai,j))
									enddo
									if (dif<=intq) then
										nachou=1 
										nfilhos=nfilhos+1
										do j=1, m 
											call random_number(rand)
											if(rand<=0.5) then
												Pn(nfilhos,j)=P(mae,j)
												else
												Pn(nfilhos,j)=P(pai,j)
											endif
											call random_number(mut)
											if (mut<=tmut) Pn(nfilhos,j)= abs(Pn(nfilhos,j)-1)
										enddo
										Pn(nfilhos,m+2)= P(mae, m+2)
										Pn(nfilhos,m+1)= P(mae,m+1)
									endif
								endif
								if (nachou==1) exit
						enddo
					enddo
					P=Pn
					! Identifying species
					if(mod(i,acada).eq.0) then
						call species (P,  num,nl,intq,ncontesp,mesp)
	    				OPEN(UNIT=60,FILE='Temporal.txt',STATUS='old',Access = 'append')
						Do ne=1, mesp
						contl1=0
						contl2=0
							Do j=1, num
								if(P(j,m+2).eq.1.and.P(j,m+1).eq.ne) then
									contl1=contl1+1
								endif
								if(P(j,m+2).eq.2.and.P(j,m+1).eq.ne) then
									contl2=contl2+1
								endif
!
							Enddo
							if(contl1.ne.0.or.contl2.ne.0)then
							write(60,*) num, contnrep,c,i, ne,contl1,contl2
							endif
						Enddo	
						close(60)
					endif
					!MIGRATION
					do j=1, num
						call random_number(conect)
						novo=0
						if (conect<c) then
							if (P(j,m+2).eq.1) then
								novo=2
							else
								novo=1
							endif	
							P(j,m+2)=novo
						endif
					enddo
				  ! SORTING
					Pn=0
					cont=1
					Do lagoa=1, 2
						Do j=1, num
							if (P(j,m+2).eq.lagoa)then
							Pn(cont,:)=P(j,:)
							cont=cont+1
							endif
						enddo
					enddo	
					P=Pn
				enddo
			enddo
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!))!!!!!
subroutine species(MaV,np,nl,nG,ncontesp,mesp) !P,Nv,num,nl,intq,ncontesp
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension::MaV(np,nl)
	Integer, ALLOCATABLE :: NES(:),Nes2(:),Nfreq(:,:),Nlo(:),Nv(:)
	Allocate(Nlo(np), Nv(np), NES(mesp))
	nfim=0
	nGes=NG  
	ngrupo=0
	Nv=0
	Nlo=0
!!!!!
	!!
	Do while(nfim.eq.0) 
		nfim=1
		Do l=1,np 
			if(Nv(l).eq.0)then 
				nfim=0
				ngrupo=ngrupo+1
				Nv(l)=ngrupo
				go to 22	
			endif
		Enddo
22      	continue
		nbus=0
		Do while (nbus.eq.0) 
			nbus=1
			Do i=1,np
				  if(Nv(i).eq.ngrupo.and.Nlo(i).eq.0) then	
					Nlo(i)=1 
					Do ii=1,np
						if(Nv(ii).eq.0) then 
							ndiff=0
							Do j=1,nl-2 !genoma
							 n1=(MaV(i,j))
							 n2=(MaV(ii,j))
							 ndiff=ndiff+abs(n1-n2)		
							Enddo
							if(ndiff.le.NG)then 
								nbus=0							
								Nv(ii)=ngrupo
								else
								if(ndiff.le.NGes) then
									if(MaV(i,nl-1).eq.MaV(ii,nl-1)) then
										nbus=0							
										Nv(ii)=ngrupo
									endif
								endif
							endif
						endif
					enddo
				endif
			Enddo
		Enddo
	Enddo
	
	!!!!!!!!!!
	ALLOCATE(Nfreq(ngrupo,2))
	Nes=0
	nmax=mesp
	ncontesp=0 
	Do n=1,nmax
		Nfreq=0
		Do i=1,np
			mva=(MaV(i,nl-1))
			if(mva.eq.n)then
				Nfreq(nv(i),2)=Nfreq(nv(i),2)+1	
				Nfreq(nv(i),1)=nv(i)	
			endif
		Enddo
		call sortc(Nfreq,ngrupo)
		ngesco=Nfreq(1,1)
		if(Nes(n).eq.0)then
			Do i=1,np
				IF(Nv(i).eq.ngesco)then
					MaV(i,nl-1)=n
					Nes(n)=1
				endif				
			Enddo
		else
		print*, 'Attention!passou no else'
			mesp=mesp+1
			Do i=1,np
				if(Nv(i).eq.Nfreq(1,1))then
					MaV(i,nl-1)=mesp
				endif
			enddo
		endif
		if(Nes(n).gt.0)	ncontesp=ncontesp+1 
		if(ngrupo.gt.1) then
			do nn=2,ngrupo
					if(Nfreq(nn,2).ge.1)then
						mesp=mesp+1
						ncontesp=ncontesp+1 !SA
						Do i=1,np
							if(Nv(i).eq.Nfreq(nn,1))then
								MaV(i,nl-1)=mesp
							endif
						enddo
					Endif
			Enddo
		endif
	enddo

      endsubroutine species
!!!!!!!!!!!!!!!!!
	subroutine sortc(Nio,nb)
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension:: Nio(nb,2)
	do i=1,nb-1
		ip=nb+1-i
		ipa=nb+1-i-1
		if(Nio(ip,2).gt.Nio(ipa,2)) then
			naux=Nio(ip,2)
			Nio(ip,2)=Nio(ipa,2)
			Nio(ipa,2)=naux
			naux=Nio(ip,1)
			Nio(ip,1)=Nio(ipa,1)
			Nio(ipa,1)=naux
		endif
	enddo	
	endsubroutine sortc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

