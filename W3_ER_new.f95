MODULE network
IMPLICIT NONE
TYPE node_list
	INTEGER :: names		  !index of the node i
	INTEGER,DIMENSION(2) :: s			  !spin value
	REAL*16 :: h                      !field value
	REAL*16,DIMENSION(2,3)::p	  !probability at time t,t-1,t-2 of P(si=+or-1)
	REAL*16,DIMENSION(1:2,1:2,2)::pp  !PAIR PROB P(si,sj) at time t and t-1 for the 4 
	INTEGER,DIMENSION(2)::s1	  ! spin +1 or -1
	INTEGER::degree			  !degree of the node i
END TYPE node_list

TYPE neighbour_list
	INTEGER,DIMENSION(2) :: sn	  !value of spin +1 or -1
	REAL*16 ::JJ			  !coupling constant
	INTEGER:: na			  !names of the neighbour (refers to the node names)
	INTEGER:: na2
	REAL*16,DIMENSION(2,2,2)::pp	  !PAIR PROB P(si,sj) at time t and t-1 for the 4 possibilities
	REAL*16,DIMENSION(2,2)::pp_t
END TYPE neighbour_list
CONTAINS
!***************************** ORDERING SUBROUTINE***************************************************
SUBROUTINE ord(x_datas)
	INTEGER                        :: n,i,j
        REAL                           ::x_temp
	INTEGER,DIMENSION(:),INTENT(INOUT)::x_datas
        n=SIZE(x_datas)
DO i=1,n-1
DO j=1,n-1
IF (x_datas(j)<x_datas(j+1)) THEN
x_temp=x_datas(j)
x_datas(j)=x_datas(j+1)
x_datas(j+1)=x_temp
END IF
END DO
END DO 
END SUBROUTINE ord
!*****************************RANDOM NUMBER GENERATOR************************************************
REAL*8 FUNCTION lcg(seed) RESULT(random)
	INTEGER*8 :: a,c,m,seed	
	INTEGER*8,SAVE::Oldseed=0 !per cambiare simulazione bisogna cambiare l'oldseed!!!!! tutto dipende da questo!!!! 
	a = 9301
	c = 49297
	m = 233280
	If (OldSeed .EQ. 0) OldSeed = Seed
	OldSeed = Mod(a*OldSeed+c, m)
	random = REAL(Oldseed)/m
END FUNCTION lcg
!*******************************************************************************************************
!*********************INITIALIZATION NODE  AND NEIGHBOUR LIST*******************************************
!*******************************************************************************************************
SUBROUTINE init_prob_ER(array,array2,matrix,a,b,pp)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix
INTEGER::i,j
REAL*16::x,xx,a,b,pp

DO i=1,SIZE(array)
	CALL random_number(x)
	xx=a+(b-a)*x
	array(i)%p(1,3)=xx; array(i)%p(2,3)=1-xx	!STAR initialization
	array2(i)%p(1,2)=xx; array2(i)%p(2,2)=1-xx	!diamond initialization
	IF(array(i)%degree .ne. 0)THEN	
	WHERE (matrix(:,:)%na == i )
		matrix(:,:)%pp(1,1,2)=xx*pp			! 		1
		matrix(:,:)%pp(2,1,2)=(1-xx)*pp			! 		2
		matrix(:,:)%pp(1,2,2)=xx*(1-pp)			! 		3
		matrix(:,:)%pp(2,2,2)=(1-xx)*(1-pp)		! 		4
	END WHERE
	END IF
	array2(i)%pp(1,1,2)=xx*pp			! 		1
	array2(i)%pp(2,1,2)=(1-xx)*pp			! 		2
	array2(i)%pp(1,2,2)=xx*(1-pp)			! 		3
	array2(i)%pp(2,2,2)=(1-xx)*(1-pp)		! 		4
END DO
DO i=1,SIZE(array2)
	IF(array(i)%degree.ne.0)THEN
		DO j=1,array2(i)%degree
			array2(i)%p(1,3)=matrix(i,j)%pp(1,1,2)+matrix(i,j)%pp(2,1,2)
			array2(i)%p(2,3)=matrix(i,j)%pp(1,2,2)+matrix(i,j)%pp(2,2,2)
		END DO
	ELSE
		array2(i)%p(1,3)=pp;array2(i)%p(2,3)=1-pp
	END IF
END DO
END SUBROUTINE init_prob_ER

SUBROUTINE init_prob_ER2(array,array2,matrix,a,b,a1,b1)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix
INTEGER::i,j
REAL*16::x,xx,y,yy,a,b,a1,b1,pp
DO i=1,SIZE(array)
	CALL random_number(y)
	CALL random_number(x)
	xx=a+(b-a)*x
	yy=a1+(b1-a1)*y

	array(i)%p(1,3)=yy; array(i)%p(2,3)=1-yy	
	array2(i)%p(1,3)=xx; array2(i)%p(2,3)=1-xx
	array2(i)%p(1,2)=yy; array2(i)%p(2,2)=1-yy
END DO

DO i=1,SIZE(array2)
	DO j=1,array2(i)%degree
		matrix(i,j)%pp(1,1,2)=array2(i)%p(1,3)*array2(matrix(i,j)%na)%p(1,2)
		matrix(i,j)%pp(2,1,2)=array2(i)%p(1,3)*(1-array2(matrix(i,j)%na)%p(1,2))
		matrix(i,j)%pp(1,2,2)=(1-array2(i)%p(1,3))*array2(matrix(i,j)%na)%p(1,2)
		matrix(i,j)%pp(2,2,2)=(1-array2(i)%p(1,3))*(1-array2(matrix(i,j)%na)%p(1,2))
	END DO
	array2(i)%pp(1,1,2)=array2(i)%p(1,3)*array2(i)%p(1,2)
	array2(i)%pp(2,1,2)=array2(i)%p(1,3)*(1-array2(i)%p(1,2))
	array2(i)%pp(1,2,2)=(1-array2(i)%p(1,3))*array2(i)%p(1,2)
	array2(i)%pp(2,2,2)=(1-array2(i)%p(1,3))*(1-array2(i)%p(1,2))
END DO
END SUBROUTINE init_prob_ER2

SUBROUTINE init_prob_ER3(array,array2,matrix,a,b,pp1)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix
INTEGER::i,j
REAL*16::x,xx,a,b
REAL*16,DIMENSION(:)::pp1

DO i=1,SIZE(array)
	CALL random_number(x)
	xx=a+(b-a)*x
	array(i)%p(1,3)=xx; array(i)%p(2,3)=1-xx	!STAR initialization
	array2(i)%p(1,3)=pp1(1)+pp1(2);array2(i)%p(2,3)=pp1(3)+pp1(4)
	array2(i)%p(1,2)=pp1(1)+pp1(3);array2(i)%p(2,2)=pp1(2)+pp1(4)
	DO j=1,array2(i)%degree
		matrix(i,j)%pp(1,1,2)=pp1(1)
		matrix(i,j)%pp(2,1,2)=pp1(2)
		matrix(i,j)%pp(1,2,2)=pp1(3)
		matrix(i,j)%pp(2,2,2)=pp1(4)
		array2(i)%pp(1,1,2)=pp1(1)			! 		1
		array2(i)%pp(2,1,2)=pp1(2)			! 		2
		array2(i)%pp(1,2,2)=pp1(3)			! 		3
		array2(i)%pp(2,2,2)=pp1(4)			! 		4
	END DO
END DO
END SUBROUTINE init_prob_ER3



SUBROUTINE copy_prob(array,array2,matrix,matrix2)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix,matrix2
array2=array
matrix2=matrix
END SUBROUTINE copy_prob


SUBROUTINE init_voterER_MC(array)
TYPE(node_list),DIMENSION(:)::array
INTEGER::i
REAL*16::x
DO i=1,SIZE(array)
CALL random_number(x)
IF (x>array(i)%p(2,3))THEN
	array(i)%s=1
ELSE
	array(i)%s=-1
END IF
END DO
END SUBROUTINE init_voterER_MC
!*************************************NEIGHBOUR LIST****************************************************
SUBROUTINE Voter_ER(array,array2,matrix,jin,hin,pin,a,b,a1,b1,pp,pp1,link)
	TYPE(node_list),DIMENSION(:)::array,array2
	TYPE(neighbour_list),DIMENSION(:,:)::matrix
	TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::temp	
	REAL*16::a,b,a1,b1,pp
	REAL*16,DIMENSION(:)::pp1
	INTEGER::i,j,l,k,n,m,si,sj,link
	REAL*16::x,jin,pin,hin
	n=SIZE(array)
	m=SIZE(matrix(1,:))

DO i=1,n
array(i)%p=0;array2(i)%p=0;array(i)%names=0;array2(i)%names=0;array(i)%pp=0;array2(i)%pp=0
DO j=1,m
matrix(i,j)%pp=0;matrix(i,j)%pp_t=0;matrix(i,j)%na=0;matrix(i,j)%jj=0
END DO 
END DO


	DO i=1,n
		array(i)%names=i;array2(i)%names=i		!initialization index
		array(i)%h=hin;array2(i)%h=hin   		!initialization field  "-0.5+lcg(seed)"
		DO j=1,m
		matrix(i,j)%jj=jin				!initialization coupling constant	
		matrix(i,j)%na=j
		END DO
		matrix(i,i)%jj=0	!no self loop condition
		matrix(i,i)%na=0
	END DO
	DO i=1,n				
		DO j=1,m		!initialization 2-times prob,coupling&index following ER rules
		CALL RANDOM_NUMBER(x)
		IF(x<pin)THEN
		matrix(i,j)%pp=0
		matrix(j,i)%pp=0
		matrix(i,j)%pp_t=0
		matrix(j,i)%pp_t=0
		matrix(i,j)%jj=0
		matrix(j,i)%jj=0
		matrix(i,j)%na=0
		matrix(j,i)%na=0
		END IF
		END DO
	END DO

	DO i=1,n			!Put all the non-zero link in front of the neigh list
		CALL ord(matrix(i,:)%na)
	END DO
	l=1
	DO i=1,n			!compute array degree
		k=1
		IF(matrix(i,k)%na.eq.0)THEN
			array(i)%degree=0;array2(i)%degree=0
		END IF
		DO WHILE(matrix(i,k)%na.ne.0)
			array(i)%degree=k;array2(i)%degree=k
			matrix(i,k)%na2=l
			l=l+1		
			k=k+1
		END DO
	END DO
	link=(l-1)/2

	ALLOCATE(temp(n,m))		!rearrange coupling constant and 2-times-probability	
	temp%jj=0;
	temp%pp(1,1,1)=0;temp%pp(1,2,1)=0;temp%pp(2,1,1)=0;temp%pp(2,2,1)=0
	temp%pp(1,1,2)=0;temp%pp(1,2,2)=0;temp%pp(2,1,2)=0;temp%pp(2,2,2)=0
	DO i=1,n
		DO j=1,array(i)%degree		
		temp(i,j)=matrix(i,matrix(i,j)%na)		
		END DO
	END DO
	DO i=1,n
	DO j=1,m
	matrix(i,j)%jj=temp(i,j)%jj
!	matrix(i,j)%pp=temp(i,j)%pp
	END DO
	END DO
	DEALLOCATE(temp)

!	CALL init_prob_ER(array,array2,matrix,a,b,pp)
	CALL init_prob_ER2(array,array2,matrix,a,b,a1,b1)
!	CALL init_prob_ER3(array,array2,matrix,a,b,pp1)
END SUBROUTINE Voter_ER
!*******************************************************************************************************
!********************************INITIALIZATION MAX NUM OF CONFIG***************************************
!*******************************************************************************************************
SUBROUTINE initialization_config(N,conf)
TYPE(node_list),DIMENSION(:)::N
INTEGER,DIMENSION(:,:),ALLOCATABLE::conf
INTEGER,DIMENSION(:),ALLOCATABLE::temp
INTEGER::dmax,i,j,k,s,a
a=SIZE(N)
ALLOCATE(temp(a))
DO i=1,a
temp(i)=N(i)%degree
END DO
CALL ord(temp)
dmax=temp(1)+1
DEALLOCATE(conf)
ALLOCATE(conf(2**dmax,dmax))
DO j=1,dmax
s=1
DO i=1,2**dmax
conf(i,j)=s
IF (0==MODULO(i,2**(j-1))) THEN 
	s=-s
END IF
END DO
END DO
DEALLOCATE(temp)
!PRINT*,conf
END SUBROUTINE initialization_config

!********************************DIAMOND APPROXIMATION***************************************************
SUBROUTINE diamond_Voter_approx(node,neig,config_matrix,pin,tot_link,tol,tmin,tmax,u,v,vv)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER,DIMENSION(:,:),ALLOCATABLE::config_matrix
INTEGER::i,j,t,u,v,vv,tmin,tmax,tot_link
REAL*16::Pi,Pi2,tol,pin,mu,mu2,weight_magn,magn_in,tau,XI
REAL*16,DIMENSION(2)::magn,ro
WRITE (u,'(a50)') 'DIAMOND APPROXIMATIONS ON VOTER 1D'
WRITE (u,'(a25,f12.8)')'Initial probability Pi:', node(5)%p(1,2)

OPEN(UNIT=50, FILE='DataER/3Factorization_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=51, FILE='DataER/3Correlation_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=54, FILE='DataER/3Marginalization_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=55, FILE='DataER/3Correlation_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')


t=1;
CALL useful_parameter_diam(XI,mu,mu2,tau,magn_in,node)

magn=(/10.0,5.0/);ro=(/10.0,5.0/)!condition just to enter in the do while
!	DO WHILE (ABS(magn(1)-magn(2))>tol .AND. ABS(ro(1)-ro(2))>tol .AND. t<tmax )	
	DO WHILE (t<tmax)		!ABS(magn-weight_magn)>tol .AND.
	WRITE (50,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (54,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (55,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (51,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (u,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (v,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (v,'(8f12.8)') node(1:8)%p(1,3)
	WRITE (v,'(8f12.8)') node(1:8)%p(2,3)

	DO i=1,SIZE(node)
	CALL PD_III(i,node,neig,config_matrix)
	CALL print_diam_config(node,neig,i,u)
	CALL PD_II(i,node,neig,config_matrix)
	Pi=PD_I(i,1,node,neig,config_matrix)
	Pi2=PD_I(i,-1,node,neig,config_matrix)
	node(i)%p(1,1)=Pi/(Pi+Pi2)
	node(i)%p(2,1)=Pi2/(Pi+Pi2)
	END DO
	t=t+1
	CALL norm_and_timestep(node,neig)
	CALL print_correlation_control(node,neig,55)
	CALL print_marginalization_control(node,neig,54)
	CALL factorization_control(node,neig,50)
	CALL print_correlation_log(node,neig,51)
	CALL magnetization(node,magn,weight_magn,tot_link)
	CALL correlation_function(node,neig,ro(2),tot_link)
	CALL correlation_function_tdiff(node,neig,ro(1),tot_link)
	magn(1)=(magn_in - weight_magn)*EXP(-2*REAL(t)/(tau)) + weight_magn

PRINT*,'time',t
PRINT*,'Active link density t diff',ro(1)
PRINT*,'Active link density',ro(2)
PRINT*,'Active link density teo',XI*(1-magn(2))
PRINT*,'Weighted Magnetization',weight_magn
PRINT*,'Average Magnetization teo',magn(1)
PRINT*,'Average Magnetization',magn(2)
WRITE (vv,'(i8,f10.4,6f12.8)') t,REAL(t)/SIZE(node),magn(2),magn(1),weight_magn,ro(2),ro(1),XI*(1-magn(2))
END DO

PRINT*,'************************DIAMOND APPROXIMATIONS******************************'
PRINT*,'Final time',t
PRINT*,'Final Average Magnetization teo',magn(1)
PRINT*,'Final Average Magnetization',magn(2)
PRINT*,'Final Active link density',ro(2)
PRINT*,'Final Active link density t diff',ro(1)
PRINT*,'Final Weighted Magnetization',weight_magn
PRINT*,'*************************************************************************'
CLOSE(50);CLOSE(51);CLOSE(54);CLOSE(55)
END SUBROUTINE diamond_Voter_approx

REAL*16 FUNCTION PD_I(i,s,ND,ED,config) RESULT(Pdi)
	INTEGER::i,j,s,si,ss,k,kmax
	TYPE(node_list),DIMENSION(:)::ND
	TYPE(neighbour_list),DIMENSION(:,:)::ED
	INTEGER,DIMENSION(:,:)::config
	REAL*16::P
	P=1
	Pdi=0
	kmax=ND(i)%degree+1	
	IF (s==1)THEN
		ss=1
	ELSE 
		ss=2
	END IF
IF(kmax==1) THEN
	Pdi=ND(i)%p(ss,2)
ELSE
	DO si=1,2			
	DO k=1,2**kmax				
		DO j=1,kmax
		IF(config(k,j)==1)THEN
			ss=1
		ELSE
			ss=2
		END IF
		IF(j==kmax)THEN				
			P=P*ND(i)%pp(ss,si,2)			
		ELSE					
			P=P*ED(i,j)%pp(ss,si,2)
		END IF	
		END DO
	Pdi = Pdi + W_voter3(i,ND,ED,s,config(k,:kmax))*ND(i)%p(si,3)**(1-kmax)*P
	P=1
	END DO
	END DO
END IF
END FUNCTION PD_I

REAL*16 FUNCTION WI_voter3(i,jj,NW,EW,sw,swj,dw) RESULT(Wi)
	INTEGER::i,jj,sw,ssw,swj,sswj,k,l,j
	INTEGER,DIMENSION(:)::dw
	REAL*16::Wint,Wtemp
	TYPE(node_list),DIMENSION(:)::NW
	TYPE(neighbour_list),DIMENSION(:,:)::EW
	IF(sw==1)THEN
		ssw=+1
	ELSE 
		ssw=-1
	END IF
	IF(swj==1)THEN
		sswj=+1
	ELSE 
		sswj=-1
	END IF

	Wi=0;Wint=0.0;Wtemp=0.0
	l=1
	IF(jj==NW(i)%degree+1)THEN
		DO k=1,NW(i)%degree
			Wint=Wint+dw(k)
		END DO
	ELSE
		DO k=1,NW(i)%degree
			IF(k==jj)THEN
				Wint=Wint+sswj
			ELSE
				Wint=Wint+dw(l)
				l=l+1
			END IF
		END DO
	END IF
	IF(jj==NW(i)%degree+1)THEN
		IF(ssw == sswj)THEN
			Wtemp=(1-1.0/SIZE(NW))
			Wi=Wi + Wtemp
		END IF
	ELSE	
		IF(ssw == dw(NW(i)%degree))THEN
			Wtemp=(1-1.0/SIZE(NW))
			Wi=Wi + Wtemp
		END IF
	END IF
	Wtemp=(0.5/SIZE(NW))*(1+ssw*Wint/NW(i)%degree)
	Wi = Wi + Wtemp
END FUNCTION WI_voter3


SUBROUTINE PD_II(i,ND,ED,config) 
	INTEGER::i,j,l,jj,s,si,sj,ss,k,kmax
	TYPE(node_list),DIMENSION(:)::ND
	TYPE(neighbour_list),DIMENSION(:,:)::ED
	INTEGER,DIMENSION(:,:)::config
	REAL*16::P,Pi
	P=1
	Pi=0
	kmax=ND(i)%degree+1
IF(kmax.ne.1)THEN
DO jj=1,kmax
	DO s=1,2	!spin si al tempo t
	DO sj=1,2	!spin sj al tempo t-1 (il primo vicino che escludo)
		DO si=1,2	!spin si al tempo t-2
		DO k=1,2**(kmax-1)	!spin sj al tempo t-1 (tutti i primi vicini eccetto 1)
		l=1
			DO j=1,kmax					
				IF (j==jj)THEN
					IF(j==kmax)THEN			
						P=P*ND(i)%pp(sj,si,2)		
					ELSE				
						P=P*ED(i,jj)%pp(sj,si,2)
					END IF				
				ELSE
					IF(config(k,l).eq.1)THEN			
						ss=1
					ELSE
						ss=2
					END IF
					IF(j==kmax)THEN			
						P=P*ND(i)%pp(ss,si,2)		
					ELSE				
						P=P*ED(i,j)%pp(ss,si,2)			
					END IF				
					l=l+1
				END IF
			END DO
		Pi = Pi + WI_voter3(i,jj,ND,ED,s,sj,config(k,:kmax-1))*ND(i)%p(si,3)**(1-kmax)*P
		P=1
		END DO
		END DO
	IF(jj==kmax)THEN
		ND(i)%pp(s,sj,1)=Pi						
	ELSE								
		WHERE(ED(ED(i,jj)%na,:)%na==i) 
			ED(ED(i,jj)%na,:)%pp(s,sj,1)=Pi
		END WHERE
	END IF								
	Pi=0
	END DO 
	END DO
END DO
ELSE 
ND(i)%pp(:,:,1)=ND(i)%pp(:,:,2)						
END IF
END SUBROUTINE PD_II

SUBROUTINE PD_III(i,ND,ED,config) 
	INTEGER::i,j,l,jj,s,si,sii,sj,ss,k,kmax
	TYPE(node_list),DIMENSION(:)::ND
	TYPE(neighbour_list),DIMENSION(:,:)::ED
	INTEGER,DIMENSION(:,:)::config
	REAL*16::P,Pi
	P=1
	Pi=0
	kmax=ND(i)%degree+1
DO jj=1,kmax-1
DO sii=1,2 		!spin si al tempo t-1
DO sj=1,2		!spin sj al tempo t-1 (il primo vicino che escludo)
	DO s=1,2	!spin si al tempo t
	DO si=1,2	!spin si al tempo t-2
		DO k=1,2**(kmax-2)	!spin sj al tempo t-1 (tutti i primi vicini eccetto 1)
			l=1
			DO j=1,kmax					
				IF (j==jj)THEN
					P=P*ED(i,jj)%pp(sj,si,2)
				ELSE
					IF(config(k,l).eq.1)THEN			
						ss=1
					ELSE
						ss=2
					END IF
					IF(j==kmax)THEN			
						P=P*ND(i)%pp(sii,si,2)		
					ELSE				
						P=P*ED(i,j)%pp(ss,si,2)			
						l=l+1
					END IF				
				END IF
			END DO
		Pi = Pi + WII_voter3(i,jj,ND,ED,s,sii,sj,config(k,:kmax-2))*ND(i)%p(si,3)**(1-kmax)*P
		P=1
		END DO
	END DO
	END DO 
!ED(i,jj)%pp_t(sii,sj)=Pi
	WHERE(ED(ED(i,jj)%na,:)%na==i) 
		ED(ED(i,jj)%na,:)%pp_t(sj,sii)=Pi
	END WHERE
	Pi=0
END DO
END DO
END DO
END SUBROUTINE PD_III

REAL*16 FUNCTION WII_voter3(i,jj,NW,EW,sw,swi,swj,dw) RESULT(Wi)
	INTEGER::i,jj,sw,ssw,swj,sswj,k,l,j,swi,sswi
	INTEGER,DIMENSION(:)::dw
	REAL*16::Wint,Wtemp
	TYPE(node_list),DIMENSION(:)::NW
	TYPE(neighbour_list),DIMENSION(:,:)::EW
	IF(sw==1)THEN
		ssw=+1
	ELSE 
		ssw=-1
	END IF
	IF(swi==1)THEN
		sswi=+1
	ELSE 
		sswi=-1
	END IF
	IF(swj==1)THEN
		sswj=+1
	ELSE 
		sswj=-1
	END IF

	Wi=0;Wint=0.0;Wtemp=0.0
	l=1
	DO k=1,NW(i)%degree
		IF(k==jj)THEN
			Wint=Wint+sswj
		ELSE
			Wint=Wint+dw(l)
			l=l+1
		END IF
	END DO
	IF(ssw == sswi)THEN
		Wtemp=(1-1.0/SIZE(NW))
		Wi=Wi + Wtemp
	END IF
	Wtemp=(0.5/SIZE(NW))*(1+ssw*Wint/NW(i)%degree)
	Wi = Wi + Wtemp
END FUNCTION WII_voter3


REAL*16 FUNCTION W_voter3(i,NW,EW,sw,dw) RESULT(Wi)	!NON NORMALIZZATO CON 1 - 1/N 
	INTEGER::i,sw,k
	REAL*16::Wint
	INTEGER,DIMENSION(:)::dw
	TYPE(node_list),DIMENSION(:)::NW
	TYPE(neighbour_list),DIMENSION(:,:)::EW
	
	Wi=0;Wint=0
	DO k=1,NW(i)%degree
		Wint = Wint + dw(k)
	END DO
	IF(sw == dw(NW(i)%degree+1))THEN
		Wi=Wi + (1-1.0/SIZE(NW))
	END IF
	Wi = Wi + (0.5/SIZE(NW))*(1+sw*Wint/NW(i)%degree)
END FUNCTION W_voter3
!************************************STAR APPROXIMATION**************************************************
REAL*16 FUNCTION P_I(i,s,NP,EP,config) RESULT(Pi)
	INTEGER::i,j,s,ss,k,kmax
	TYPE(node_list),DIMENSION(:)::NP
	TYPE(neighbour_list),DIMENSION(:,:)::EP
	INTEGER,DIMENSION(:,:)::config
	REAL*16::P
	P=1
	Pi=0
	kmax=NP(i)%degree+1
IF(kmax==1) THEN
	IF (s==1)THEN
		ss=1
	ELSE 
		ss=2
	END IF
	Pi = NP(i)%p(ss,3)
ELSE
	DO k=1,2**kmax	
		DO j=1,kmax
		IF(config(k,j)==1)THEN
			ss=1
		ELSE
			ss=2
		END IF
		IF(j==kmax)THEN					
			P=P*NP(i)%p(ss,3)			
		ELSE						
			P=P*NP(EP(i,j)%na)%p(ss,3)
		END IF		
		END DO
	Pi = Pi + W_voter3(i,NP,EP,s,config(k,:kmax))*P
	P=1
	END DO
END IF
END FUNCTION P_I

SUBROUTINE P_I2(i,NP,EP,config)
	INTEGER::i,j,s,ss,si,sj,jj,k,kmax,l
	TYPE(node_list),DIMENSION(:)::NP
	TYPE(neighbour_list),DIMENSION(:,:)::EP
	INTEGER,DIMENSION(:,:)::config
	REAL*16::P,Pi
	P=1
	Pi=0
	kmax=NP(i)%degree+1
DO jj=1,kmax-1
	DO si=1,2	!si a t-1
	DO sj=1,2	!sj a t-1 (vicino escluso)
		DO s=1,2	! si a t
		DO k=1,2**(kmax-2)	!sj a t-1 (tutti tranne l'escluso)
			l=1
			DO j=1,kmax
			IF (j==jj)THEN
				P=P*NP(EP(i,j)%na)%p(sj,3)
			ELSE
				IF(config(k,l).eq.1)THEN			
					ss=1
				ELSE
					ss=2
				END IF
				IF(j==kmax)THEN					
					P=P*NP(i)%p(si,3)			
				ELSE						
					P=P*NP(EP(i,j)%na)%p(ss,3)
				l=l+1
				END IF		
			END IF
			END DO
		Pi = Pi + WII_voter3(i,jj,NP,EP,s,si,sj,config(k,:kmax-2))*P
		P=1
		END DO
		END DO
	WHERE(EP(EP(i,jj)%na,:)%na==i) 
		EP(EP(i,jj)%na,:)%pp_t(sj,si)=Pi
	END WHERE
	Pi=0
	END DO
	END DO
END DO
END SUBROUTINE P_I2


SUBROUTINE star_Voter_approx(node,neig,config_matrix,pin,tot_link,tol,tmin,tmax,u,v,vv)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER,DIMENSION(:,:),ALLOCATABLE::config_matrix
INTEGER::i,j,t,u,v,vv,tmin,tmax,tot_link,NK
REAL*16::Pi,Pi2,tol,pin,mu,XI,weight_magn,magn_in,mu2,tau,mk,mk_in
REAL*16,DIMENSION(2)::magn,sigma,ro
OPEN(UNIT=52, FILE='DataER/3Correlation_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=53, FILE='DataER/3Correlation_control_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
WRITE (u,'(a50)') 'STAR APPROXIMATION ON VOTER MODEL ON ER GRAPH'

	CALL useful_parameter_star(XI,mu,mu2,tau,magn_in,node)
	t=1
	magn=(/ 10 , 5 /);ro=(/ 10 , 5 /);
	DO WHILE (t<tmax)	!ABS(magn(2) - weight_magn)>tol .AND. 
		WRITE (53,'(a3,i6,a50)') 't=',t, '    *************************************'
		WRITE (52,'(a3,i6,a50)') 't=',t, '    *************************************'
		WRITE (u,'(a3,i6,a50)') 't=',t, '    *************************************'
		WRITE (u,'(8f12.8)') node(1:8)%p(1,3)
		WRITE (u,'(8f12.8)') node(1:8)%p(2,3)
		DO i=1,SIZE(node)
			Pi=P_I(i,1,node,neig,config_matrix)
			Pi2=P_I(i,-1,node,neig,config_matrix)
			node(i)%p(1,2)=Pi/(Pi+Pi2)
			node(i)%p(2,2)=Pi2/(Pi+Pi2)
		CALL P_I2(i,node,neig,config_matrix)
		END DO
		CALL print_correlation_control(node,neig,53)
	neig%pp_t(1,1)=neig%pp_t(1,1)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(2,1)=neig%pp_t(2,1)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(1,2)=neig%pp_t(1,2)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(2,2)=neig%pp_t(2,2)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
		node%p(1,3)=node%p(1,2)
		node%p(2,3)=node%p(2,2)
		t=t+1
		CALL print_correlation_log(node,neig,52)	
		CALL correlation_function(node,neig,ro(2),tot_link)
		CALL magnetization(node,magn,weight_magn,tot_link)
		magn(1)=0
		magn(1)=(magn_in - weight_magn)*EXP(-2*REAL(t)/(tau)) + weight_magn
		ro(1)=XI*(1-magn(2))
PRINT*,'time',t
PRINT*,'Weighted Magnetization',weight_magn
PRINT*,'Average Magnetization',magn(2)
PRINT*,'Average Magnetization teo',magn(1)
PRINT*,'Active link density',ro(2)
PRINT*,'Active link density 2 - metastable',ro(1)
	WRITE (v,'(i6,f12.8)') t,node(5)%p(1,2)
	WRITE (vv,'(i6,7f12.8)') t,REAL(t)/SIZE(node),pin,magn(2),magn(1),weight_magn,ro(2),ro(1)
	END DO
PRINT*,'************************STAR APPROXIMATIONS******************************'
PRINT*,'Final time',t
PRINT*,'Final Average Magnetization',magn(2)
PRINT*,'Final Average Magnetization teo',magn(1)
PRINT*,'Final Weighted Magnetization',weight_magn
PRINT*,'Final Active link density',ro(2)
PRINT*,'Active link density 2 - metastable',ro(1)
PRINT*,'*************************************************************************'
CLOSE(52);CLOSE(53)
END SUBROUTINE star_Voter_approx


!************************************USEFUL FUNCTION*************************************
SUBROUTINE print_diam_config(node,neig,i,u)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,u
	WRITE (u,'(a10,i4)') 'NODE:', i
	WRITE (u,'(2a12)') 'P(si=+1)','P(si=-1)'
	WRITE (u,'(2f12.8)') node(i)%p(1,3),node(i)%p(2,3)
	WRITE (u,'(a15)') 'P(si=+1,sj=+1)'
	WRITE (u,'(7f12.8)') neig(i,1:7)%pp(1,1,2)
	WRITE (u,'(a15)') 'P(si=-1,sj=+1)'
	WRITE (u,'(7f12.8)') neig(i,1:7)%pp(1,2,2)
	WRITE (u,'(a15)') 'P(si=+1,sj=-1)'
	WRITE (u,'(7f12.8)') neig(i,1:7)%pp(2,1,2)
	WRITE (u,'(a15)') 'P(si=-1,sj=-1)'
	WRITE (u,'(7f12.8)') neig(i,1:7)%pp(2,2,2)
END SUBROUTINE print_diam_config

SUBROUTINE print_correlation_log(node,neig,uu)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,u,uu
DO i=1,SIZE(node)
	WRITE (uu,'(a10,i4)') 'NODE:', i
!	WRITE (uu,'(2a12)') 'P(si=+1)','P(si=-1)'
!	WRITE (uu,'(2f12.8)') node(i)%p(1,3),node(i)%p(2,3)
	WRITE (uu,'(a15)') 'P(si=+1,sj=+1)'
	WRITE (uu,'(7f12.8)') neig(i,1:7)%pp_t(1,1)
	WRITE (uu,'(a15)') 'P(si=-1,sj=+1)'
	WRITE (uu,'(7f12.8)') neig(i,1:7)%pp_t(1,2)
	WRITE (uu,'(a15)') 'P(si=+1,sj=-1)'
	WRITE (uu,'(7f12.8)') neig(i,1:7)%pp_t(2,1)
	WRITE (uu,'(a15)') 'P(si=-1,sj=-1)'
	WRITE (uu,'(7f12.8)') neig(i,1:7)%pp_t(2,2)
END DO
END SUBROUTINE print_correlation_log



SUBROUTINE magnetization(node,magn,weight_magn,tot_link)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,tot_link
REAL*16::weight_magn
REAL*16,DIMENSION(:)::magn
magn(2)=0;weight_magn=0
DO i=1,SIZE(node)
	magn(2)=magn(2)+(node(i)%p(1,2)-node(i)%p(2,2))
	weight_magn=weight_magn+node(i)%degree*(node(i)%p(1,2)-node(i)%p(2,2))
END DO
weight_magn=weight_magn/(2*tot_link)
magn(2)=magn(2)/SIZE(node)
END SUBROUTINE magnetization

SUBROUTINE norm_and_timestep(node,neig)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
	node%pp(1,1,2)=node%pp(1,1,1)/(node%pp(1,1,1)+node%pp(2,1,1)+node%pp(1,2,1)+node%pp(2,2,1))
	node%pp(2,1,2)=node%pp(2,1,1)/(node%pp(1,1,1)+node%pp(2,1,1)+node%pp(1,2,1)+node%pp(2,2,1))
	node%pp(1,2,2)=node%pp(1,2,1)/(node%pp(1,1,1)+node%pp(2,1,1)+node%pp(1,2,1)+node%pp(2,2,1))
	node%pp(2,2,2)=node%pp(2,2,1)/(node%pp(1,1,1)+node%pp(2,1,1)+node%pp(1,2,1)+node%pp(2,2,1))

	neig%pp(1,1,2)=neig%pp(1,1,1)/(neig%pp(1,1,1)+neig%pp(2,1,1)+neig%pp(1,2,1)+neig%pp(2,2,1))
	neig%pp(2,1,2)=neig%pp(2,1,1)/(neig%pp(1,1,1)+neig%pp(2,1,1)+neig%pp(1,2,1)+neig%pp(2,2,1))
	neig%pp(1,2,2)=neig%pp(1,2,1)/(neig%pp(1,1,1)+neig%pp(2,1,1)+neig%pp(1,2,1)+neig%pp(2,2,1))
	neig%pp(2,2,2)=neig%pp(2,2,1)/(neig%pp(1,1,1)+neig%pp(2,1,1)+neig%pp(1,2,1)+neig%pp(2,2,1))

	neig%pp_t(1,1)=neig%pp_t(1,1)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(2,1)=neig%pp_t(2,1)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(1,2)=neig%pp_t(1,2)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))
	neig%pp_t(2,2)=neig%pp_t(2,2)/(neig%pp_t(1,1)+neig%pp_t(2,1)+neig%pp_t(1,2)+neig%pp_t(2,2))

	node%p(1,3)=node%p(1,2)
	node%p(2,3)=node%p(2,2)
	node%p(1,2)=node%p(1,1)
	node%p(2,2)=node%p(2,1)
END SUBROUTINE norm_and_timestep

SUBROUTINE print_marginalization_control(node,neig,u)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,u
REAL*16::a,b,d,e,f

WRITE(u,'(2a6,4a20)') 'node i','node j','sum P(sj,si)','P(sj^t-1/si^t-2)','Difference','Norm 2-t/1-t'

DO i=1,SIZE(node)
DO j=1,node(i)%degree
a=neig(i,j)%pp(1,1,2)+neig(i,j)%pp(1,2,2)
b=neig(i,j)%pp(1,1,2)+neig(i,j)%pp(2,1,2)
d=node(neig(i,j)%na)%p(1,2)
e=neig(i,j)%pp(1,1,2)+neig(i,j)%pp(1,2,2)+neig(i,j)%pp(2,1,2)+neig(i,j)%pp(2,2,2)
f=node(i)%p(1,3)+node(i)%p(2,3)
WRITE(u,'(2i6,4f20.8)') i,neig(i,j)%na,a,d,d-a,e
WRITE(u,'(2i6,4f20.8)') i,i,b,node(i)%p(1,3),node(i)%p(1,3)-b,f
END DO
END DO

END SUBROUTINE print_marginalization_control

SUBROUTINE factorization_control(node,neig,u)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,u
REAL*16::a,b

WRITE (u,'(2a6,3a30)')'i','j','P(sj^t-1=+1,si^t-2=+1)','P(si^t-2=+1)*P(sj^t-1=+1)','Difference'
DO i=1,SIZE(node)
DO j=1,node(i)%degree
a=neig(i,j)%pp(1,1,1)
b=node(i)%p(1,3)*node(neig(i,j)%na)%p(1,2)
WRITE ( 50,'(2i6,3f30.10)')i,neig(i,j)%na,a,b,a-b
END DO
END DO
END SUBROUTINE factorization_control

SUBROUTINE print_correlation_control(node,neig,u)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,u
REAL*16::a,b,d,e
WRITE(u,'(2a6,4a15)') 'node i','node j','sum corr','P(sj^t-1)','Difference','Norm'
DO i=1,SIZE(node)
DO j=1,node(i)%degree
a=neig(i,j)%pp_t(1,1)+neig(i,j)%pp_t(2,1)
b=neig(i,j)%pp_t(1,1)+neig(i,j)%pp_t(1,2)
d=node(neig(i,j)%na)%p(1,3)
e=neig(i,j)%pp_t(1,1)+neig(i,j)%pp_t(2,1)+neig(i,j)%pp_t(1,2)+neig(i,j)%pp_t(2,2)
WRITE(u,'(2i6,4f15.8)') i,neig(i,j)%na,a,d,d-a,e
WRITE(u,'(2i6,3f15.8)') i,i,b,node(i)%p(1,3),node(i)%p(1,3)-b
END DO
END DO

END SUBROUTINE print_correlation_control

SUBROUTINE correlation_function_tdiff(node,neig,ro,tot_link)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,si,sj,ssi,ssj,tot_link
REAL*16::ro
REAL*16::ro_temp
ro=0;!ro_temp=0
DO i=1,SIZE(node)
	DO j=1,node(i)%degree+1
	ro_temp=0
		DO si=-1,1,2
		DO sj=-1,1,2
			IF (si==1)THEN
			ssi=1
			ELSE
			ssi=2
			END IF
			IF (sj==1)THEN
			ssj=1
			ELSE
			ssj=2
			END IF
			IF(j==node(i)%degree+1)THEN
			ro_temp=ro_temp+si*sj*node(i)%pp(ssi,ssj,2)
			ELSE
			ro_temp=ro_temp+si*sj*neig(i,j)%pp(ssi,ssj,2)
			END IF
!PRINT*,'si,sj',si,sj,ro_temp
		END DO
		END DO
	ro=ro+ro_temp
!PRINT*,'j',j,ro_temp,ro
	END DO
!PRINT*,'i',i,ro_temp
END DO
ro=(1-ro/(2*tot_link))/2
!PRINT*,ro
END SUBROUTINE correlation_function_tdiff


SUBROUTINE correlation_function(node,neig,ro,tot_link)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,si,sj,ssi,ssj,tot_link
REAL*16::ro
REAL*16::ro_temp
ro=0
DO i=1,SIZE(node)
	DO j=1,node(i)%degree
	ro_temp=0
		DO si=-1,1,2
		DO sj=-1,1,2
			IF (si==1)THEN
			ssi=1
			ELSE
			ssi=2
			END IF
			IF (sj==1)THEN
			ssj=1
			ELSE
			ssj=2
			END IF
			ro_temp=ro_temp+si*sj*neig(i,j)%pp_t(ssi,ssj)
!PRINT*,'si,sj',si,sj,ro_temp
		END DO
		END DO
!PRINT*,'j',j,ro_temp,ro(2)
	ro=ro+ro_temp
	END DO
!PRINT*,'i',i,ro(2)
END DO
ro=(1-ro/(2*tot_link))/2
END SUBROUTINE correlation_function

SUBROUTINE useful_parameter_diam(XI,mu,mu2,tau,magn_in,node)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
REAL*16::mu,mu2,tau,magn_in,XI
INTEGER::i
	mu=degree_distribution(node)
	mu2=second_moment(node)
	XI=(mu-2)/(2*(mu-1))
	tau=(((mu-1)*mu**2)*SIZE(node))/((mu-2)*mu2)
	DO i=1,SIZE(node)		!initial magnetization
		magn_in=magn_in+(node(i)%p(1,2)-node(i)%p(2,2))
	END DO
	magn_in=magn_in/SIZE(node)
END SUBROUTINE useful_parameter_diam


SUBROUTINE useful_parameter_star(XI,mu,mu2,tau,magn_in,node)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
REAL*16::mu,mu2,tau,magn_in,XI
INTEGER::i
	mu=degree_distribution(node)
	mu2=second_moment(node)
	XI=(mu-2)/(2*(mu-1))
	tau=(((mu-1)*mu**2)*SIZE(node))/((mu-2)*mu2)
	DO i=1,SIZE(node)		!initial magnetization
		magn_in=magn_in+(node(i)%p(1,3)-node(i)%p(2,3))
	END DO
	magn_in=magn_in/SIZE(node)
END SUBROUTINE useful_parameter_star

REAL*16 FUNCTION degree_distribution(node) RESULT(mu) 
	TYPE(node_list),DIMENSION(:)::node
	REAL,DIMENSION(:),ALLOCATABLE::NK
	INTEGER::i

	ALLOCATE(NK(SIZE(node)))
	NK=0
	DO i=1,SIZE(node)
		NK(node(i)%degree+1)=NK(node(i)%degree+1)+1
	END DO

	NK=NK/SIZE(node)
	mu=0
	DO i=1,SIZE(node)
		mu=mu+(i-1)*NK(i)
	END DO
	DEALLOCATE(NK)
END FUNCTION degree_distribution

REAL*16 FUNCTION second_moment(node) RESULT(mu2) 
	TYPE(node_list),DIMENSION(:)::node
	REAL,DIMENSION(:),ALLOCATABLE::NK
	INTEGER::i

	ALLOCATE(NK(SIZE(node)))
	NK=0
	DO i=1,SIZE(node)
		NK(node(i)%degree+1)=NK(node(i)%degree+1)+1
	END DO

	NK=NK/SIZE(node)
	mu2=0
	DO i=1,SIZE(node)
		mu2=mu2+((i-1)**2)*NK(i)
	END DO
	DEALLOCATE(NK)
END FUNCTION second_moment

!********************************************************************************************************
!**************************************VOTER DYNAMICS****************************************************
!********************************************************************************************************
SUBROUTINE VoterER_step(VN,VE,M)
	INTEGER::i,j,M,M1
	INTEGER*8::seed=0
	REAL*16::x
	TYPE(node_list),DIMENSION(:)::VN
	TYPE(neighbour_list),DIMENSION(:,:)::VE
	CALL random_number(x)
	i=1+SIZE(VN)*x!lcg(seed) !x
IF(VN(i)%degree/=0)THEN
	CALL random_number(x)
	j=1+VN(i)%degree*x!lcg(seed) !x
	IF( VN(i)%s(1) /= VN(VE(i,j)%na)%s(1) )THEN
!PRINT*,i,j,VN(i)%s(1),M
		VN(i)%s(1)=VN(VE(i,j)%na)%s(1)	!direct voter model
		M=M+2*VN(i)%s(1)
!		VN(VE(i,j)%na)%s(1)=VN(i)%s(1)	!reverse voter model
!		M=M+2*VN(VE(i,j)%na)%s(1)
!PRINT*,VN(i)%s(1),M
	END IF
!ELSE
END IF
!M1=0				!VERIFICA MAGNETIZZAZIONE
!DO i=1,SIZE(VN)
!M1=M1+VN(i)%s(1)
!END DO
!PRINT*,M1
!PRINT*,'****************************'
END SUBROUTINE VoterER_step

SUBROUTINE VoterER_step_linkupdate(VN,VE,link,M)
	INTEGER::i,j,k,l,l2,link,M,M1
	INTEGER*8::seed=0
	REAL*16::x
	TYPE(node_list),DIMENSION(:)::VN
	TYPE(neighbour_list),DIMENSION(:,:)::VE

	CALL random_number(x)
	l=1+link*x
	CALL random_number(x)
	l2=1+2*x

k=1
DO i=1,SIZE(VN)
DO j=1,VN(i)%degree
	k=k+1
	IF ( k == l )THEN
		IF(VN(i)%s(1)/=VN(VE(i,j)%na)%s(1))THEN
			IF( l2==1 )THEN
!PRINT*,i,j,VN(i)%s(1),M
				VN(i)%s(1)=VN(VE(i,j)%na)%s(1)
				M=M+2*VN(i)%s(1)
!PRINT*,VN(i)%s(1),M
				EXIT
			ELSE
!PRINT*,i,j,VN(VE(i,j)%na)%s(1),M
				VN(VE(i,j)%na)%s(1)=VN(i)%s(1)
				M=M+2*VN(VE(i,j)%na)%s(1)
!PRINT*,VN(VE(i,j)%na)%s(1),M
				EXIT
			END IF	
		END IF
	END IF
END DO
END DO
!M1=0				!VERIFICA MAGNETIZZAZIONE
!DO i=1,SIZE(VN)
!M1=M1+VN(i)%s(1)
!END DO
!PRINT*,M1
!PRINT*,'*********************'
END SUBROUTINE VoterER_step_linkupdate

!***************************MAGNETIZATION DATA FOR VOTER DYNAMICS**********************************
SUBROUTINE Magn_VoterER(NG,NE,pin,tol,tmin,tmax,tmax2,u,uu)
	TYPE(node_list),DIMENSION(:)::NG
	TYPE(neighbour_list),DIMENSION(:,:)::NE
	REAL*16::tol,M_av,M_av2,M_temp,M_temp2,pin
	INTEGER::i,j,t,t1,t2,t3,u,uu,tmin,tmax,tmax2,M,link

	M_av=0.0;M_av2=0.0;t=1;t1=1;t2=1;t3=1;M_temp=0;M_temp2=0;link=0
	DO i=1,SIZE(NG)
		link=link+NG(i)%degree
	END DO

	DO WHILE (t2<tmax)! ABS(M_av-M_av2)>tol .AND.			
		CALL init_voterER_MC(NG)
		M=0
		DO i=1,SIZE(NG)	!compute initial magnetization
			M = M + NG(i)%s(1)
		END DO
		t1=1
PRINT*,'num sim',t2-1,'correct',t-1,'failed',t3-1
PRINT*,'initial magnetization',M
		DO WHILE ( ABS(M) /= SIZE(NG) .AND. t1 < tmax2)
			CALL VoterER_Step(NG,NE,M)
!			CALL VoterER_Step_linkupdate(NG,NE,link,M)
			t1 = t1 + 1	
		END DO
		M_temp2= M_temp2 + REAL(M)/SIZE(NG)
		M_av2= M_temp2 / t2
PRINT*,'final general magnetization',M_av2
		t2=t2+1
		IF(t1<tmax2)THEN 	
			M_temp= M_temp + REAL(M)/SIZE(NG)
			M_av= M_temp / t
PRINT*,'final correct magnetization',M_av	
			WRITE (u,'(2i10,2f15.10)') t,t1,REAL(M)/SIZE(NG),M_av
			t=t+1
		ELSE
			t3=t3+1
			IF(t3>tmax**2)THEN
				PRINT*,'consensus never reached'
				EXIT
			END IF
		END IF
	END DO
WRITE (uu,'(f15.10,i10,f15.10,i10,f15.10)') Pin,t2,M_av2,t,M_av
END SUBROUTINE Magn_VoterER

SUBROUTINE VoterER_Step_1sim(VN,VE,M,ro,pp)
	INTEGER::i,j,k,M,M1,pp,ro
	INTEGER*8::seed=0
	REAL*16::x
	TYPE(node_list),DIMENSION(:)::VN
	TYPE(neighbour_list),DIMENSION(:,:)::VE
	CALL random_number(x)
	i=1+SIZE(VN)*x!lcg(seed) !x
IF(VN(i)%degree/=0)THEN
	CALL random_number(x)
	j=1+VN(i)%degree*x!lcg(seed) !x
	IF( VN(i)%s(1) /= VN(VE(i,j)%na)%s(1) )THEN
!PRINT*,i,j,VN(i)%s(1),M
		VN(i)%s(1)=VN(VE(i,j)%na)%s(1)	!direct voter model
		M=M+2*VN(i)%s(1)
		pp = pp + VN(i)%s(1)
		DO k=1,VN(i)%degree
			IF(VN(i)%s(1) == VN(VE(i,k)%na)%s(1))THEN
				ro=ro-2
			ELSE
				ro=ro+2
			END IF
		END DO
!		VN(VE(i,j)%na)%s(1)=VN(i)%s(1)	!reverse voter model
!		M=M+2*VN(VE(i,j)%na)%s(1)
!PRINT*,VN(i)%s(1),M
	END IF
END IF
!M1=0				!VERIFICA MAGNETIZZAZIONE
!DO i=1,SIZE(VN)
!M1=M1+VN(i)%s(1)
!END DO
!PRINT*,M1
!PRINT*,'****************************'
END SUBROUTINE VoterER_step_1sim


SUBROUTINE Magn_VoterER_1sim(NG,NE,tmax2,u)
	TYPE(node_list),DIMENSION(:)::NG
	TYPE(neighbour_list),DIMENSION(:,:)::NE
	REAL*16::tol,M_av,ro2,p2,XI,XI2,mu
	INTEGER::i,j,t,u,tmax2,M,link,ro,p

	M_av=0.0;t=1;link=0;M=0;p=0;ro=0;p2=0;ro2=0;M_av=0
	DO i=1,SIZE(NG)
		link=link+NG(i)%degree
	END DO
	mu=degree_distribution(NG)
	CALL init_voterER_MC(NG)
	DO i=1,SIZE(NG)	!compute initial magnetization and active link
		M = M + NG(i)%s(1)
		IF(NG(i)%s(1) == 1 ) THEN	! P(si=+1) and P(si=-1)
			p = p + 1
		END IF
		DO j=1,NG(i)%degree
			IF(NG(i)%s(1) .ne. NG(NE(i,j)%na)%s(1) ) THEN
			ro=ro+1
			END IF
		END DO
	END DO
	p2 = REAL(p) / SIZE(NG)
	ro2 = REAL(ro) / link
	M_av= REAL(M) / SIZE(NG)
	XI2 = (2.0*(mu-2))/(mu-1)
PRINT*,'********************INITIAL CONFIGURATION*******************'
PRINT*,'active link',ro2
PRINT*,'magn',M
PRINT*,'average magn',M_av
PRINT*,'P(si=+1)',p,'P(si=-1)',SIZE(NG)-p
PRINT*,'************************************************************'

	WRITE (u,'(a10,6a12)') 'time','XI_emp','XI_teo','active link','P(si=+1)','P(si=-1)','Av. Magn'
	DO WHILE ( ABS(M) /= SIZE(NG) .AND. t < tmax2)
		XI = ro2 / (p2*(1-p2))

!		IF (0==MODULO(t,100)) THEN 
			WRITE (u,'(i10,6f12.8)') t,XI,XI2,ro2,p2,1-p2,M_av
PRINT*,'time',t
PRINT*,'XI',XI
PRINT*,'active link',ro2
PRINT*,'magn',M
PRINT*,'average magn',M_av
PRINT*,'P(si=+1)',p2
PRINT*,'*******************************************'
!		END IF
		CALL VoterER_Step_1sim(NG,NE,M,ro,p)
!		CALL VoterER_Step_linkupdate(NG,NE,link,M)
		t = t + 1	
		p2 = REAL(p) / SIZE(NG)
		ro2 = REAL(ro) / link
		M_av= REAL(M) / SIZE(NG)
	END DO
END SUBROUTINE Magn_VoterER_1sim
!******************************************************************************************************
!********************************PRINT NETWORK CONFIGURATION*******************************************
!****************************************************************************************************** 
SUBROUTINE config_print(node,neig,link,mu_in,u)
TYPE(node_list),DIMENSION(:)::node
TYPE(neighbour_list),DIMENSION(:,:)::neig
INTEGER::i,n,u,j,link
REAL*16::mu_in,mu,mu2,XI
n=SIZE(node)
mu=degree_distribution(node)
mu2=second_moment(node)
XI=(mu-2)/(2*(mu-1))
WRITE (u,'(a20,i6)')'number of link',link
WRITE (u,'(2(a20,f12.8))')'degree in',mu_in,'degree out', mu
WRITE (u,'(2(a20,f12.8))')'2nd moment',mu2,'tau',(((mu-1)*mu**2)*SIZE(node))/((mu-2)*mu2)
WRITE (u,'(a20,f12.8)')'Xi=(mu-2)/(2*(mu-1))',XI
WRITE (u,'(a50)')'**********************************************'

DO i=1,n
WRiTE (u,'(a15,i6,2(a15,f8.5))')'node',node(i)%names,'P(si_t-2=+1)',node(i)%p(1,3),'P(si_t-2=-1)',node(i)%p(2,3)
WRiTE (u,'(a15,i6,2(a15,f8.5))')'node',node(i)%names,'P(si_t-1=+1)',node(i)%p(1,2),'P(si_t-1=-1)',node(i)%p(2,2)
WRITE (u,'(a15)') 'neighbour list'
	DO j=1,node(i)%degree
	WRITE (u,'(a4,i4,a24,4f8.5)') 'node',neig(i,j)%na,'P(++),P(-+),P(+-),P(--)',neig(i,j)%pp(:,:,2)	
	END DO	
WRITE (u,'(a4,i4,a24,4f8.5)') 'node',i,'P(++),P(-+),P(+-),P(--)',node(i)%pp(:,:,2)	
WRITE (u,'(a20,i4)')'Degree of the node:',node(i)%degree
WRITE (u,'(a50)')'**********************************************'
END DO
END SUBROUTINE config_print

END MODULE network


!********************************************************************************************************
!********************************************************************************************************
!********************************************MAIN PROGRAM************************************************
!********************************************************************************************************
!********************************************************************************************************

PROGRAM VoterER
USE network
IMPLICIT NONE
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node_star,node_diam,node_MC
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig_star,neig_diam,neig_MC
INTEGER,DIMENSION(:,:),ALLOCATABLE::config_matrix
INTEGER::n,m,u,uu,uuu,v,vv,vvv,w,ww,www,z,zz,tmin,tmax,tmax2,k,deg,prob,link,tmax1sim
REAL*16::h_in,j_in,tol,p_in,mu,a,b,a1,b1,pp
REAL*16,DIMENSION(4)::pp1
link=0
h_in=1;j_in=1				!initialization config!!!!
n=5;	 				!initialization number of nodes and neighbours
tmin=1; tmax=1000			!minimum and max number of time step
tol=1E-3				!tolerance
tmax2=5*n**4				!max number for montecarlo
tmax1sim=1000
mu=3					!average degree
p_in=1-SQRT(mu/n)			!p_in to have or not a link
a=0.0d0;b=1d0				!range for p(si_t-2=+1) at random
a1=0.0d0;b1=1d0				!range for p(si_t-2=+1) at random
pp=0.5d0
pp1=(/.1d0,0.1d0,0.1d0,0.7d0/)		!p(si_t-2=+1)
!**************************DIAMOND APPROXIMATIONS ON ERDOS RENYI GRAPH***********************************
u=34;uu=35;uuu=36;v=37;vv=38;vvv=39;w=40;ww=41;www=42;z=43;zz=44
OPEN(UNIT=u, FILE='DataER/3VotDirER_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=uu, FILE='DataER/3Psi_VotDirER_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=uuu, FILE='DataER/3Magn_VotDirER_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')

OPEN(UNIT=v, FILE='DataER/3VotDirER_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=vv, FILE='DataER/3Psi_VotDirER_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=vvv, FILE='DataER/3Magn_VotDirER_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')

OPEN(UNIT=w, FILE='DataER/3MagnLog_VotDirER_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=ww, FILE='DataER/3Magn_VotDirER_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=www, FILE='DataER/3Voter_1SimVOT_n5.dat',STATUS='REPLACE',ACTION='WRITE')

OPEN(UNIT=z, FILE='DataER/3data_VotDirER_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=zz, FILE='DataER/3data_VotDirER_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')


WRITE (uuu,'(a6,7a12)')'t','t/N','Av M','Av M teo','Weight Magn','Act Link','Link t diff', 'Link teo'
WRITE (vvv,'(a6,7a12)') 't','t/N   ','Plink  ','Average Magn','Av Magn teo','Weight Magn','Active link','Link teo'
WRITE (w,'(2a10,2a15)')'#sim','Exit time','Final Magn','Average magn'
WRITE (ww,'(a15,a10,a15,a10,a15)') 'Plink   ','num sim','MagnTot   ','n_sim corr','Magn correct'


ALLOCATE(config_matrix(1,1))
PRINT*,'****************************************************************************'
PRINT*,'INITIAL AVERAGE PROBABILITY',(a+b)/2
PRINT*,'INITIAL Plink',mu,p_in
PRINT*,'****************************************************************************'
	ALLOCATE(node_star(n));ALLOCATE(node_diam(n));ALLOCATE(node_MC(n));
	ALLOCATE(neig_star(n,n));ALLOCATE(neig_diam(n,n));ALLOCATE(neig_MC(n,n))
	CALL Voter_ER(node_star,node_diam,neig_diam,h_in,j_in,p_in,a,b,a1,b1,pp,pp1,link)
	CALL copy_prob(node_star,node_MC,neig_diam,neig_MC)
	CALL init_voterER_MC(node_MC)
	CALL initialization_config(node_star,config_matrix)
	neig_star=neig_diam
	CALL config_print(node_star,neig_star,link,mu,z)
	CALL config_print(node_diam,neig_diam,link,mu,zz)
	CALL diamond_Voter_approx(node_diam,neig_diam,config_matrix,p_in,link,tol,tmin,tmax,u,uu,uuu)
!	CALL star_Voter_approx(node_star,neig_star,config_matrix,p_in,link,tol,tmin,tmax,v,vv,vvv)
!	CALL Magn_VoterER(node_MC,neig_MC,p_in,tol,tmin,tmax,tmax2,w,ww)
!	CALL Magn_VoterER_1sim(node_MC,neig_MC,tmax1sim,www)
	DEALLOCATE(node_star);DEALLOCATE(node_diam);DEALLOCATE(node_MC);
	DEALLOCATE(neig_star);DEALLOCATE(neig_diam);DEALLOCATE(neig_MC)

CLOSE(u);CLOSE(uu);CLOSE(uuu);
CLOSE(v);CLOSE(vv);CLOSE(vvv);
CLOSE(w);CLOSE(ww);CLOSE(www);
CLOSE(z);CLOSE(zz);

END PROGRAM VoterER
