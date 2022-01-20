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
REAL*16::x,xx,a,b
REAL*16::pp

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




SUBROUTINE init_prob1D(array,array2,matrix,a,b,a1,b1)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix
INTEGER::i,j
REAL*16::x,xx,y,yy,a,b,a1,b1,pp
DO i=1,SIZE(array)
	CALL random_number(y)
	CALL random_number(x)
	xx=a+(b-a)*x
	yy=a1+(b1-a1)*y

	array(i)%p(1,3)=yy; array(i)%p(2,3)=1-yy	!STAR initialization
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
END SUBROUTINE init_prob1D


SUBROUTINE copy_prob(array,array2,matrix,matrix2)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix,matrix2
array2=array
matrix2=matrix
END SUBROUTINE copy_prob


SUBROUTINE init_voter1D_MC(array)
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
END SUBROUTINE init_voter1D_MC
!*************************************NEIGHBOUR LIST****************************************************
SUBROUTINE init_Voter1D(array,array2,matrix,hin,jin,a,b,a1,b1,pp,pp1)
TYPE(node_list),DIMENSION(:)::array,array2
TYPE(neighbour_list),DIMENSION(:,:)::matrix
INTEGER::i,j,si,sj,n,k
REAL*16::x,y,hin,jin,a,b,a1,b1,pp
REAL*16,DIMENSION(:)::pp1
n=SIZE(array)
DO i=1,n
	array(i)%names=i;array2(i)%names=i
	array(i)%degree=2;array2(i)%degree=2	!initialization degree
	array(i)%p=0;array2(i)%p=0		!initialization 1-time-probability
	array(i)%h=hin;array2(i)%h=hin		!initialization field
	DO j=1,array(i)%degree
		DO si=1,2
		DO sj=1,2	
			matrix(i,j)%pp(sj,si,:)=0
		END DO 
		END DO 
	END DO
END DO
DO i=2,n-1
	matrix(i,1)%na=i-1	!initialization node indexes 
	matrix(i,2)%na=i+1
END DO
matrix%jj=jin			!initialization coupling constant	
matrix(1,1)%na=n; matrix(1,2)%na=2; matrix(n,1)%na=n-1; matrix(n,2)%na=1   !periodic boundary condition

CALL init_prob1D(array,array2,matrix,a,b,a1,b1)
!CALL init_probER(array,array2,matrix,a,b,pp)
!CALL init_prob_ER3(array,array2,matrix,a,b,pp1)

END SUBROUTINE init_Voter1D
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
SUBROUTINE diamond_Voter_approx(node,neig,config_matrix,tol,tmin,tmax,u,v,vv)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER,DIMENSION(:,:),ALLOCATABLE::config_matrix
INTEGER::i,j,t,u,v,vv,tmin,tmax,tot_link
REAL*16::Pi,Pi2,tol,pin
REAL*16,DIMENSION(2)::magn,ro
WRITE (u,'(a50)') 'DIAMOND APPROXIMATIONS ON VOTER 1D'
WRITE (u,'(a25,f12.8)')'Initial probability Pi:', node(5)%p(1,2)

OPEN(UNIT=50, FILE='Data1DRandom/3Prob_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=51, FILE='Data1DRandom/3Correlation_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=54, FILE='Data1DRandom/3Correlation_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=55, FILE='Data1DRandom/3Marginalization_control_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')

WRITE ( 50,'(2a6,2a20)')'i','j','P(sj^t-1,si^t-2)','P(si^t-2)*P(sj^t-1)'

t=1;

magn=(/10.0,5.0/);ro=(/10.0,5.0/)!condition just to enter in the do while
!	DO WHILE (ABS(magn(1)-magn(2))>tol .AND. ABS(ro(1)-ro(2))>tol .AND. t<tmax )	
	DO WHILE (t<tmax)		!ABS(magn-weight_magn)>tol .AND.
	WRITE (50,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (51,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (55,'(a3,i6,a50)') 't=',t, '    *************************************'
	WRITE (54,'(a3,i6,a50)') 't=',t, '    *************************************'
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
	CALL print_correlation_control(node,neig,54)
	CALL print_marginalization_control(node,neig,55)
	CALL factorization_control(node,neig,50)
	CALL print_correlation_log(node,neig,51)
	CALL magnetization(node,magn)
        CALL correlation_function_tdiff(node,neig,ro(1))
	CALL correlation_function(node,neig,ro)

PRINT*,'time',t
PRINT*,'Active link density tdiff',ro(1)
PRINT*,'Active link density',ro(2)
PRINT*,'Average Magnetization',magn(2)
	WRITE (vv,'(i10,f20.10,3f40.37)') t,REAL(t)/SIZE(node),magn(2),ro(2),ro(1)
END DO

PRINT*,'************************DIAMOND APPROXIMATIONS******************************'
PRINT*,'Final time',t
PRINT*,'Final Average Magnetization',magn(2)
PRINT*,'Final Active link density t diff',ro(1)
PRINT*,'Final Active link density',ro(2)
PRINT*,'*************************************************************************'
CLOSE(50);CLOSE(51);CLOSE(54);CLOSE(55);
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
	DO si=1,2	!spin i a t-2			
	DO k=1,2**kmax	!spin j a t-1	
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
					IF(j==kmax)THEN			
						P=P*ND(i)%pp(sii,si,2)		
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
				IF(j==kmax)THEN			
					P=P*NP(i)%p(si,3)		
				ELSE				
					P=P*NP(EP(i,j)%na)%p(sj,3)
				END IF				
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


SUBROUTINE star_Voter_approx(node,neig,config_matrix,tol,tmin,tmax,u,v,vv)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER,DIMENSION(:,:),ALLOCATABLE::config_matrix
INTEGER::i,j,t,u,v,vv,tmin,tmax
REAL*16::Pi,Pi2,tol
REAL*16,DIMENSION(2)::magn,ro
OPEN(UNIT=52, FILE='Data1DRandom/3Correlation_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=53, FILE='Data1DRandom/3Correlation_control_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
WRITE (u,'(a50)') 'STAR APPROXIMATION ON VOTER MODEL ON ER GRAPH'

	t=1
	magn=(/ 10 , 5 /);ro=(/ 10 , 5 /);
	DO WHILE (t<tmax)	!ABS(magn(2) - weight_magn)>tol .AND. 
	IF (t>tmin)THEN
	ro(1)=ro(2)
	END IF
		WRITE (52,'(a3,i6,a50)') 't=',t, '    *************************************'
		WRITE (53,'(a3,i6,a50)') 't=',t, '    *************************************'
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
		node%p(1,3)=node%p(1,2)
		node%p(2,3)=node%p(2,2)
		t=t+1
		CALL print_correlation_log(node,neig,52)	
		CALL correlation_function(node,neig,ro)
		CALL magnetization(node,magn)
PRINT*,'time',t
PRINT*,'Average Magnetization',magn(2)
PRINT*,'Active link density',ro(2)
	WRITE (v,'(i6,f12.8)') t,node(5)%p(1,2)
	WRITE (vv,'(i6,3f12.8)') t,REAL(t)/SIZE(node),magn(2),ro(2)
	END DO
PRINT*,'************************STAR APPROXIMATIONS******************************'
PRINT*,'Final time',t
PRINT*,'Final Average Magnetization',magn(2)
PRINT*,'Final Active link density',ro(2)
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

SUBROUTINE magnetization(node,magn)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i
REAL*16,DIMENSION(:)::magn
magn(2)=0
DO i=1,SIZE(node)
	magn(2)=magn(2)+(node(i)%p(1,2)-node(i)%p(2,2))
END DO
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

SUBROUTINE correlation_function(node,neig,ro)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,si,sj,ssi,ssj
REAL*16,DIMENSION(:)::ro
REAL*16::ro_temp
ro(2)=0
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
	ro(2)=ro(2)+ro_temp
	END DO
!PRINT*,'i',i,ro(2)
END DO
ro(2)=(1-ro(2)/(2*SIZE(node)))/2
!ro(2)=ro(2)/(2*SIZE(node))
END SUBROUTINE correlation_function

SUBROUTINE correlation_function_tdiff(node,neig,ro)
TYPE(node_list),DIMENSION(:),ALLOCATABLE::node
TYPE(neighbour_list),DIMENSION(:,:),ALLOCATABLE::neig
INTEGER::i,j,si,sj,ssi,ssj
REAL*16::ro
REAL*16::ro_temp
ro=0
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
!PRINT*,'j',j,ro_temp,ro(2)
	ro=ro+ro_temp
	END DO
!PRINT*,'i',i,ro(2)
END DO
ro=(1-ro/(2*SIZE(node)))/2
END SUBROUTINE correlation_function_tdiff
!********************************************************************************************************
!**************************************VOTER DYNAMICS****************************************************
!********************************************************************************************************
SUBROUTINE voter1D_step(VN,VE,M)
	INTEGER::i,j,M
	INTEGER*8::seed=0
	REAL*16::x
	TYPE(node_list),DIMENSION(:)::VN
	TYPE(neighbour_list),DIMENSION(:,:)::VE
	CALL random_number(x)
	i=1+SIZE(VN)*x!lcg(seed)
	CALL random_number(x)
	j=1+VN(i)%degree*x!lcg(seed)
	IF( VN(i)%s(1) /= VN(VE(i,j)%na)%s(1) )THEN
		VN(i)%s(1)=VN(VE(i,j)%na)%s(1)	!direct voter model
!		VN(VE(i,j)%na)%s(1)=VN(i)%s(1)	!reverse voter model
		M=M+2*VN(i)%s(1)
!PRINT*,i,j,2*REAL(VN(i)%s(1))/SIZE(VN),M
	END IF
END SUBROUTINE voter1D_step
!***************************MAGNETIZATION DATA FOR VOTER DYNAMICS**********************************
SUBROUTINE Magn_Voter1D(NG,NE,tol,tmin,tmax,tmax2,u)
	TYPE(node_list),DIMENSION(:)::NG
	TYPE(neighbour_list),DIMENSION(:,:)::NE
	REAL*16::tol,M_av,M_av2,M_temp
	INTEGER::i,j,t,t1,t2,u,tmin,tmax,tmax2,M
	WRITE (u,'(2a10,2a15)')'#sim','Exit time','     ','Final Magn','Average magn'
	M_av=10.0;M_av2=5.0;t=1;t2=0;M_temp=0
	DO WHILE (t<tmax)! ABS(M_av-M_av2)>tol .AND.			
		IF(t>tmin)THEN		
		M_av2=M_av
		END IF
		CALL init_voter1D_MC(NG)
		M=0
		DO i=1,SIZE(NG)	!compute initial magnetization
			M = M + NG(i)%s(1)
		END DO
		t1=1
PRINT*,'num sim',t,'num sim failed',t2,'initial magn',M
		DO WHILE ( ABS(M) /= SIZE(NG) .AND. t1 < tmax2 )
			CALL Voter1D_Step(NG,NE,M)
			t1 = t1 + 1	
		END DO
		IF(t1<tmax2)THEN 	
			M_temp= M_temp + REAL(M)/SIZE(NG)
			M_av= M_temp / t
PRINT*,t,M_av	
			WRITE (u,'(2i10,2f15.10)') t,t1,REAL(M)/SIZE(NG),M_av
		t=t+1
		ELSE
		t2=t2+1
		PRINT*,t2
			IF(t2>5*tmax)THEN
			PRINT*,'consensus never reached'
			EXIT
			END IF
		END IF
	END DO
END SUBROUTINE Magn_Voter1D
!******************************************************************************************************
!********************************PRINT NETWORK CONFIGURATION*******************************************
!****************************************************************************************************** 
SUBROUTINE config_print(node,neig,u)
TYPE(node_list),DIMENSION(:)::node
TYPE(neighbour_list),DIMENSION(:,:)::neig
INTEGER::i,n,u,j
n=SIZE(node)

DO i=1,n
WRiTE (u,'(a15,i6,2(a15,f8.5))')'node',node(i)%names,'P(si_t-2=+1)',node(i)%p(1,3),'P(si_t-2=-1)',node(i)%p(2,3)
WRiTE (u,'(a15,i6,2(a15,f8.5))')'node',node(i)%names,'P(si_t-1=+1)',node(i)%p(1,2),'P(si_t-1=-1)',node(i)%p(2,2)
	WRITE (u,'(a15)') 'neighbour list'
	DO j=1,node(i)%degree
	WRITE (u,'(a4,i4,a24,4f8.5)') 'node',neig(i,j)%na,'P(++),P(-+),P(+-),P(--)',neig(i,j)%pp(:,:,2)	
	END DO	
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

PROGRAM Voter1D
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
n=5;m=2 				!initialization number of nodes and neighbours
tmin=1; tmax=100000			!minimum and max number of time step
tol=1E-3				!tolerance
tmax2=5*n**4				!max number for montecarlo
a=0.0d0;b=1d0				!range for p(si_t-2=+1) at random
a1=0.0d0;b1=1d0				!range for p(si_t-2=+1) at random
pp=0.5
pp1=(/0.25d0,0.25d0,0.25d0,0.25d0/)
!***********************************PRINT NODE AND NEIGHBOUR LIST****************************************
u=30
ALLOCATE(node_star(n));ALLOCATE(node_diam(n));ALLOCATE(node_MC(n));
ALLOCATE(neig_star(n,m));ALLOCATE(neig_diam(n,m));ALLOCATE(neig_MC(n,m))
ALLOCATE(config_matrix(1,1))
CALL init_Voter1D(node_star,node_diam,neig_diam,h_in,j_in,a,b,a1,b1,pp,pp1)
CALL copy_prob(node_star,node_MC,neig_diam,neig_MC)
CALL init_voter1D_MC(node_MC)
CALL initialization_config(node_star,config_matrix)
neig_star=neig_diam

OPEN(UNIT=u, FILE='Data1DRandom/3data_Voter1D_n40_pp025.dat', STATUS='REPLACE', ACTION='WRITE')
CALL config_print(node_diam,neig_diam,u)
CLOSE(u)

!**************************DIAMOND APPROXIMATIONS ON ERDOS RENYI GRAPH***********************************
u=34;uu=35;uuu=36;v=37;vv=38;vvv=39;w=40
OPEN(UNIT=u, FILE='Data1DRandom/3Voter1D_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=uu, FILE='Data1DRandom/3Psi_Voter1D_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=uuu, FILE='Data1DRandom/3Magn_Voter1D_diam_n5.dat', STATUS='REPLACE', ACTION='WRITE')

OPEN(UNIT=v, FILE='Data1DRandom/3Voter1D_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=vv, FILE='Data1DRandom/3Psi_Voter1D_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')
OPEN(UNIT=vvv, FILE='Data1DRandom/3Magn_Voter1D_star_n5.dat', STATUS='REPLACE', ACTION='WRITE')

OPEN(UNIT=w, FILE='Data1DRandom/3Magn_Voter1D_n5.dat', STATUS='REPLACE', ACTION='WRITE')

WRITE (uuu,'(a8,4a20)') 't   ','t/N   ','Average Magn','Link Magn','Link tdiff'
WRITE (vvv,'(a8,3a12)') 't   ','t/N   ','Average Magn','Link Magn'

!Do i=1,10
!n=n+10
CALL diamond_Voter_approx(node_diam,neig_diam,config_matrix,tol,tmin,tmax,u,uu,uuu)
!END DO

!CALL star_Voter_approx(node_star,neig_star,config_matrix,tol,tmin,tmax,v,vv,vvv)


!CALL Magn_Voter1D(node_MC,neig_MC,tol,tmin,tmax,tmax2,w)



CLOSE(u)
CLOSE(uu)
CLOSE(uuu)
CLOSE(v)
CLOSE(vv)
CLOSE(vvv)

END PROGRAM Voter1D
