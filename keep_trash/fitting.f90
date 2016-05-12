module fitting

implicit none

integer, parameter      :: sp = 8, lgt=kind(.true.), i4b=4, spc = 8, dpc = 16, dp = 8

    INTERFACE assert_eq
        MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
    END INTERFACE
    
    INTERFACE swap
        MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
            swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
            masked_swap_rs,masked_swap_rv,masked_swap_rm
    END INTERFACE
    
    INTERFACE diagmult
        MODULE PROCEDURE diagmult_rv,diagmult_r
    END INTERFACE
    
!     INTERFACE outerprod
!         MODULE PROCEDURE outerprod_r,outerprod_d
!     END INTERFACE

 contains

    SUBROUTINE mrqmin(x,y,sig,a,maska,covar,alpha,chisq,funcs,alamda)
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: x,y,sig
        REAL(SP), DIMENSION(:), INTENT(INOUT) :: a
        REAL(SP), DIMENSION(:,:), INTENT(OUT) :: covar,alpha
        REAL(SP), INTENT(OUT) :: chisq
        REAL(SP), INTENT(INOUT) :: alamda
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
        INTEGER(I4B) :: ma,ndata
        INTEGER(I4B), SAVE :: mfit
        INTERFACE
            SUBROUTINE funcs(x,a,yfit,dyda)
            REAL(8), DIMENSION(:), INTENT(IN) :: x,a
            REAL(8), DIMENSION(:), INTENT(OUT) :: yfit
            REAL(8), DIMENSION(:,:), INTENT(OUT) :: dyda
            END SUBROUTINE funcs
        END INTERFACE
        call mrqmin_private
        CONTAINS
    !BL
        SUBROUTINE mrqmin_private
            REAL(SP), SAVE :: ochisq
            REAL(SP), DIMENSION(:), ALLOCATABLE, SAVE :: atry,beta
            REAL(SP), DIMENSION(:,:), ALLOCATABLE, SAVE :: da
            ndata=assert_eq(size(x),size(y),size(sig),'mrqmin: ndata')
            ma=assert_eq((/size(a),size(maska),size(covar,1),size(covar,2),&
                size(alpha,1),size(alpha,2)/),'mrqmin: ma')
            mfit=count(maska)
            if (alamda < 0.d0) then
                allocate(atry(ma),beta(ma),da(ma,1))
                alamda=0.001d0
                call mrqcof(a,alpha,beta)
                ochisq=chisq
                atry=a
            end if
            covar(1:mfit,1:mfit)=alpha(1:mfit,1:mfit)
            call diagmult(covar(1:mfit,1:mfit),1.d0+alamda)
            da(1:mfit,1)=beta(1:mfit)
            call gaussj(covar(1:mfit,1:mfit),da(1:mfit,1:1))
            if (alamda == 0.d0) then
                call covsrt(covar,maska)
                deallocate(atry,beta,da)
                RETURN
            end if
            atry=a+unpack(da(1:mfit,1),maska,0.d0)
            call mrqcof(atry,covar,da(1:mfit,1))
            if (chisq < ochisq) then
                alamda=0.1d0*alamda
                ochisq=chisq
                alpha(1:mfit,1:mfit)=covar(1:mfit,1:mfit)
                beta(1:mfit)=da(1:mfit,1)
                a=atry
            else
                alamda=10.d0*alamda
                chisq=ochisq
            end if
        END SUBROUTINE mrqmin_private
    !BL
        SUBROUTINE mrqcof(a,alpha,beta)
            REAL(SP), DIMENSION(:), INTENT(IN) :: a
            REAL(SP), DIMENSION(:), INTENT(OUT) :: beta
            REAL(SP), DIMENSION(:,:), INTENT(OUT) :: alpha
            INTEGER(I4B) :: j,k,l,m
            REAL(SP), DIMENSION(size(x),size(a)) :: dyda
            REAL(SP), DIMENSION(size(x)) :: dy,sig2i,wt,ymod
            call funcs(x,a,ymod,dyda)
            sig2i=1.d0/(sig**2)
            dy=y-ymod
            j=0
            do l=1,ma
                if (maska(l)) then
                    j=j+1
                    wt=dyda(:,l)*sig2i
                    k=0
                    do m=1,l
                        if (maska(m)) then
                            k=k+1
                            alpha(j,k)=dot_product(wt,dyda(:,m))
                            alpha(k,j)=alpha(j,k)
                        end if
                    end do
                    beta(j)=dot_product(dy,wt)
                end if
            end do
            chisq=dot_product(dy**2,sig2i)
        END SUBROUTINE mrqcof
    END SUBROUTINE mrqmin
    
    SUBROUTINE covsrt(covar,maska)
        REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: covar
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: maska
        INTEGER(I4B) :: ma,mfit,j,k
        ma=assert_eq(size(covar,1),size(covar,2),size(maska),'covsrt')
        mfit=count(maska)
        covar(mfit+1:ma,1:ma)=0.d0
        covar(1:ma,mfit+1:ma)=0.d0
        k=mfit
        do j=ma,1,-1
            if (maska(j)) then
                call swap(covar(1:ma,k),covar(1:ma,j))
                call swap(covar(k,1:ma),covar(j,1:ma))
                k=k-1
            end if
        end do
    END SUBROUTINE covsrt
    
    SUBROUTINE gaussj(a,b)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
    LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
    REAL(SP) :: pivinv
    REAL(SP), DIMENSION(size(a,1)) :: dumc
    INTEGER(I4B), TARGET :: irc(2)
    INTEGER(I4B) :: i,l,n
    INTEGER(I4B), POINTER :: irow,icol
    n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n
        lpiv = (ipiv == 0)
        irc=maxloc(abs(a),outerand(lpiv,lpiv))
        ipiv(icol)=ipiv(icol)+1
        if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
        if (irow /= icol) then
            call swap(a(irow,:),a(icol,:))
            call swap(b(irow,:),b(icol,:))
        end if
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol) == 0.0) &
            call nrerror('gaussj: singular matrix (2)')
        pivinv=1.d0/a(icol,icol)
        a(icol,icol)=1.d0
        a(icol,:)=a(icol,:)*pivinv
        b(icol,:)=b(icol,:)*pivinv
        dumc=a(:,icol)
        a(:,icol)=0.d0
        a(icol,icol)=pivinv
        a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
        b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
        a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
        b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do
    do l=n,1,-1
        call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
    END SUBROUTINE gaussj
    
    
!BL
    FUNCTION outerprod(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod
    
    FUNCTION outerand(a,b)
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
        spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerand

    
    SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
    END SUBROUTINE nrerror
    
    FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
        assert_eq2=n1
    else
        write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq2'
    end if
    END FUNCTION assert_eq2
!BL
    FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
        assert_eq3=n1
    else
        write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq3'
    end if
    END FUNCTION assert_eq3
!BL
    FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
        assert_eq4=n1
    else
        write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eq4'
    end if
    END FUNCTION assert_eq4
!BL
    FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
        assert_eqn=nn(1)
    else
        write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
        STOP 'program terminated by assert_eqn'
    end if
    END FUNCTION assert_eqn
    
        SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_i
!BL
    SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_r
!BL
    SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_rv
!BL
    SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_c
!BL
    SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_cv
!BL
    SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_cm
!BL
    SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_z
!BL
    SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_zv
!BL
    SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_zm
!BL
    SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
        swp=a
        a=b
        b=swp
    end if
    END SUBROUTINE masked_swap_rs
!BL
    SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    END SUBROUTINE masked_swap_rv
!BL
    SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
        swp=a
        a=b
        b=swp
    end where
    END SUBROUTINE masked_swap_rm
    
    SUBROUTINE diagmult_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
        mat(j,j)=mat(j,j)*diag(j)
    end do
    END SUBROUTINE diagmult_rv
!BL
    SUBROUTINE diagmult_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
        mat(j,j)=mat(j,j)*diag
    end do
    END SUBROUTINE diagmult_r
end module fitting