	allocate (Gmid(1:Nvals-2))
	
	if (Gam(3)/=1.d20) then
		intSum=lFD(0.d0,kT)/(1.d0+exp(Gam(1)))
		fmax=intSum
		do j=1,Nvals-2!integrate from Ef to Umax
			if (mode=='f') then
				Ej=j*dE
			else
				Ej=Um-j*dE
			endif
			Umax=Um-Ej
			Gj=Gamow_general(F,W-Ej,R,gamma,Umax,xm)
			Gmid(j)=Gj(1)
			if ( .not. isnan(Gj(1)) ) then
				fj=(lFD(Ej,kT))/(1.d0+exp(Gj(1)))
			endif
			intSum=intSum+fj
			if (fj>fmax) then
				fmax=fj
			endif
			if (abs(fj/fmax)<cutoff) then
				exit
			endif
		enddo
		dG=Gam(3)
	else
		j=Nvals-1
		Gj=0.d0
		dG=0.d0
		Ej=Um
		intSum=lFD(Um,kT)/2.d0
		fmax=intSum
	endif

	Nmid=j-1
	allocate(Gmidnew(1:Nmid))
	if (mode=='t') then
		Gmidnew=Gmid(Nmid:1:-1)
	else
		Gmidnew=Gmid(1:Nmid)
	endif

	if (j==(Nvals-1) .or. mode=='t') then
		Emax=10.d0*kT
		Nmax=ceiling(Emax/dE)
		allocate (Ghigh(1:Nmax))
		do j=1,Nmax!integrate above Umax
			Gj(1)=Gj(1)-dG*dE
			Ghigh(j)=Gj(1)
			Ej=Ej+dE
			fj=lFD(Ej,kT)/(1.d0+exp(Gj(1)))
			intSum=intSum+fj
			if (fj>fmax) then
				fmax=fj
			endif
			if (abs(fj/fmax)<cutoff) then
				exit
			endif
		enddo
		Nhigh=j-1
!		print *, 'Ghigh=', Ghigh
	else
		Nhigh=0
		allocate(Ghigh(1))
	endif

	if (mode=='f') then
		Emin=-6.d0/Gam(2)
	else
		Emin=0.d0
	endif
	Nmin=ceiling(abs(Emin)/dE)
	allocate (Glow(1:Nmin))
	do j=1,Nmin!integrate below Ef
		Ej=j*dE
		Umax=Um+Ej
		Gj=Gamow_general(F,W+Ej,R,gamma,Umax,xm)
		Glow(j)=Gj(1)
		fj=lFD(-Ej,kT)/(1.d0+exp(Gj(1)))
		intSum=intSum+fj
		if (fj>fmax) then
			fmax=fj
		endif
		if (abs(fj/fmax)<cutoff)  then !.or. isnan(fj))
			exit
		endif
	enddo
	Nlow=j-1

	allocate(G(Nhigh+Nlow+Nmid))
	G=[Glow(Nlow:1:-1),Gmid(1:Nmid),Ghigh(1:Nhigh)]
	deallocate(Glow,Gmid,Ghigh,Gmidnew)
	Jcur=zs*kT*intSum*dE
	
	contains
	
	pure function lFD(E,kT) result(L)
		real(dp), intent(in)::E,kT
		real(dp) :: L
		if (E>6.d0*kT) then
			L=exp(-E/kT)
		else
			L=log(1.d0+exp(-E/kT))
		endif
	end function lFD



function Notting_num(F,W,R,gamma,T,G,E0,dE) result(heat)

	real(dp), intent(in)::F,W,R,T,gamma,E0,dE !F: local field, 
	!W: work function, R: Radius of curvature, kT: boltzmann*temperature
	real(dp),intent(in) :: G(:)
	
	
	real(dp)::outsum,insum,E,kT,heat,integrand
	integer:: i,N

	outsum=0.d0
	insum=0.d0
	N=size(G)
	E=E0
	kT=kBoltz*T
	do i=1,N
		insum=insum+dE/(1.d0+exp(G(i)))
		integrand=insum*E/(1.d0+exp(E/kT))
		outsum=outsum+integrand
		print *, E, integrand
		E=E+dE
		
	enddo
	heat=outsum*dE*zs
end function Notting_num
