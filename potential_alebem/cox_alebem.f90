!===============================================================================
!	solve 2d water wave problem with one side wavemaker and oneside free by bem
!	wen-haui tsao 2021 lsu
!===============================================================================
      program water_wave
      implicit none
      integer smtyp,nnrl,mdeg
	  integer i,j,nt,nga,narc,ntim,outstep,nfield,npl,wtyp,nnode,nelm,icon,niter,nwg,e1loc,outyp,aletyp      
	  integer,allocatable::nelem(:),me(:),ns(:),buntyp(:),ln(:,:),ipiv(:)
	  real*8 dep,clen,d_out,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,width,tho,delttime,wgx(20)
	  real*8 time,dis,vel,acc,p_atm,for,ddis,dtemp,wc,e1,e2,etol
      real*8 wave_period,wave_height,psi,amp,omega,wave_k,wave_c
	  real*8,allocatable::coor(:,:),side_l(:)
	  real*8,allocatable::node(:,:),norm(:,:),jcb(:),leng(:),phi(:),pphi(:),phit(:),pphit(:)
	  real*8,allocatable::ker1(:,:),ker2(:,:),dp(:,:),dpds(:),dpdss(:),pr(:),dpdt(:),accmo(:)
	  real*8,allocatable::phit_temp(:)
	  real*8,allocatable::wt(:),rt(:),sha1(:),sha2(:),sh(:,:),awt(:),art(:)
      real*8,allocatable::g1(:,:),h1(:,:),eye(:)
	  real*8,allocatable::di(:),ve(:),ac(:)
      real*8,allocatable::dpm(:,:),dpdtm(:),cv(:,:),d2pm(:,:),d2pdtm(:)

      open(unit=1,file='1.ipt',status='old')
      open(unit=2,file='2.ipt',status='old')
      open(unit=3,file='3.ipt',status='old')
      open(unit=4,file='4.ipt',status='old')
      open(unit=5,file='io.dat')
      open(unit=6,file='s.dat')
      open(unit=7,file='wg.dat')
!      open(unit=8,file='p.dat')
!      open(unit=9,file='f.dat')
!      open(unit=10,file='e.dat')
!      open(unit=11,file='domain.dat')
        open(unit=21,file='err.dat')
        open(unit=22,file='abort.txt')
        open(unit=23,file='cfl.dat')
      open(unit=99,file='test.txt')

!---smoothing parameters
    call input_3(smtyp,nnrl,mdeg,narc,aletyp)
    allocate(awt(narc),art(narc))
    awt=0.d0
    art=0.d0
    call gauss(awt,art,narc)
    
!---topograghy and wave type
	call input_2(npl,wtyp,outyp,nwg,wgx,wave_period,wave_height,psi)
	allocate(nelem(npl),me(npl),ns(npl),buntyp(npl),coor(npl,2),side_l(npl))
	nelem=0
	ns=0
	buntyp=0
	coor=0.d0
	side_l=0.d0

!---input all kinds of parameters
	call input_1(npl,coor,nfield,nnode,nelm,nelem,me,ns,buntyp,dep,clen,nga,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,width,tho,&
	&ntim,delttime,outstep,icon,niter,etol)
	allocate(ln(nelm,2),node(nnode,2),norm(nelm,2),jcb(nelm),leng(nelm))
	allocate(phi(nnode),pphi(nnode),phit(nnode),pphit(nnode),phit_temp(nnode))
	allocate(ker1(nnode,nnode),ker2(nnode,nnode),dp(nnode,2))
	allocate(dpds(nnode),dpdss(ns(1)),dpdt(nnode),pr(nnode),accmo(nnode))
	allocate(wt(nga),rt(nga),sha1(nga),sha2(nga),sh(2,nga))
	allocate(di(ntim),ve(ntim),ac(ntim))
    allocate(g1(nnode,nnode),h1(nnode,nnode),eye(nnode),ipiv(nnode))
    allocate(dpm(ns(1),2),dpdtm(ns(1)),cv(ns(1),2),d2pm(ns(1),2),d2pdtm(ns(1)))
	ln=0
	node=0.d0
	norm=0.d0
	jcb=0.d0
	leng=0.d0
	phi=0.d0
	pphi=0.d0
	phit=0.d0
	pphit=0.d0
	ker1=0.d0
	ker2=0.d0
	dp=0.d0
	dpds=0.d0
    dpdss=0.d0
	dpdt=0.d0
	pr=0.d0
	accmo=0.d0
	phit_temp=0.d0
	wt=0.d0
	rt=0.d0
	sha1=0.d0
	sha2=0.d0
	sh=0.d0
	di=0.d0
	ve=0.d0
	ac=0.d0
    g1=0.d0
    h1=0.d0
    eye=1.d0
    ipiv=0
    dpm=0.d0
    dpdtm=0.d0
    cv=0.d0
    d2pm=0.d0
    d2pdtm=0.d0
    
!---compute mu1 and mu2 for the manfrove section
        call input_4(mu1,mu2)

!---prepare if outlet is a wall (if a wall, no iteration needed)
    if (outyp==0)then
        niter=1
    end if
        
!---shape function and mesh
    call gauss(wt,rt,nga)
    call shap(sha1,sha2,sh,nga,rt)
	call length(npl,coor,side_l)
	call mesh(npl,nnode,nelm,nelem,ln,coor,side_l,node)
    write(*,*) 'pass mesh'

!---generate waves
	select case (wtyp)
    case(1)
        call periodic(dep,wave_period,wave_height,amp)
        omega=2.d0*dacos(-1.d0)/wave_period
		do i=1,ntim
		 time=(i-1)*delttime
		 di(i)=amp*dcos(omega*time+psi)-amp !amp*sin(omega*time+psi)
		 ve(i)=-amp*omega*dsin(omega*time+psi) !amp*omega*cos(omega*time+psi)
		 ac(i)=-amp*omega**2*dcos(omega*time+psi) !-amp*omega**2*sin(omega*time+psi)
		 end do
	case(2)
		nt=omega/delttime
		do i=1,nt+1
		time=(i-1)*delttime
		di(i)=0.5d0*amp-0.5d0*amp*dcos(dacos(-1.d0)/omega*time) !str/dur*time
		ve(i)=dacos(-1.d0)/omega*0.5d0*amp*dsin(dacos(-1.d0)/omega*time) !str/dur
		ac(i)=(dacos(-1.d0)/omega)**2*0.5d0*amp*dcos(dacos(-1.d0)/omega*time) !0.d0
		end do
		di(nt+2:ntim)=di(nt+1)
    case(3)
        nt=10.d0/delttime
        call solitary(dep,wave_height,grav,1,1,nt,delttime,di,ve,ac)
        di(nt+1:ntim)=di(nt)
        ve(nt+1:ntim)=ve(nt)
        ac(nt+1:ntim)=ac(nt)

    end select  
    
do nt=1,ntim
	write(*,*) nt,'th'
     time=(nt-1)*delttime
!---wavemaker dis, vel, acc, phase lag (in radius)
	 dis=di(nt)
     vel=ve(nt)
     acc=ac(nt)
	 write(5,"(7(e15.8,1x))") time,dis,vel,acc

!---calculate wave speed for radiation condition
	d_out=node(ns(1),2)-coor(npl-2,2)
	call wave_spd(grav,omega,d_out,wc)

!---remesh location and potential of a free-surface node at current (nt-1)th time
!---the following calculation is for the next time-step)
      call remesh(npl,nnode,nelem,node,ns,time,dep,amp,nwg,wgx)
      write(*,*) 'pass remesh'
!---output boundary
      if(mod(nt,outstep)==1)then
          write(6,"(5000(e15.8,1x))") node(:,1) 
          write(6,"(5000(e15.8,1x))") node(:,2)
      end if
      
!---build kernel of zone1, zone2, zone3 and assemble them into big kernel
	   call kernel(ker1,ker2,node,norm,jcb,leng,ln,nnode,nelm,nga,sha1,sha2,sh,wt,eye)
       write(*,*) 'pass kernel'

!---calculate pressure on the boundary
	  call pressure(icon,time,tho,grav,dep,npl,nnode,ns,node,phit,dp,pr,p_atm)

!---calculate pressure and velocity in the domain
!	  call domain(npl,nga,nfield,nnode,nelm,nelem,ns,ln,node,norm,jcb,phi,pphi,phit,pphit,&
!				 &sha1,sha2,sh,wt,tho,grav,dep,p_atm,dp,pr)

!**************implicit solver for phit and pphit on the absorption side**************
do j=1,niter

	  phit_temp=phit
!---apply bc for solving free-surface pphi
       call bound(npl,nnode,nelem,ns,buntyp,outyp,phi,pphi,phit,vel,wc)
!	   call solve_lapack(npl,phi,pphi,ker1,ker2,nnode,nelem,buntyp)
!	   call solve_back(npl,phi,pphi,ker1,ker2,nnode,nelem,buntyp)
        call solve_lapack2_1(npl,phi,pphi,ker1,ker2,g1,h1,nnode,nelem,buntyp,ipiv)

!---first-order taylor series expansion (also get tangential velocity on zone boundary)
!	  call tayes1(npl,nnode,nelm,nelem,me,ns,ln,node,norm,jcb,phi,pphi,dpds,dp,phit,dpdt,&
!				 &dep,clen,grav,mu,ardzone,vel,wc)
	  call tayes1_cs(npl,nnode,nelm,nelem,me,ns,ln,node,norm,jcb,phi,pphi,dpds,dpdss,dp,phit,dpdt,&
				 &dep,clen,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,vel,wc,awt,art,narc)

!---apply bc for solving free-surface pphit
      call accbc(npl,nnode,nelm,nelem,me,ln,norm,jcb,pphi,dp,accmo)
      call boundt(npl,nnode,nelm,ns,nelem,buntyp,outyp,ln,phit,pphit,dpds,jcb,wc,acc,accmo)
!	  call solve_lapack(npl,phit,pphit,ker1,ker2,nnode,nelem,buntyp)
!	  call solve_back(npl,phit,pphit,ker1,ker2,nnode,nelem,buntyp)
      call solve_lapack2_2(npl,phit,pphit,g1,h1,nnode,nelem,buntyp,ipiv)

if (outyp==1)then
call converge(nelem(2)+1,phit_temp(ns(1)+1:ns(2)),phit(ns(1)+1:ns(2)),e1loc,e1,e2)

	if(e1<=etol.and.e2<etol)then
		write(*,*) 'pass converged',j
		write(21,*) time,j,e1,e2
		goto 205
	else if(j>=niter)then
		write(22,*) 'cg fail',time,e1loc,e1,e2
		write(*,*) 'converge fail'
		stop
	end if
end if

end do
!**************************************************************************************

205 continue

!---check the cfl number	
!	  call courant(time,delttime,nnode,nelm,ln,node,dp,jcb)
    
!**************free surface updating in el(aletyp=0) or ale(aletyp=1) frame**************    
    if (aletyp==0)then
!---second-order taylor series expansion in el frame
        call tayes2(phi,pphi,phit,pphit,dpds,dpdt,dp,npl,nnode,nelem,node,nelm,norm,jcb,ln,delttime,grav,acc)
        
    else if (aletyp==1)then
!---calculate mesh velocity
        call ale_vel(npl,nnode,nelm,nelem,ns,ln,node,norm,dp,dpdt,dpm,dpdtm,cv)
    
!---calculate mesh acceleration
        call ale_acc(delttime,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,clen,npl,nnode,nelem,ns,nelm,ln,node,norm,jcb,&
                    &phi,pphi,phit,pphit,dp,dpds,dpdss,dpdt,cv,d2pm,d2pdtm)
        write(*,*) 'pass ale'        
        
!---tse sum in ale frame
        call ale_tse_sum(delttime,npl,nnode,ns,node,phi,dpm,dpdtm,d2pm,d2pdtm)
    end if
    
write(*,*) 'pass tayes'

    
    
    
    
!---smoothing for free-surface node and phi
      if (smtyp==1)then
          call fs_smooth(nnrl,mdeg,ns(1),node(1:ns(1),2))
          call fs_smooth(nnrl,mdeg,ns(1),phi(1:ns(1)))
          write(*,*) 'pass smooth'
      end if
end do

stop
    end


!**********************************************************************
      subroutine headline(id,iread)
!**********************************************************************
      character*2 id
      id =  '*'
      do while (id .eq. '*')
      read(iread,'(a1)') id
      end do
      return
      end
!**********************************************************************
subroutine input_1(npl,coor,nfield,nnode,nelm,nelem,me,ns,buntyp,dep,clen,nga,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,width,tho,&
				&ntim,delttime,outstep,icon,niter,etol)
!**********************************************************************
      implicit none
	  integer i,j,npl,nfield,nnode,nelm,nga,ntim,icon,niter,outstep
	  integer nelem(npl),me(npl),ns(npl),buntyp(npl)
	  real*8 coor(npl,2),dep,clen,endtime,delttime,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,width,tho,etol
      character*2 id
         id = '*'
!---nodal coordinates
         call headline(id,1)
         read(1,*)  ((coor(i,j),j=1,2),i=1,npl)
!---element mesh number
         call headline(id,1)
        read(1,*) (nelem(i),i=1,npl)
        nnode = 0
        nelm  = 0
		me(1) = nelem(1)
		ns	  = 0
        do i=1,npl
         nelm = nelm+nelem(i)
         nnode = nnode+(nelem(i)+1)
        end do

		do i=2,npl
		 me(i)=me(i-1)+nelem(i)
		end do

      ns(1)=nelem(1)+1
	  do i=2,npl
      ns(i)=ns(i-1)+nelem(i)+1
	  end do

	  nfield=(nelem(1)-1)*(nelem(npl)-1)

!---boundary type
         call headline(id,1)
        read(1,*) (buntyp(i),i=1,npl)
!---read the gaussian integration point no.
         call headline(id,1)
         read(1,*)  nga
!---read grav acc, max of wave damping coef, zone of wave damping, tank width, fluid density
         call headline(id,1)
         read(1,*)  grav,mu1,mu2,porous_x1,porous_x2,porous_x3,width,tho
!---read time,time interval,output step,index of bernoulli constant, number if iteration,tolerance error
         call headline(id,1)
         read(1,*)  endtime,delttime,outstep,icon,niter,etol
		ntim=endtime/delttime+1
		dep=maxval(coor(:,2))
        clen=maxval(coor(:,1))

      return
      end
!**********************************************************************
      subroutine input_2(npl,wtyp,outyp,nwg,wgx,wave_period,wave_height,psi)
!**********************************************************************
      implicit none
	  integer npl,wtyp,outyp,nwg
	  real*8 wave_period,wave_height,psi,wgx(20)
      character*2 id
         id = '*'
!---number of planes
         call headline(id,2)
         read(2,*) npl
!---wave generation: 1=periodic; 2=solitary
         call headline(id,2)
         read(2,*) wtyp
!---read wave_period,wave_height
         call headline(id,2)
         read(2,*)  wave_period,wave_height,psi
!---wave generation: 0=wall; 1=radiation
         call headline(id,2)
         read(2,*) outyp
!---wave gauge number
         call headline(id,2)
         read(2,*) nwg
		 if (nwg>20)then
		 write(*,*) "need more allocation for wave gauge"
		 write(22,*) "need more allocation for wave gauge"
		 stop
		 end if
!---wave gauge location
         call headline(id,2)
         read(2,*) wgx(1:nwg)

      return
    end
!**********************************************************************
      subroutine input_3(smtyp,nnrl,mdeg,narc,aletyp)
!**********************************************************************
      implicit none
	  integer smtyp,nnrl,mdeg,narc,aletyp
      character*2 id
         id = '*'
!---do you need free-surface smoothing: 1 = yes; 0 = no
         call headline(id,3)
         read(3,*) smtyp
!---number of neighboring node and degree of polynomial for sg filter
!---note that nl = nr and nl + nr < m
         call headline(id,3)
         read(3,*) nnrl,mdeg

!---number of gaussian quadrature for calculating arc length
         call headline(id,3)
         read(3,*) narc
         
!---do you use ale approach? 0 = no; 1 = yes
         call headline(id,3)
         read(3,*) aletyp
         
      return
    end
!**********************************************************************
      subroutine input_4(mu1,mu2)
!**********************************************************************
      implicit none
      integer i,nlayer
      real*8 mu1,mu2,a(28),b(28)
      character*2 id
      id = '*'
      nlayer = 24
!---read alpha
         call headline(id,4)
         read(4,*) a
!---read beta
         call headline(id,4)
         read(4,*) b

         do i=1,nlayer
         mu1=mu1+a(i)
         mu2=mu2+b(i)
         end do
         mu1=mu1/nlayer
         mu2=mu2/nlayer

      return
    end
!**********************************************************************
      subroutine length(npl,coor1,side_l1)
!**********************************************************************
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
	  integer npl
      real*8  coor1(npl,2),side_l1(npl)
      do i=1,npl-1
        side_l1(i)=dsqrt((coor1(i+1,1)-coor1(i,1))**2+(coor1(i+1,2)-coor1(i,2))**2)
      end do
        side_l1(npl)=dsqrt((coor1(npl,1)-coor1(1,1))**2+(coor1(npl,2)-coor1(1,2))**2)

      return
      end
!**********************************************************************
      subroutine mesh(npl,nnode,nelm,nelem,ln,coor,leng,node)
!********************************************************************
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      integer npl,nel(npl),nnode,nelm,nelem(npl),ln(nelm,2)
      real*8 sx,sy,norm,delt,leng(npl),coor(npl,2),node(nnode,2)

k=0
do i=1,npl-1
	j=npl-i
    nel(i) = nelem(i)+1
    delt=leng(j)/nelem(i)
	sx=coor(j,1)-coor(j+1,1)
	sy=coor(j,2)-coor(j+1,2)
	norm=dsqrt(sx**2+sy**2)
	sx=sx/norm
	sy=sy/norm
    do l=1,nelem(i)+1
       node(l+k,1)=coor(j+1,1)+(l-1)*delt*sx
       node(l+k,2)=coor(j+1,2)+(l-1)*delt*sy
    end do
    k=k+nel(i)
end do

    nel(npl) = nelem(npl)+1
    delt=leng(npl)/nelem(npl)
	sx=coor(npl,1)-coor(1,1)
	sy=coor(npl,2)-coor(1,2)
	norm=dsqrt(sx**2+sy**2)
	sx=sx/norm
	sy=sy/norm
    do i=1,nelem(npl)+1
      node(i+k,1)=coor(1,1)+(i-1)*delt*sx
      node(i+k,2)=coor(1,2)+(i-1)*delt*sy
    end do

!----to give the local element node number
      l=1
	  n=1
      do i=1,npl
       do j=1,nelem(i)
        ln(n,1)=l
        ln(n,2)=l+1
        l=l+1
		n=n+1
       end do
       l=l+1
      end do

      return
      end
!********************************************************************
      subroutine shap(sha1,sha2,sh,nga,rt)
!********************************************************************
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      integer nga
      real*8 rt(nga),sha1(nga),sha2(nga),sh(2,nga)
      do m=1,nga
        sha1(m)=0.5d0*(1-rt(m))
        sha2(m)=0.5d0*(1+rt(m))

        sh(1,m)=sha1(m)
        sh(2,m)=sha2(m)
      end do
	return
	end 
!**********************************************************************
      subroutine remesh(npl,nnode,nelem,node,ns,time,dep,amp,nwg,wgx)
!**********************************************************************
      implicit integer(i-n)
      implicit real*8(a-h,o-z)
      integer npl,nnode,nelem(npl),ns(npl),nwg,il,ir
      real*8 time,dep,amp,leng,node(nnode,2),sx,sy,norm,wgx(20),wgy(20)
      
!------ensure duplicate point on end node of free surface-----
    node(nnode,:)=node(1,:)
	node(ns(1)+1,:)=node(ns(1),:)

!------bottom end node goes with free surface -----
	node(ns(npl-1),1)=node(1,1)
	node(ns(2),1)=node(ns(1),1)

!----- remesh for all planes except the free surface
do i=2,npl
	sx=node(ns(i),1)-node(ns(i-1),1)
	sy=node(ns(i),2)-node(ns(i-1),2)
	norm=dsqrt(sx**2+sy**2)
	sx=sx/norm
	sy=sy/norm
	delt=norm/nelem(i)
	  k=1
      do l=ns(i-1)+1,ns(i)
        node(l,1)=node(ns(i-1),1)+delt*(k-1)*sx
        node(l,2)=node(ns(i-1),2)+delt*(k-1)*sy
		k=k+1
      end do
end do

!==========output wave elevation at the wave gauges (in cm)==========
do i=1,nwg
call bwloc(-wgx(i),ns(1),-node(1:ns(1),1),0,ir,il)
temp=(wgx(i)-node(il,1))/(node(ir,1)-node(il,1))
wgy(i)=node(il,2)+temp*(node(ir,2)-node(il,2))-dep
end do
write(7,"(20(e15.8,1x))") time,wgy(1:nwg)*100.d0

      return
    end
!********************************************************************
      subroutine bound(npl,nnode,nelem,ns,buntyp,outyp,phi,pphi,phit,vel,wc)
!********************************************************************
       implicit integer (i-n)
       implicit real*8 (a-h,o-z)
       integer npl,nnode,nelem(npl),ns(npl),buntyp(npl),outyp
	   real*8 r,vel,wc
       real*8 phi(nnode),pphi(nnode),phit(nnode)

       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
			   phi(j)=phi(j)
			   pphi(j)=0.d0
            else
			  if (i==2)then
                  if(outyp==0)then
                      pphi(j)=0.d0
                      phi(j)=0.d0
                  else
                      pphi(j)=-phit(j)/wc
                      phi(j)=0.d0
                  end if
			  else if (i==npl)then
			   pphi(j)=-vel
			   phi(j)=0.d0
			  else
			   pphi(j)=0.d0
			   phi(j)=0.d0
			  end if
            end if
          end do
          n=n+(nelem(i)+1)
       end do

      return
      end
!********************************************************************
      subroutine kernel(ker1,ker2,node,norm,jcb,leng,ln,nnode,nelm,nga,sha1,sha2,sh,wt,eye)
!********************************************************************
      implicit none
      integer i,j,k,l,m,n,nnode,nelm,nga,ln(nelm,2)
      real*8 rd,sigma(nnode),eye(nnode)
      real*8  ker1(nnode,nnode),ker2(nnode,nnode)
      real*8  nx,ny,h(2),g(2),xfunc(10),yfunc(10),pxi1(2)
      real*8  leng(nelm),norm(nelm,2),jcb(nelm),node(nnode,2)
      real*8  wt(nga),sha1(nga),sha2(nga),sh(2,nga)

        ker1=0.d0
        ker2=0.d0
!**** calculate the jacobian
      do j=1,nelm
      leng(j)=dsqrt((node(ln(j,1),1)-node(ln(j,2),1))**2+(node(ln(j,1),2)-node(ln(j,2),2))**2)
		do l=1,2
          pxi1(l)=(-0.5d0)*node(ln(j,1),l)+ 0.5d0*node(ln(j,2),l)
        end do
        nx=pxi1(2)
        ny=pxi1(1)
        jcb(j)=dsqrt(nx**2+ny**2)
        norm(j,1)=-nx/jcb(j)
        norm(j,2)=ny/jcb(j)
      end do
      
!$omp parallel do private(i,j,k,m,xfunc,yfunc,rd,g,h)
!***the surface kernels
      do i = 1,nnode
       do j=1,nelm
       do m=1,nga
          xfunc(m)=sha1(m)*node(ln(j,1),1)+sha2(m)*node(ln(j,2),1)
          yfunc(m)=sha1(m)*node(ln(j,1),2)+sha2(m)*node(ln(j,2),2)
       end do
        do k=1,2
         g(k)=0.d0
         h(k)=0.d0
         rd=dsqrt((node(i,1)-node(ln(j,k),1))**2+(node(i,2)-node(ln(j,k),2))**2)
        if (rd .le. 0.0000001d0) then
        h(k)=0.d0
	    g(k)=leng(j)/2*(1.5d0-dlog(leng(j)))
		else !---non dsinguler term
         g(k)=0.d0
         do m=1,nga
            h(k)=h(k)+(-1.d0)/((xfunc(m)-node(i,1))**2+(yfunc(m)-node(i,2))**2)*&
                    &((xfunc(m)-node(i,1))*norm(j,1)+(yfunc(m)-node(i,2))*norm(j,2))&
                    &*jcb(j)*sh(k,m)*wt(m)
            g(k)=g(k)+dlog(1.d0/((xfunc(m)-node(i,1))**2+(yfunc(m)-node(i,2))**2)**0.5d0)*jcb(j)*sh(k,m)*wt(m)
         end do
      end if
         ker1(i,ln(j,k))=ker1(i,ln(j,k))+h(k)
         ker2(i,ln(j,k))=ker2(i,ln(j,k))+g(k)
      end do
      end do
      end do
!$omp end parallel do
      
!***dsingular of ker1
    do i=1,nnode
    ker1(i,i)=0.d0
    end do      
      call dgemm('n','n',nnode,1,nnode,1.d0,ker1,nnode,eye,nnode,0.d0,sigma,nnode)
    do i=1,nnode
    ker1(i,i)=-sigma(i)
    end do  
!       do i=1,nnode
!          sigma=0.d0
!          do j=1,nnode
!             sigma=ker1(i,j)+sigma
!          end do
!          ker1(i,i)=-sigma
!       end do

      return
      end
!**********************************************************************
       subroutine solve_back(npl,phi,pphi,ker1,ker2,nnode,nelem,buntyp)
!**********************************************************************
!      to solve ker1*phi=ker2*pphi
!      pphi=partial phi over partial n
!      using gaussian elimination with backsubstitution
!======================================================================
       implicit integer (i-n)
       implicit real*8    (a-h,o-z)
       integer npl,nnode,nelem(npl),buntyp(npl)
       real*8    phi(nnode),pphi(nnode)
       real*8    ker1(nnode,nnode),ker2(nnode,nnode)
       real*8    h1(nnode,nnode),q1(nnode),temp(nnode)
       real*8  sum,a,sig,g1(nnode,nnode),p1(nnode)

!********** to move the ker1 and ker2**********************************
!-----phi pphi----         
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               q1(j)=phi(j)
            else
               q1(j)=pphi(j)
            end if
          end do
          n=n+(nelem(i)+1)
       end do
!-------------------
         do i=1,nnode
           n=0
           do l=1,npl
             do j=k+n,(nelem(l)+1)+n
             if (buntyp(l) .eq. 1) then
               h1(i,j)=ker1(i,j)  !g*unknown=h*known
               g1(i,j)=ker2(i,j)
              else
               h1(i,j)=-ker2(i,j)
               g1(i,j)=-ker1(i,j)
             end if
             end do
            n=n+(nelem(l)+1)
          end do
        end do

       do i=1,nnode
          temp(i)=0.d0
          do j=1,nnode
          temp(i)=temp(i)+h1(i,j)*q1(j)
          end do
       end do
!*************gaussian elimination with backsubstitution*********
      do i=1,nnode
        g1(i,nnode+1)=temp(i)
      end do
      do i=1,nnode
         sum=0.d0
         do k=i,nnode
            if (g1(k,i) .ne. 0) then
               if (k .ne. i) then
               if (g1(i,i) .eq. 0.d0) then
                 write(*,*) 'some of the diag-terms are zero'
                 write(22,*) 'some of the diag-terms are zero'
                 stop
               end if
               a=g1(k,i)/g1(i,i)
               do j=i,nnode+1
                  g1(k,j)=g1(k,j)-a*g1(i,j)
               end do
               end if
            end if
            sum=sum+g1(k,i)
          end do
          if (sum .eq. 0.d0) then
          write(*,*) 'no unique solution exists  stop at gausseli'
          write(22,*) 'no unique solution exists  stop at gausseli'
          stop
          end if
      end do

      p1(nnode)=g1(nnode,nnode+1)/g1(nnode,nnode)
      do i=nnode-1,1,-1
         sig=0.d0
         do j=i+1,nnode
            sig=g1(i,j)*p1(j)+sig
          end do
         p1(i)=(g1(i,nnode+1)-sig)/g1(i,i)
      end do
!=================================
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               pphi(j)=p1(j)
            else
               phi(j)=p1(j)
            end if
          end do
          n=n+nelem(i)+1
       end do

      return
      end
!**********************************************************************
       subroutine solve_lapack(npl,phi,pphi,ker1,ker2,nnode,nelem,buntyp)
!**********************************************************************
!      to solve ker1*phi=ker2*pphi
!      pphi=partial phi over partial n
!======================================================================
       implicit integer (i-n)
       implicit real*8    (a-h,o-z)
       integer npl,nnode,nelem(npl),buntyp(npl)
	   integer info,ipiv(nnode)
       real*8    phi(nnode),pphi(nnode)
       real*8    ker1(nnode,nnode),ker2(nnode,nnode)
       real*8    h1(nnode,nnode),q1(nnode),g1(nnode,nnode),p1(nnode)

!********** to move the ker1 and ker2**********************************
!-----phi pphi----         
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               q1(j)=phi(j)
            else
               q1(j)=pphi(j)
            end if
          end do
          n=n+(nelem(i)+1)
       end do
!-------------------
         do i=1,nnode
           n=0
           do l=1,npl
             do j=k+n,(nelem(l)+1)+n
             if (buntyp(l) .eq. 1) then
               h1(i,j)=ker1(i,j)  !g*unknown=h*known
               g1(i,j)=ker2(i,j)
              else
               h1(i,j)=-ker2(i,j)
               g1(i,j)=-ker1(i,j)
             end if
             end do
            n=n+(nelem(l)+1)
          end do
        end do

       do i=1,nnode
          p1(i)=0.d0
          do j=1,nnode
          p1(i)=p1(i)+h1(i,j)*q1(j)
          end do
       end do

!*************solve by calling lapack*********
call dgesv(nnode,1,g1,nnode,ipiv,p1,nnode,info)

!=================================
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               pphi(j)=p1(j)
            else
               phi(j)=p1(j)
            end if
          end do
          n=n+nelem(i)+1
       end do

      return
    end
!**********************************************************************
       subroutine solve_lapack2_1(npl,phi,pphi,ker1,ker2,g1,h1,nnode,nelem,buntyp,ipiv)
!**********************************************************************
       implicit none
       integer i,j,k,l,n,npl,nnode,nelem(npl),buntyp(npl)
	   integer info,ipiv(nnode)
       real*8 phi(nnode),pphi(nnode),ker1(nnode,nnode),ker2(nnode,nnode)
       real*8 h1(nnode,nnode),q1(nnode),g1(nnode,nnode),p1(nnode)
	   character*1 trans
       trans = 'n'
       
!-----move phi and pphi----         
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               q1(j)=phi(j)
            else
               q1(j)=pphi(j)
            end if
          end do
          n=n+(nelem(i)+1)
       end do
       
!-----move ker1 and ker2---- 
         do i=1,nnode
           n=0
           do l=1,npl
             do j=k+n,(nelem(l)+1)+n
             if (buntyp(l) .eq. 1) then
               h1(i,j)=ker1(i,j)  !g*unknown=h*known
               g1(i,j)=ker2(i,j)
              else
               h1(i,j)=-ker2(i,j)
               g1(i,j)=-ker1(i,j)
             end if
             end do
            n=n+(nelem(l)+1)
          end do
        end do

       p1=matmul(h1,q1)
       
!*************solve by calling lapack*********
call dgetrf(nnode,nnode,g1,nnode,ipiv,info)
call dgetrs(trans,nnode,1,g1,nnode,ipiv,p1,nnode,info)
       
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               pphi(j)=p1(j)
            else
               phi(j)=p1(j)
            end if
          end do
          n=n+nelem(i)+1
       end do
                     
      return
    end
!**********************************************************************
       subroutine solve_lapack2_2(npl,phi,pphi,g1,h1,nnode,nelem,buntyp,ipiv)
!**********************************************************************
       implicit none
       integer i,j,k,l,n,npl,nnode,nelem(npl),buntyp(npl)
	   integer info,ipiv(nnode)
       real*8 phi(nnode),pphi(nnode)
       real*8 h1(nnode,nnode),q1(nnode),g1(nnode,nnode),p1(nnode)
	   character*1 trans
       trans = 'n'
       
!-----move phi and pphi----         
       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               q1(j)=phi(j)
            else
               q1(j)=pphi(j)
            end if
          end do
          n=n+(nelem(i)+1)
       end do

       p1=matmul(h1,q1)

!*************solve by calling lapack*********
call dgetrs(trans,nnode,1,g1,nnode,ipiv,p1,nnode,info)

       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
               pphi(j)=p1(j)
            else
               phi(j)=p1(j)
            end if
          end do
          n=n+nelem(i)+1
       end do

      return
      end
!********************************************************************
subroutine tayes1(npl,nnode,nelm,nelem,me,ns,ln,node,norm,jcb,phi,pphi,&
	  &dpds,dp,phit,dpdt,dep,clen,grav,mu,ardzone,vel,wc)
!********************************************************************
    implicit none
    integer i,j,k,npl,nnode,nelm,nelem(npl),me(npl),ns(npl),ln(nelm,2)
    real*8 wc,dep,clen,grav,mu,ardzone,vel
	real*8 node(nnode,2),norm(nelm,2),jcb(nelm),phi(nnode),pphi(nnode),dp(nnode,2)
	real*8 dpds(nnode),phit(nnode),dpdt(nnode)

	dpds=0.d0
!*********************on free surface*********************
    do i=1,me(1)
    do j=1,2
		if(ln(i,j).eq.1) then
            dp(ln(i,j),1)=-pphi(nnode)
            dpds(ln(i,j))=(dp(ln(i,j),1)-pphi(ln(i,j))*norm(i,1))/norm(i,2)
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(1)) then
		    dp(ln(i,j),1)=pphi(ns(1)+1)
            dpds(ln(i,j))=(dp(ln(i,j),1)-pphi(ln(i,j))*norm(i,1))/norm(i,2)
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
        else
			dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do
        
    do i=1,ns(1)
        if (node(i,1)>=ardzone)then
          dpdt(i)=0.5d0*(dp(i,1)**2+dp(i,2)**2)-grav*(node(i,2)-dep)-mu*(node(i,1)-ardzone)/(clen-ardzone)*phi(i)
        else
          dpdt(i)=0.5d0*(dp(i,1)**2+dp(i,2)**2)-grav*(node(i,2)-dep)
        end if
        phit(i)=dpdt(i)-(dp(i,1)**2+dp(i,2)**2)
    end do

!*********************on right wall*********************
    do i=me(1)+1,me(2)
    do j=1,2
		if(ln(i,j).eq.ns(1)+1) then
			dpds(ln(i,j))=-dp(ns(1),2)
            dp(ln(i,j),1)=pphi(ln(i,j))
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(2)) then
			dpds(ln(i,j))=0.d0
            dp(ln(i,j),1)=pphi(ln(i,j))
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
			dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
			dp(ln(i,j),1)=pphi(ln(i,j))
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

!*********************left wall*********************
    do i=me(npl-1)+1,me(npl)
    do j=1,2
		if(ln(i,j).eq.ns(npl-1)+1) then
			dpds(ln(i,j))=0.d0
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.nnode) then
			dpds(ln(i,j))=dp(1,2)
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
			dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

!*********************bottom*********************
    do i=me(2)+1,me(npl-1)
    do j=1,2
		if(ln(i,j).eq.ns(2)+1) then
          dpds(ln(i,j))=-pphi(ns(2))
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
          dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(npl-1)) then
          dpds(ln(i,j))=pphi(ns(npl-1)+1)
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
          dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
		  dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
		  dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

      return
    end
!********************************************************************
subroutine tayes1_cs(npl,nnode,nelm,nelem,me,ns,ln,node,norm,jcb,phi,pphi,&
	  &dpds,dpdss,dp,phit,dpdt,dep,clen,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,vel,wc,awt,art,narc)
!********************************************************************
    implicit none
    integer i,j,k,narc,npl,nnode,nelm,nelem(npl),me(npl),ns(npl),ln(nelm,2)
    real*8 wc,dep,clen,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,vel,v_sqr
	real*8 node(nnode,2),norm(nelm,2),jcb(nelm),phi(nnode),pphi(nnode),dp(nnode,2)
	real*8 dpds(nnode),phit(nnode),dpdt(nnode),awt(narc),art(narc),dpdss(ns(1))

	dpds=0.d0
!*********************on free surface*********************
    call fs_csdiff(ns(1),node(1:ns(1),1),node(1:ns(1),2),phi(1:ns(1)),awt,art,narc,dpds(1:ns(1)),dpdss)
    
    do i=1,me(1)
    do j=1,2
		if(ln(i,j).eq.1) then
            dp(ln(i,j),1)=-pphi(nnode)
            dpds(ln(i,j))=(dp(ln(i,j),1)-pphi(ln(i,j))*norm(i,1))/norm(i,2)
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(1)) then
		    dp(ln(i,j),1)=pphi(ns(1)+1)
            dpds(ln(i,j))=(dp(ln(i,j),1)-pphi(ln(i,j))*norm(i,1))/norm(i,2)
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
        else
			!dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i) already got from cubic spline interpolation
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do
    
    do i=1,ns(1)
    v_sqr = dp(i,1)**2+dp(i,2)**2
        if (node(i,1)>=porous_x1 .and. node(i,1)<=porous_x2)then
           dpdt(i)=0.5d0*v_sqr-grav*(node(i,2)-dep)-mu1*phi(i)-mu2*phi(i)*dsqrt(v_sqr)
        else if (node(i,1)>=porous_x3)then
           dpdt(i)=0.5d0*v_sqr-grav*(node(i,2)-dep)-mu1*phi(i)
!          dpdt(i)=0.5d0*v_sqr-grav*(node(i,2)-dep)-mu1*(node(i,1)-ardzone)/(clen-ardzone)*phi(i)
        else
          dpdt(i)=0.5d0*v_sqr-grav*(node(i,2)-dep)
        end if
        phit(i)=dpdt(i)-v_sqr
    end do

!*********************on right wall*********************
    do i=me(1)+1,me(2)
    do j=1,2
		if(ln(i,j).eq.ns(1)+1) then
			dpds(ln(i,j))=-dp(ns(1),2)
            dp(ln(i,j),1)=pphi(ln(i,j))
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(2)) then
			dpds(ln(i,j))=0.d0
            dp(ln(i,j),1)=pphi(ln(i,j))
            dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
			dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
			dp(ln(i,j),1)=pphi(ln(i,j))
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

!*********************left wall*********************
    do i=me(npl-1)+1,me(npl)
    do j=1,2
		if(ln(i,j).eq.ns(npl-1)+1) then
			dpds(ln(i,j))=0.d0
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.nnode) then
			dpds(ln(i,j))=dp(1,2)
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
			dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
			dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
			dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

!*********************bottom*********************
    do i=me(2)+1,me(npl-1)
    do j=1,2
		if(ln(i,j).eq.ns(2)+1) then
          dpds(ln(i,j))=-pphi(ns(2))
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
          dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else if(ln(i,j).eq.ns(npl-1)) then
          dpds(ln(i,j))=pphi(ns(npl-1)+1)
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
          dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		else
		  dpds(ln(i,j))=(-0.5d0*phi(ln(i,1))+0.5d0*phi(ln(i,2)))/jcb(i)
		  dp(ln(i,j),1)=dpds(ln(i,j))*norm(i,2)+pphi(ln(i,j))*norm(i,1)
		  dp(ln(i,j),2)=-dpds(ln(i,j))*norm(i,1)+pphi(ln(i,j))*norm(i,2)
		end if
    end do
    end do

      return
    end
!**********************************************************************
      subroutine fs_csdiff(n,x,y,phi,awt,art,narc,dpds,dpdss)
!**********************************************************************
      implicit none
      integer i,j,n,ne,narc
      real*8 x(n),y(n),phi(n),dpds(n),dpdss(n),awt(narc),art(narc),s(n),darc,u
      real*8 b(n),c(n),d(n),y0,y1,y2

!---calculate the arc length and tangential coordinate
      call spline(n,x,y,b,c,d)
      s(1)=0.d0
      do i=2,n
          darc=0.d0
          do j=1,narc
              u=0.5d0*(1-art(j))*x(i-1)+0.5d0*(1+art(j))*x(i)
              call seval(n,u,x,y,b,c,d,y0,y1,y2)
              darc=darc+dsqrt(1.d0+y1**2)*0.5d0*(x(i)-x(i-1))*awt(j)
          end do
          s(i)=s(i-1)+darc      
      end do

!---calculate dphi/ds
      call spline(n,s,phi,b,c,d)
      do i=1,n
          call seval(n,s(i),s,phi,b,c,d,y0,dpds(i),dpdss(i))
      end do

      return
    end
!********************************************************************
      subroutine accbc(npl,nnode,nelm,nelem,me,ln,norm,jcb,pphi,dp,accmo)
!********************************************************************
    implicit none
    integer i,j,k,npl,nnode,nelm,nelem(npl),me(npl),ln(nelm,2)
    real*8 dpndx,dpndy,dpnds
	real*8 norm(nelm,2),jcb(nelm),pphi(nnode),dp(nnode,2),accmo(nnode)

do k=1,npl-1
  do i=me(k)+1,me(k+1)
    do j=1,2
		dpnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
		dpndx=dpnds*norm(i,2)	! phinn=phiss=0 for linear element
		dpndy=-dpnds*norm(i,1)
		accmo(ln(i,j))=dp(ln(i,j),1)*dpndx+dp(ln(i,j),2)*dpndy
    end do
  end do
end do

      return
      end
!********************************************************************
      subroutine boundt(npl,nnode,nelm,ns,nelem,buntyp,outyp,ln,phit,pphit,dpds,jcb,wc,acc,accmo)
!********************************************************************
       implicit integer (i-n)
       implicit real*8    (a-h,o-z)
       integer nnode,nelm,ns(npl),nelem(npl),buntyp(npl),ln(nelm,2),outyp
	   real*8 r,acc,wc
       real*8 jcb(nelm),phit(nnode),pphit(nnode),dpds(nnode),accmo(nnode),dpdss(nnode)

	pphit=0.d0

    do i=nelem(1)+1,nelem(1)+nelem(2)
    do j=1,2
		if(ln(i,j).eq.ns(2)) then
		  dpdss(ln(i,j))=0.d0 ! assume a flat ground so vn=0
		else
		  dpdss(ln(i,j))=0.5d0*(-dpds(ln(i,1))+dpds(ln(i,2)))/jcb(i)
		end if
    end do
    end do

       k=1
       n=0
       do i=1,npl
          do j=k+n,(nelem(i)+1)+n
            if (buntyp(i) .eq. 1) then
			 phit(j)=phit(j)
			 pphit(j)=0.d0
            else
			  if (i==2)then
                  if (outyp==0)then
                      pphit(j)=0.d0
                      phit(j)=0.d0 
                  else
                      pphit(j)=wc*dpdss(j)
                      phit(j)=0.d0
                  end if
			  else if (i==npl)then
			    pphit(j)=-acc-accmo(j)
			    phit(j)=0.d0
			  else
			    pphit(j)=0.d0-accmo(j)
			    phit(j)=0.d0
			  end if
            end if
          end do
          n=n+(nelem(i)+1)
       end do

      return
      end
!**********************************************************************
subroutine tayes2(phi,pphi,phit,pphit,dpds,dpdt,dp,npl,nnode,nelem,node,nelm,norm,&
				 &jcb,ln,delttime,grav,acc)
!**********************************************************************
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      integer npl,nnode,nelm,nelem(npl),ln(nelm,2)
      real*8 delttime,grav,acc,d2pdt,dptds,dpdnds
      real*8 norm(nelm,2),jcb(nelm)
      real*8 node(nnode,2),phi(nnode),pphi(nnode),phit(nnode),pphit(nnode)
	  real*8 dp(nnode,2),dpds(nnode),dpdt(nnode)
      real*8 new_node(nelem(1)+1,2),new_phi(nelem(1)+1),d2p(2)

	do i=1,nelem(1)
	  do j=1,2
		if(ln(i,j) .eq. 1) then
				d2p(1)=-pphit(nnode)
				dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
				dptds=(d2p(1)-pphi(1)*dpdnds*norm(i,2)-(dpds(1)*dpdnds+pphit(1))*norm(1,1))/norm(1,2)
				d2p(2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
		else if (ln(i,j) .eq. nelem(1)+1) then
				d2p(1)=pphit(nelem(1)+2)
				dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
				dptds=(d2p(1)-pphi(1)*dpdnds*norm(i,2)-(dpds(1)*dpdnds+pphit(1))*norm(1,1))/norm(1,2)
				d2p(2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
		else
			dptds=(-0.5d0*phit(ln(i,1))+0.5d0*phit(ln(i,2)))/jcb(i)
			dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
			d2p(1)=(dptds+pphi(ln(i,j))*dpdnds)*norm(i,2)-(-dpds(ln(i,j))*dpdnds-pphit(ln(i,j)))*norm(i,1)
			d2p(2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
		end if

		d2pdt=dp(ln(i,j),1)*d2p(1)+dp(ln(i,j),2)*d2p(2)-grav*dp(ln(i,j),2)
		new_phi(ln(i,j))=phi(ln(i,j))+delttime*dpdt(ln(i,j))+d2pdt*delttime**2/2
		new_node(ln(i,j),1)=node(ln(i,j),1)+dp(ln(i,j),1)*delttime+0.5d0*d2p(1)*delttime**2
		new_node(ln(i,j),2)=node(ln(i,j),2)+dp(ln(i,j),2)*delttime+0.5d0*d2p(2)*delttime**2
	  end do
	end do

	do i=1,nelem(1)+1
	node(i,:)=new_node(i,:)
	phi(i)=new_phi(i)
	end do

      return
    end
!**********************************************************************
      subroutine ale_vel(npl,nnode,nelm,nelem,ns,ln,node,norm,dp,dpdt,dpm,dpdtm,cv)
!**********************************************************************
      implicit none
      integer i,j,npl,nnode,nelm,nelem(npl),ns(npl),ln(nelm,2)
      real*8 node(nnode,2),norm(nelm,2),dp(nnode,2),dpdt(nnode),dpm(ns(1),2),dpdtm(ns(1)),cv(ns(1),2)
      
!---calculate mesh velocity      
      do i=1,nelem(1)
          do j=1,2
              dpm(ln(i,j),1)=0.d0
              dpm(ln(i,j),2)=dp(ln(i,j),2)+(dp(ln(i,j),1)-dpm(ln(i,j),1))*norm(i,1)/norm(i,2)
          end do
      end do
      
!---calculate convective velocity and dpdt in ale frame
      do i=1,ns(1)
          cv(i,1)=dp(i,1)-dpm(i,1)
          cv(i,2)=dp(i,2)-dpm(i,2)
          dpdtm(i)=dpdt(i)-cv(i,1)*dp(i,1)-cv(i,2)*dp(i,2)
      end do
      
      return
    end
!**********************************************************************
subroutine ale_acc(delttime,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,clen,npl,nnode,nelem,ns,nelm,ln,node,norm,jcb,&
                  &phi,pphi,phit,pphit,dp,dpds,dpdss,dpdt,cv,d2pm,d2pdtm)
!**********************************************************************
      implicit none
      integer i,j,npl,nnode,nelm,nelem(npl),ns(npl),ln(nelm,2)
      real*8 delttime,grav,mu1,mu2,porous_x1,porous_x2,porous_x3,clen,dptds,dpdnds,temp
      real*8 norm(nelm,2),jcb(nelm)
      real*8 node(nnode,2),phi(nnode),pphi(nnode),phit(nnode),pphit(nnode)
	  real*8 dp(nnode,2),dpds(nnode),dpdt(nnode),d2p(nnode,2),d2pdt(nnode)
      real*8 cv(ns(1),2),d2pm(ns(1),2),d2pdtm(ns(1)),dpdss(ns(1))
      real*8 uxx(ns(1)),uxy(ns(1)),uyy(ns(1)),dpdtdx(ns(1),2)

!---calculate the second order term of tse in el frame (dpdss=dpdnn=0 for linear element)
	do i=1,nelem(1)
	  do j=1,2
		if(ln(i,j) .eq. 1) then
            d2p(ln(i,j),1)=-pphit(nnode)
            dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
            dptds=(d2p(ln(i,j),1)-pphi(1)*dpdnds*norm(i,2)-(dpds(1)*dpdnds+pphit(1))*norm(1,1))/norm(1,2)
            d2p(ln(i,j),2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
        else if (ln(i,j) .eq. nelem(1)+1) then
            d2p(ln(i,j),1)=pphit(nelem(1)+2)
            dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
            dptds=(d2p(ln(i,j),1)-pphi(1)*dpdnds*norm(i,2)-(dpds(1)*dpdnds+pphit(1))*norm(1,1))/norm(1,2)
            d2p(ln(i,j),2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
		else
			dptds=(-0.5d0*phit(ln(i,1))+0.5d0*phit(ln(i,2)))/jcb(i)
			dpdnds=(-0.5d0*pphi(ln(i,1))+0.5d0*pphi(ln(i,2)))/jcb(i)
			d2p(ln(i,j),1)=(dptds+pphi(ln(i,j))*dpdnds)*norm(i,2)-(-dpds(ln(i,j))*dpdnds-pphit(ln(i,j)))*norm(i,1)
			d2p(ln(i,j),2)=(pphit(ln(i,j))+dpds(ln(i,j))*dpdnds)*norm(i,2)-(dptds+pphi(ln(i,j))*dpdnds)*norm(i,1)
        end if
        uxx(ln(i,j))=2.d0*dpdnds*norm(i,1)*norm(i,2)+dpdss(ln(i,j))*(norm(i,2)**2-norm(i,1)**2) !2.d0*dpdnds*norm(i,1)*norm(i,2)
        uxy(ln(i,j))=dpdnds*(norm(i,2)**2-norm(i,1)*2)-2.d0*dpdss(ln(i,j))*norm(i,2)*norm(i,1) !dpdnds*(norm(i,2)*norm(i,2)-norm(i,1)*norm(i,1))
        uyy(ln(i,j))=-2.d0*dpdnds*norm(i,1)*norm(i,2)-dpdss(ln(i,j))*(norm(i,2)**2-norm(i,1)**2) !-2.d0*dpdnds*norm(i,1)*norm(i,2)
		d2pdt(ln(i,j))=dp(ln(i,j),1)*d2p(ln(i,j),1)+dp(ln(i,j),2)*d2p(ln(i,j),2)-grav*dp(ln(i,j),2)
        dpdtdx(ln(i,j),1)=dp(ln(i,j),1)*uxx(ln(i,j))+dp(ln(i,j),2)*uxy(ln(i,j))+grav*norm(i,1)/norm(i,2)
        dpdtdx(ln(i,j),2)=dp(ln(i,j),1)*uxy(ln(i,j))+dp(ln(i,j),2)*uyy(ln(i,j))
	  end do
    end do

    do i=1,ns(1)
        if (node(i,1)>=porous_x1 .and. node(i,1)<=porous_x2)then
            d2pdt(i)=d2pdt(i)-mu1*dpdt(i)
            dpdtdx(i,1)=dpdtdx(i,1)-mu1*dp(i,1) !-mu2*dp(i,1) not yet corrected
            dpdtdx(i,2)=dpdtdx(i,2)-mu1*dp(i,2)
        else if (node(i,1)>=porous_x3)then
            d2pdt(i)=d2pdt(i)-mu1*dpdt(i)
            dpdtdx(i,1)=dpdtdx(i,1)-mu1*dp(i,1)
            dpdtdx(i,2)=dpdtdx(i,2)-mu1*dp(i,2)
        end if
    end do

!---calculate the mesh acceleration
      do i=1,nelem(1)+1
          d2pm(i,1)=d2p(i,1)-cv(i,1)*uxx(i)-cv(i,2)*uxy(i)
          d2pm(i,2)=d2p(i,2)-cv(i,1)*uxy(i)-cv(i,2)*uyy(i)
          d2pdtm(i)=d2pdt(i)-cv(i,1)*dpdtdx(i,1)-cv(i,2)*dpdtdx(i,2)
      end do
    
      return
    end
!**********************************************************************
subroutine ale_tse_sum(delttime,npl,nnode,ns,node,phi,dpm,dpdtm,d2pm,d2pdtm)
!**********************************************************************
      implicit none 
      integer i,npl,nnode,ns(npl)
      real*8 delttime
      real*8 node(nnode,2),phi(nnode),new_node(ns(1),2),new_phi(ns(1))
	  real*8 dpm(ns(1),2),dpdtm(ns(1)),d2pm(ns(1),2),d2pdtm(ns(1))

	do i=1,ns(1)
		new_phi(i)=phi(i)+delttime*dpdtm(i)+0.5d0*delttime**2*d2pdtm(i)
		new_node(i,1)=node(i,1)+delttime*dpm(i,1)+0.5d0*delttime**2*d2pm(i,1)
		new_node(i,2)=node(i,2)+delttime*dpm(i,2)+0.5d0*delttime**2*d2pm(i,2)
	end do

	do i=1,ns(1)
	    node(i,:)=new_node(i,:)
	    phi(i)=new_phi(i)
	end do

      return
    end
!**********************************************************************
      subroutine fs_smooth(nnrl,mdeg,n,y)
!**********************************************************************
      implicit none
      integer n,nnrl,mdeg,flag
      real*8 y(n)
      
      call savgol_filter(nnrl,nnrl,0,mdeg,n,y,flag)
      
      return
    end
!********************************************************************
subroutine pressure(icon,time,tho,grav,dep,npl,nnode,ns,node,phit,dp,pr,p_atm)
!********************************************************************
      implicit none
      integer  i,icon,npl,nnode,ns(npl)
      real*8 time,dep,tho,grav,p1,p2,p3,p_atm
      real*8 node(nnode,2),phit(nnode),dp(nnode,2),pr(nnode)
	  real*8 cp1(ns(1)),cp2(ns(1)),cp3(ns(1)),cp(ns(1))

!----atom pressure (bernoulli constant) on the free surface
	do i=icon,icon ! ns(1) !
	cp1(i)=tho*phit(i)
	cp2(i)=tho*0.5d0*(dp(i,1)**2+dp(i,2)**2)
	cp3(i)=tho*grav*(node(i,2)-dep)
	cp(i)=cp1(i)+cp2(i)+cp3(i)
	enddo
	p_atm=cp(icon)

!----pressure on zone 1 boundary
	do i=ns(1)+2,nnode-1
	p1=tho*phit(i)
	p2=tho*0.5d0*(dp(i,1)**2+dp(i,2)**2)
	p3=tho*grav*(node(i,2)-dep)
	pr(i)=p_atm-(p1+p2+p3)
	end do

      return
      end
!********************************************************************
subroutine domain(npl,nga,nfield,nnode,nelm,nelem,ns,ln,node,norm,jcb,phi,pphi,phit,pphit,&
				 &sha1,sha2,sh,wt,tho,grav,dep,p_atm,dp,pr)
!********************************************************************
      implicit none
      integer  i,j,k,l,m,il,ir,npl,nga,nfield,nnode,nelm,nelem(npl),ns(npl),ln(nelm,2)
	  real*8 tho,grav,dep,p_atm
	  real*8 hb,dx,dy,pi2,temp,p1,p2,p3
	  real*8 node(nnode,2),norm(nelm,2),jcb(nelm)
	  real*8 phi(nnode),pphi(nnode),phit(nnode),pphit(nnode),dp(nnode,2),pr(nnode)
	  real*8 dnode(nfield,2),dvx(nfield),dvy(nfield),dphit(nfield),dpr(nfield)
	  real*8 ker1(nfield,nnode),ker2(nfield,nnode)
      real*8 h(2),g(2),xfunc(10),yfunc(10),pxi1(2)
	  real*8 wt(nga),sha1(nga),sha2(nga),sh(2,nga)
	
	pi2=2.d0*dacos(-1.d0)

!----create domain point
	l=1
	do i=2,ns(1)-1
	  call bwloc(node(i,1),ns(npl-1)-ns(2),node(ns(2)+1:ns(npl-1),1),ns(2),il,ir)
	  hb=node(il,2)+(node(i,1)-node(il,1))/(node(ir,1)-node(il,1))*(node(ir,2)-node(il,2))+0.01d0 ! keep it a little far away from the boundary
	  dy=-node(i,2)/nelem(npl)
		do j=2,nelem(npl)
		  temp=node(i,2)+dy*(j-1)
			if(temp>hb)then
			dnode(l,1)=node(i,1)
			dnode(l,2)=temp !node(i,2)+dy*(j-1)
			l=l+1
			end if
		end do
	end do

!--- set a dummy node to useless dnode
	do i=l,nfield
	dnode(i,:)=dnode(1,:)
	end do

	ker1=0.d0
	ker2=0.d0
	dvx=0.d0
!----calculate x velocity by bie
	do i = 1,nfield
       do j=1,nelm
		 do m=1,nga
			xfunc(m)=sha1(m)*node(ln(j,1),1)+sha2(m)*node(ln(j,2),1)
			yfunc(m)=sha1(m)*node(ln(j,1),2)+sha2(m)*node(ln(j,2),2)
		 end do
		do k=1,2
         g(k)=0.d0
         h(k)=0.d0
          do m=1,nga
			temp=jcb(j)*sh(k,m)*wt(m)
            h(k)=h(k)+(((yfunc(m)-dnode(i,2))**2-(xfunc(m)-dnode(i,1))**2)*norm(j,1)-&
					&2.d0*(yfunc(m)-dnode(i,2))*(xfunc(m)-dnode(i,1))*norm(j,2))/&
					&((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)**2*temp
            g(k)=g(k)+(xfunc(m)-dnode(i,1))/((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)*temp
         end do
	     ker1(i,ln(j,k))=ker1(i,ln(j,k))+h(k)
         ker2(i,ln(j,k))=ker2(i,ln(j,k))+g(k)
		end do
       end do
      end do
	dvx=(matmul(ker2,pphi)-matmul(ker1,phi))/pi2

	ker1=0.d0
	ker2=0.d0
	dvy=0.d0
!----calculate y velocity by bie
	do i = 1,nfield
       do j=1,nelm
		 do m=1,nga
			xfunc(m)=sha1(m)*node(ln(j,1),1)+sha2(m)*node(ln(j,2),1)
			yfunc(m)=sha1(m)*node(ln(j,1),2)+sha2(m)*node(ln(j,2),2)
		 end do
		do k=1,2
         g(k)=0.d0
         h(k)=0.d0
          do m=1,nga
			temp=jcb(j)*sh(k,m)*wt(m)
            h(k)=h(k)+(((xfunc(m)-dnode(i,1))**2-(yfunc(m)-dnode(i,2))**2)*norm(j,2)-&
					&2.d0*(yfunc(m)-dnode(i,2))*(xfunc(m)-dnode(i,1))*norm(j,1))/&
					&((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)**2*temp
            g(k)=g(k)+(yfunc(m)-dnode(i,2))/((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)*temp
         end do
	     ker1(i,ln(j,k))=ker1(i,ln(j,k))+h(k)
         ker2(i,ln(j,k))=ker2(i,ln(j,k))+g(k)
		end do
       end do
      end do
	dvy=(matmul(ker2,pphi)-matmul(ker1,phi))/pi2
	
	ker1=0.d0
	ker2=0.d0
	dphit=0.d0
!----calculate partial potential over time by bie
      do i = 1,nfield
       do j=1,nelm
		 do m=1,nga
			xfunc(m)=sha1(m)*node(ln(j,1),1)+sha2(m)*node(ln(j,2),1)
			yfunc(m)=sha1(m)*node(ln(j,1),2)+sha2(m)*node(ln(j,2),2)
		 end do
		do k=1,2
         g(k)=0.d0
         h(k)=0.d0
          do m=1,nga
			temp=jcb(j)*sh(k,m)*wt(m)
            h(k)=h(k)+(-1.d0)/((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)*&
					&((xfunc(m)-dnode(i,1))*norm(j,1)+(yfunc(m)-dnode(i,2))*norm(j,2))*temp
            g(k)=g(k)+dlog(1.d0/((xfunc(m)-dnode(i,1))**2+(yfunc(m)-dnode(i,2))**2)**0.5d0)*temp	 
		 end do
	     ker1(i,ln(j,k))=ker1(i,ln(j,k))+h(k)
         ker2(i,ln(j,k))=ker2(i,ln(j,k))+g(k)
		end do
       end do
      end do
	dphit=(matmul(ker2,pphit)-matmul(ker1,phit))/pi2

!----calculate pressure distribution in domain
	do i=1,nfield
	p1=tho*dphit(i)
	p2=tho*0.5d0*(dvx(i)**2+dvy(i)**2)
	p3=tho*grav*(dnode(i,2)-dep)
	dpr(i)=p_atm-(p1+p2+p3)
	end do

	write(11,'(3000(1x,f15.7))') node(:,1),dnode(:,1)
	write(11,'(3000(1x,f15.7))') node(:,2),dnode(:,2)
	write(11,'(3000(1x,f15.7))') dp(:,1),dvx
	write(11,'(3000(1x,f15.7))') dp(:,2),dvy
	write(11,'(3000(1x,f15.7))') pr,dpr

      return
      end
!********************************************************************
      subroutine gauss(wt,rt,nga)
!********************************************************************
      integer nga
      real*8   wt(nga),rt(nga)

      select case(nga)
       case(3)
        wt(1)=0.55555555
        wt(2)=0.88888889
        wt(3)=0.55555555
        rt(1)=0.77459667
        rt(2)=0.d0
        rt(3)=-0.77459667
       case(4)
        wt(1)=0.65214515
        wt(2)=0.34785484
        wt(3)=0.34785484
        wt(4)=0.65214515
        rt(1)=0.33998104
        rt(2)=0.86113631
        rt(3)=-0.86113631
        rt(4)=-0.33998104
       case(5)
        wt(1)=0.23692689
        wt(2)=0.47862867
        wt(3)=0.56888889
        wt(4)=0.47862867
        wt(5)=0.23692689
        rt(1)=0.90617985
        rt(2)=0.53846931
        rt(3)=0.d0
        rt(4)=-0.53846931
        rt(5)=-0.90617985
	 case(6)
	  wt(1)=0.17132449
	  wt(2)=0.36076157
	  wt(3)=0.46791393
	  wt(4)=0.46791393
	  wt(5)=0.36076157
	  wt(6)=0.17132449
	  rt(1)=0.93246951
	  rt(2)=0.66120938
	  rt(3)=0.23861918
	  rt(4)=-0.23861918
	  rt(5)=-0.66120938
	  rt(6)=-0.9346951
       case(8)
        wt(1)=0.1012285362903763d0
        wt(2)=0.2223810344533745d0
        wt(3)=0.3137066458778873d0
        wt(4)=0.3626837833783620d0
        wt(8)=0.1012285362903763d0
        wt(7)=0.2223810344533745d0
        wt(6)=0.3137066458778873d0
        wt(5)=0.3626837833783620d0
        rt(1)=0.9602898564975363d0
        rt(2)=0.7966664774136267d0
        rt(3)=0.5255324099163290d0
        rt(4)=0.1834346424956498d0
        rt(8)=-0.9602898564975363d0
        rt(7)=-0.7966664774136267d0
        rt(6)=-0.5255324099163290d0
        rt(5)=-0.1834346424956498d0
       case(10)
        wt(1)=0.d06667134
        wt(2)=0.14945134
        wt(3)=0.21908636
        wt(4)=0.26926671
        wt(5)=0.29552422
        wt(10)=0.d06667134
        wt(9)=0.14945134
        wt(8)=0.21908636
        wt(7)=0.26926671
        wt(6)=0.29552422
        rt(1)=0.97390652
        rt(2)=0.86506336
        rt(3)=0.67940956
        rt(4)=0.43339539
        rt(5)=0.14887433
        rt(10)=-0.97390652
        rt(9)=-0.86506336
        rt(8)=-0.67940956
        rt(7)=-0.43339539
        rt(6)=-0.14887433
      end select

      return
      end
!********************************************************************
	subroutine converge(n,p1,p2,e1loc,e1,e2)
!********************************************************************
	implicit none
	integer i,n,k(1),e1loc
	real*8 p1(n),p2(n),e1,e2

	e1=maxval(dabs(p1-p2))
	e1loc=maxloc(dabs(p1-p2),1)

	e2=0.d0
	do i=1,n
	e2=e2+(p1(i)-p2(i))**2
	end do
	e2=dsqrt(e2/n)

	return
	end
!********************************************************************
	subroutine bwloc(px,n,x,ist,il,ir)
!********************************************************************
	implicit none
	integer i,n,ist,il,ir
	real*8 px,x(n)

	do i=1,n-1
		if(x(i)>px.and.x(i+1)<=px)then
		ir=ist+i
		il=ist+i+1
		goto 777
		end if
	end do
777 continue

	return
	end
!********************************************************************
	subroutine wave_spd(grav,omega,d,c)
!********************************************************************
	implicit none
	integer i
	real*8 k,k2,grav,omega,d,c,pi,f0,f1
	pi=dacos(-1.d0)

	k=1.d0
	do i=1,100
	f0=k*dtanh(k*d)-omega**2/grav
	f1=dtanh(k*d)+k-k*(dtanh(k*d)**2)
	k2=k-(f0)/(f1)
		if((k2-k)/k<=0.000001d0) then
		goto 717
		end if
	k=k2
	end do
	717 continue

	c=dsqrt(grav*dtanh(k*d)/k)

	return
	end
!********************************************************************
subroutine courant(time,delttime,nnode,nelm,ln,node,dp,jcb)
!********************************************************************
    implicit none
    integer i,j,cfloc,nnode,nelm,ln(nelm,2)
    real*8 time,delttime,u,v,ve
	real*8 node(nnode,2),dp(nnode,2),jcb(nelm),cn(nelm),cfl

  do i=1,nelm
    u=dsqrt(dp(ln(i,1),1)**2+dp(ln(i,1),2)**2)
    v=dsqrt(dp(ln(i,2),1)**2+dp(ln(i,2),2)**2)
    ve=max(u,v)
    cn(i)=0.5d0*ve*delttime/jcb(i)
  end do
  cfl=maxval(cn)
  cfloc=maxloc(cn,1)
  
	write(23,*) time,cfl

    if (cfl>=250.d0)then
    write(22,*) time,"cfl=",cfl,"@ element",cfloc
    stop
    end if    

	return
	end
!**********************************************************************
subroutine lobatto(n,x,w,a,b)
!**********************************************************************
implicit none
integer i,iter,n
real a,b,c,d,x(n),w(n)
real x0,x1,e,p1,p2,p3
real,parameter::pi=3.1415926
real,parameter::emax=1.0e-5
real,external::legendre

x(1)=-1.0
w(1)=2.0/n/(n-1)
x(n)=1.0
w(n)=w(1)

do i=2,n-1
iter=1
e=100.0
x0=(1.0-3.0*(n-2)/8.0/(n-1)**3)*cos((4.0*i-3)*pi/(4.0*n-3))
    do while (e>=emax.and.iter<=1000)
    p1=(legendre(n-2,x0)-x0*legendre(n-1,x0))*(n-1)/(1.0-x0**2)
    p2=(2.0*x0*p1-n*(n-1)*legendre(n-1,x0))/(1.0-x0**2)
    p3=(2.0*x0*p2-(n*(n-1)-2)*p1)/(1.0-x0**2)
    x1=x0-2.0*p1*p2/(2.0*p2**2-p1*p3)
    e=abs(x1-x0)
    iter=iter+1
    x0=x1
    end do
x(n-i+1)=x1
end do

do i=2,n-1
w(i)=2.0/n/(n-1)/legendre(n-1,x(i))**2
end do

c=(b-a)/2.0
d=(b+a)/2.0
do i=1,n
x(i)=c*x(i)+d
w(i)=c*w(i)
end do

return
end
!**********************************************************************
function legendre(n,x)
!**********************************************************************
real*8 fi,pi,pim1,pim2
real legendre
integer i
if (n.eq.0) then
    legendre=1
    elseif (n.eq.1) then
    legendre=x
    else
    pim1=1
    pi=x
        do i=2,n
        fi=i
        pim2=pim1
        pim1=pi
        pi=((i+i-1)*x*pim1-(i-1)*pim2)/fi
        end do
    legendre=pi
endif
end
      
subroutine savgol_filter(nl,nr,ld,m,n1,y,flag)
!-----------------------------------------------------------------------------------
! this routine is used to perform the savitzky-golay algorithm.
!-----------------------------------------------------------------------------------
!    nl:: input, integer, the number of leftward data points used.
!    nr:: input, integer, the number of rightward data points used.
!    ld:: input, integer, the order of the derivative desired.
!     m:: input, integer, the order of the smoothing polynomial.
!    n1:: input, integer, the number of data points.
! y(n1):: input/output, real values, the data to be smoothed.
!  flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! author: peng jun, 2019.03.26.
!-----------------------------------------------------------------------------------
! dependence:: subroutine savgol.
! -----------------------------------------------------------------------------------

    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m, n1
    real   (kind=8), intent(inout):: y(n1)
    integer(kind=4), intent(out):: flag
    ! local variables.
    integer(kind=4):: i, j, xl(nl+nr+1)
    real   (kind=8):: y0(n1), coef(nl+nr+1)

    xl(1) = 0

    y0 = y

    do i=1, nl

        xl(i+1) = -i

    end do

    do i=1, nr

        xl(1+nl+i) = nr-i+1

    end do

    call savgol(nl,nr,ld,m,coef,flag)

    if (flag/=0) return

    do i=1, n1-nr

        y(i) = 0.0

        do j=1, nl+nr+1

            if (i+xl(j) .gt. 0) then

                y(i) = y(i) + coef(j)*y0(i+xl(j))

            end if

        end do

    end do

    if (ld==0) then

        y(1:nl) = y0(1:nl)

        y(n1-nr+1:n1) = y0(n1-nr+1:n1)

    else 

        y(1:nl) = y(nl+1)

        y(n1-nr+1:n1) = y(n1-nr)
 
    end if

    return

end subroutine savgol_filter

subroutine savgol(nl,nr,ld,m,coef,flag)
!-----------------------------------------------------------------------------------
! this routine is used to calculate a set of savitzky-golay filter coefficients.
!-----------------------------------------------------------------------------------
!            nl:: input, integer, the number of leftward data points used.
!            nr:: input, integer, the number of rightward data points used.
!            ld:: input, integer, the order of the derivative desired.
!             m:: input, integer, the order of the smoothing polynomial.
! coef(nl+nr+1):: output, real values, calculated coefficents in wrap-around order.
!          flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! author: peng jun, 2019.03.20.
!-----------------------------------------------------------------------------------
! dependence:: subroutine ludcmp;
!              subroutine lubksb.
!-----------------------------------------------------------------------------------
! reference: press et al, 1986. numberic recipes in fortran 77, 
!            the art of scientific computing, second edition. 
! note: this subroutine is remodified from page.646 in press et al.
! -----------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m
    real   (kind=8), intent(inout):: coef(nl+nr+1)
    integer(kind=4), intent(out):: flag
    ! local variables.
    integer(kind=4):: imj, ipj, k, kk, mm, indx(m+1)

    real   (kind=8):: d, fac, summ, a(m+1,m+1), b(m+1)

    flag = 0

    if (nl < 0 .or. nr < 0 .or. ld > m .or. nl+nr < m) then

        flag = 1

        return

    end if

    do ipj=0, 2*m

        summ = 0.0

        if (ipj .eq. 0) summ = 1.0

        do k=1, nr

            summ = summ + (float(k))**ipj

        end do

        do k=1, nl

            summ = summ + (float(-k))**ipj

        end do

        mm = min(ipj, 2*m-ipj)

        do imj=-mm, mm, 2

            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = summ

        end do

    end do

    call ludcmp(a,m+1,indx,d,flag)

    if (flag .ne. 0) return

    b = 0.0

    b(ld+1) = 1.0

    call lubksb(a,m+1,indx,b)

    coef = 0.0

    do k=-nl, nr

        summ = b(1)

        fac = 1.0

        do mm=1, m

            fac = fac * k

            summ = summ + b(mm+1) * fac

        end do

        kk = mod(nl+nr+1-k, nl+nr+1) + 1

        coef(kk) = summ

    end do

    return

end subroutine savgol

subroutine lubksb(a,n,indx,b)
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, the lu decomposition of a matrix.
!       n:: input, integer, the dimenstion of the matrix.
! indx(n):: input,  integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!    b(n):: output, real values, the solution vector x for 
!                   linear equations a*x=b.
!-------------------------------------------------------------------------
! author: peng jun, 2019.03.18.
!-------------------------------------------------------------------------
! dependence:: no.--------------------------------------------------------
!-------------------------------------------------------------------------
! reference: press et al, 1986. numberic recipes in fortran 77, 
!            the art of scientific computing, second edition. 
! note: this subroutine is remodified from page.39 in press et al.
! -------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n, indx(n)
    real   (kind=8), intent(in):: a(n,n)
    real   (kind=8), intent(inout):: b(n)
   ! local variables.
    integer(kind=4):: i, ii, j, ll
    real   (kind=8):: summ

    ii = 0

    do i=1, n

        ll = indx(i)

        summ = b(ll)

        b(ll) = b(i)

        if (ii .ne. 0) then

            do j=ii, i-1

                summ = summ - a(i,j) * b(j)

            end do

        else if (summ .ne. 0.0) then
            
            ii = i

        end if

        b(i) = summ

    end do

    do i=n, 1, -1

        summ = b(i)

        do j=i+1, n

            summ = summ - a(i,j) * b(j)

        end do

        b(i) = summ / a(i,i)

    end do

    return

end subroutine lubksb

subroutine ludcmp(a,n,indx,d,flag)
!-------------------------------------------------------------------------
!this routine is used in combination with lubksb to solve 
!linear equations or invert a matrix.
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, a matrix to be decomposed.
!       n:: input, integer, the dimension of the matrix.
! indx(n):: output, integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!       d:: output, integer, output as 1 or -1 depending on whether 
!           the number of row interchanges was even or odd.
!    flag:: output, integer, error message, 0=success, 1=singular matrix.
!-------------------------------------------------------------------------
! author: peng jun, 2019.03.20.
!-------------------------------------------------------------------------
! dependence:: no.--------------------------------------------------------
!-------------------------------------------------------------------------
! reference: press et al, 1986. numberic recipes in fortran 77, 
!            the art of scientific computing, second edition. 
! note: this subroutine is remodified from page.38 in press et al.
! ------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n
    integer(kind=4), intent(out):: indx(n), flag
    real   (kind=8), intent(inout):: a(n,n)
    real   (kind=8), intent(out):: d
    ! local variables.
    integer(kind=4):: i, j, k, imax
    real   (kind=8):: aamax, dum, summ, vv(n)
   
    indx = 0

    flag = 0

    d = 1.0

    do i=1, n

        aamax = 0.0

        do j=1, n

            if (abs(a(i,j)) .gt. aamax)  aamax = abs(a(i,j))

        end do

        if (aamax .eq. 0.0) then 
 
            flag = 1

            return

        end if

        vv(i) = 1.0/aamax

    end do

    do j=1, n

        do i=1, j-1

            summ = a(i,j)

            do k=1, i-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

        end do

        aamax = 0.0

        do i=j, n

            summ = a(i,j)

            do k=1, j-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

            dum = vv(i) * abs(summ)

            if (dum .ge. aamax) then
       
                imax = i
          
                aamax = dum

            end if

        end do

        if (j .ne. imax) then

            do k=1, n

                dum = a(imax,k)

                a(imax,k) = a(j,k)

                a(j,k) = dum
  
            end do

            d = -d

            vv(imax) = vv(j)

        end if

        indx(j) = imax

        if (a(j,j) .eq. 0.0) a(j,j) = tiny(0.0d+00)

        if (j .ne. n) then

            dum = 1.0 / a(j,j)

            do i=j+1, n

                a(i,j) = a(i,j) * dum

            end do

        end if

    end do

    return

    end subroutine ludcmp    
      subroutine seval (n,u,x,y,b,c,d,v_0,v_1,v_2)
!------------------------------------------------------------------------
!     evaluate a cubic spline interpolation of a discrete function f(x),
!     given in n points x(i), y(i). the b, c and d coefficients defining
!     the best cubic spline for the given points, are calculated before
!     by the spline subroutine.
!
!     inputs:
!     n       number of points of curve y = f(x)
!     u       abscissa of point to be interpolated
!     x,y     tables of dimension n, storing the coordinates
!             of curve f(x)
!     b,c,d   tables storing the coefficients defining the
!             cubic spline
!
!     outputs:
!     seval   interpolated value
!             = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
!             with dx = u-x(i), u between x(i) and x(i+1)
!
!     reference :
!     forsythe,g.e. (1977) computer methods for mathematical
!     computations. prentice-hall,inc.
!------------------------------------------------------------------------
      real *8 b(n),c(n),d(n),x(n),y(n),u,dx
      real *8 v_0,v_1,v_2
      data i/1/

!     binary search

      if (i.ge.n) i = 1
      if (u.lt.x(i)) go to 101
      if (u.le.x(i+1)) go to 301
101 i = 1
      j = n+1
201 k = (i+j)/2
      if (u.lt.x(k)) j = k
      if (u.ge.x(k)) i = k
      if (j.gt.i+1) go to 201

!     spline evaluation

   301 dx = u-x(i)
!      seval = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      v_0 = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      v_1 = b(i)+2.d0*c(i)*dx+3.d0*d(i)*dx**2
      v_2 = 2.d0*c(i)+6.d0*d(i)*dx
      return
      end

      subroutine spline (n,x,y,b,c,d)
!---------------------------------------------------------------------
!     this subroutine calculates the coefficients b,c,d of a cubic
!     spline to best approximate a discrete function given by n points
!
!     inputs:
!     n       number of given points
!     x,y     vectors of dimension n, storing the coordinates
!             of function f(x)
!
!     outputs:
!     a,b,c   vectors of dimension n, storing the coefficients
!             of the cubic spline
!
!     reference:
!     forsythe,g.e. (1977) computer methods for mathematical
!     computations. prentice-hall,inc.
!---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      dimension b(n),c(n),d(n),x(n),y(n)
      nm1 = n-1
      if (n.lt.2) return
      if (n.lt.3) go to 501

!     build the tridiagonal system
!     b (diagonal), d (upperdiagonal) , c (second member)

      d(1) = x(2)-x(1)
      c(2) = (y(2)-y(1))/d(1)
      do i = 2,nm1
      d(i) = x(i+1)-x(i)
      b(i) = 2.d0*(d(i-1)+d(i))
      c(i+1) = (y(i+1)-y(i))/d(i)
      c(i) = c(i+1)-c(i)
      end do

!     conditions at limits
!     third derivatives obtained by divided differences

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.d0
      c(n) = 0.d0
      if (n.eq.3) go to 151
      c(1) = c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)*d(1)/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))

!     forward elimination

151 do i = 2,n
      t = d(i-1)/b(i-1)
      b(i) = b(i)-t*d(i-1)
      c(i) = c(i)-t*c(i-1)
      end do

!     back substitution

      c(n) = c(n)/b(n)
      do l = 1,nm1
      i = n-l
      c(i) = (c(i)-d(i)*c(i+1))/b(i)
      end do

!     coefficients of 3rd degree polynomial

      b(n) = (y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.d0*c(n))
      do i = 1,nm1
      b(i) = (y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.d0*c(i))
      d(i) = (c(i+1)-c(i))/d(i)
      c(i) = 3.d0*c(i)
      end do
      c(n) = 3.d0*c(n)
      d(n) = d(nm1)
      return

!     cas n = 2

501 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.d0
      d(1) = 0.d0
      b(2) = b(1)
      c(2) = 0.d0
      d(2) = 0.d0
      return
    end
!********************************************************************
 subroutine periodic(h,wave_period,wave_height,wave_amplitude)
!********************************************************************
! the subroutine is used to give the stroke of the wavemaker
 implicit none
 integer i
 real*8 pi,k0,k1,hs,s,e
 real*8 h,wave_period,wave_height,wave_amplitude
 pi=dacos(-1.d0)
    k0 = 1.d0
    e = 1.d0
    do while (e > 1.e-7)
        k1 = k0-(k0*dtanh(k0*h)-(2.d0*pi/wave_period)**2/9.81)/(dtanh(k0*h)+k0-k0*(dtanh(k0*h))**2)
        e = dabs((k1-k0)/k0)
        hs = 2.d0*(dcosh(2.d0*k1*h)-1.d0)/(dsinh(2.d0*k1*h)+2.d0*k1*h)
        s = wave_height/hs
        wave_amplitude = 0.5d0*s
        k0 = k1
    end do
    
    return
    end
!********************************************************************
    subroutine solitary(h0,wave_h,g,npp,total_nc,total_step,dt,dis,vb_x,acc)
!********************************************************************
    implicit none
    integer i,npp,step,total_nc,total_step,no(npp)
    real*8 h0,wave_h,wave_k,wave_c,g,dt,x0
    real*8 vb_x(total_nc,0:total_step-1),vb_z(total_nc,0:total_step-1),x(total_nc,0:total_step-1)
    real*8 dis(total_nc,0:total_step-1),acc(total_nc,0:total_step-1)
    
    call fenton_parameters(h0,wave_h,wave_k,wave_c,g)
    
    x0=-5.d0*wave_c
    x=0.d0
    do i=1,npp
        no(i)=i
    end do

    do step = 0,total_step-1
        call wmbc(total_nc,total_step,vb_x,vb_z,step,dt,wave_h,wave_k,wave_c,h0,x0,x,g,no,npp)
    end do
    
    call get_wmk_dis_acc(total_nc,total_step,vb_x,dt,dis,acc)

    return
    end
!********************************************************************
    subroutine get_wmk_dis_acc(total_nc,total_step,vb_x,dt,dis,acc)
!********************************************************************
    implicit none
    integer i,j,total_nc,total_step
    real*8 vb_x(total_nc,0:total_step-1),dt,dis(total_nc,0:total_step-1),acc(total_nc,0:total_step-1)
    
    dis=0.d0
    acc=0.d0
    do i=1,total_nc
        do j=1,total_step-1
            dis(i,j)=dis(i,j-1)+0.5d0*dt*(vb_x(i,j-1)+vb_x(i,j))
            if (j==total_step-1)then
                acc(i,j)=acc(i,j-1)
            else if (j==0)then
                acc(i,j)=(vb_x(i,j+1)-vb_x(i,j))/dt
            else
                acc(i,j)=0.5d0*(vb_x(i,j+1)-vb_x(i,j-1))/dt
            end if
        end do
    end do

    return
    end
!********************************************************************
    subroutine wmbc(total_nc,total_step,vb_x,vb_z,step,dt,wave_h,wave_k,wave_c,h0,x0,x,g,no,npp)
!********************************************************************
    ! the subroutine is used to give the velocity of the wavemaker
    ! x(total_nc,total_step) is the coordinate of the observation point at time t
    ! x0 = initial position of wave crest
    implicit none
   integer j,total_nc,total_step,step
   integer::zz,npp,no(npp)
   real*8 vb_x(total_nc,0:total_step-1),vb_z(total_nc,0:total_step-1),dt,wave_h,wave_k,wave_c,h0,x0
   real*8 x(total_nc,0:total_step-1),zeta,kx,t,alpha,kk,cc,g,xx,s,xxx
   
   t=step*dt
   alpha=wave_h/h0
   kk=wave_k*h0
   cc=wave_c/(g*h0)**0.5
   
   do j=1,total_nc
      do zz=1,npp
          
         if(j.eq.no(zz)) then
             xxx=x(j,step) ! xxx: xi = position of the wave paddle at time t (step)
             xx=(xxx-x0-wave_c*t)/h0 ! xx: x = xi - ct -x0
             s=1.d0/cosh(kk*xx) ! sech(k*x)
              zeta=s**2*alpha
              zeta=zeta+(-.75*s**2+.75*s**4)*alpha**2
              zeta=zeta+(.625*s**2-1.8875*s**4+1.2625*s**6)*alpha**3
           zeta=zeta+(-1.36817*s**2+3.88033*s**4-4.68304*s**6+2.17088*s**8)*alpha**4
        zeta=zeta+(1.86057*s**2-7.45136*s**4+12.7637*s**6-11.4199*s**8+4.24687*s**10)*alpha**5
    zeta=zeta+(-2.57413*s**2+13.2856*s**4-31.1191*s**6+40.1068*s**8-28.4272*s**10+8.728*s**12)*alpha**6
zeta=zeta+(3.4572*s**2-22.782*s**4+68.258*s**6-116.974*s**8+120.49*s**10-71.057*s**12+18.608*s**14)*alpha**7
zeta=zeta+(-4.6849*s**2+37.67*s**4-139.28*s**6+301.442*s**8-411.416*s**10+355.069*s**12-180.212*s**14+41.412*s**16)*alpha**8
zeta=zeta+(6.191*s**2-60.57*s**4+269.84*s**6-712.125*s**8+1217.98*s**10-1384.37*s**12+1023.07*s**14-450.29*s**16+90.279*s**18)*alpha**9
              zeta=zeta*h0
              vb_x(j,step)=wave_c*zeta/(h0+zeta) ! the velocity of the wavemaker
              goto 111
              
         else
            vb_x(j,step)=0.d0
         end if
         
      end do
     111  vb_z(j,step)=0.d0
   end do
   
    end
!********************************************************************
 subroutine fenton_parameters(h0,wave_h,wave_k,wave_c,g)
!********************************************************************
! the subroutine is used to give the parameters of the solitary wave by fenton theory
! wave_h = wave height
! h0 = still water depth
! g = gravity
      real*8 wave_h,wave_k,wave_c,h0,g,alpha
      alpha=wave_h/h0
      wave_k=sqrt(3./4.d0*alpha)
      wave_k=wave_k*(1-.625*alpha+.554688*alpha**2-.561535*alpha**3+.567095*alpha**4-.602969*alpha**5 &
                  +.624914*alpha**6-.670850*alpha**7+.700371*alpha**8)
      wave_c=1+alpha
      wave_c=wave_c-.05*alpha**2-.0428571*alpha**3-.0342857*alpha**4-.0315195*alpha**5-.0292784*alpha**6 &
          -.0268451*alpha**7-.0302634*alpha**8-.0219347*alpha**9
      wave_k=wave_k/h0
      wave_c=wave_c**0.5*(g*h0)**0.5
    end

