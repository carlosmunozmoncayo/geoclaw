! =====================================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================================
!
!solve Riemann problems for the 1D shallow water equations
!    with source-term resulting from variable topography b(x,t)
!    (h)_t + (u h)_x = 0
!    (uh)_t + (uuh + 0.5*gh**2)_x = -g*h*b_x
!
!
!
!     On input,
!     ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
!     On output, wave contains the fwaves/s,
!                s the speeds,
!                amdq the  left-going flux difference  A**- \Delta q
!               apdq the right-going flux difference  A**+ \Delta q
!
!     Note that the i'th Riemann problem has left state qr(i-1,:)
!     #                                    and right state ql(i,:)
!     From the basic clawpack routine step1, rp1 is called with ql=qr=q.
!
!      This is for use with the Riemann solver(s) used in GeoClaw
!        but for 1D problems. That is, it calls the same solver for left and
!        right states that is used for 2d problems.
!
!        This routine deals with dry-state problems over topography
!        by testing a wall boundary condition, like that done in GeoClaw 2d
!
!        to call other point-wise Riemann solvers alter the call on line 177

    use geoclaw_module, only: dry_tolerance, grav, rho

    !Parameters for Boussinesq equations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !For the moment, the only system that requires changes in the RS
    !is EDC HypRel (boussEquations = -1)
    use bouss_module, only: boussEquations, useBouss, boussMinDepth
    use bouss_module, only: riemann_normal_HypRel_EDC, riemann_normal_HypRel_BBBD
    use bouss_module, only: c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
    use bouss_module, only: gamma => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    ! Input arguments
    integer, intent(in) :: maxmx,meqn,mwaves,mbc,mx,maux

    double precision, intent(in), dimension(meqn, 1-mbc:maxmx+mbc) :: ql,qr
    double precision, intent(in), dimension(maux, 1-mbc:maxmx+mbc) :: auxl,auxr

    ! Output arguments
    double precision, intent(out) :: s(mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out) :: fwave(meqn, mwaves, 1-mbc:maxmx+mbc)
    double precision, intent(out), dimension(meqn, 1-mbc:maxmx+mbc) :: amdq,apdq

    !Local
    integer :: m,i,mw,maxiter
    double precision :: hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
    double precision :: bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    double precision :: hstartest,hstarHLL,sLtest,sRtest
    double precision :: wall(2), fw(3,3), sw(3)
    double precision :: g,drytol,pL,pR

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Local EDC HypRel Boussinesq equations
    double precision :: wL, wR
    double precision :: fwEDC(5,5), swEDC(5) !To use with 2D normal RS
    !In 1D, these should have shapes (4,4) and (4),
    !therefore the following 2D information should be discarded (should be 0 anyway):
    ! The 3rd component of all the f-waves (fwEDC(3,:)) 
    ! The 4th f-wave (fwEDC(:,4)) and the 4th wave speed (swEDC(4))
    integer :: skip_fourth_idx(6)
    integer :: skip_fifth_idx(6)
    !Local BBBD HypRel Boussinesq equations
    double precision :: sgl, sgr, pbl, pbr
    double precision :: fwBBBD(7,7), swBBBD(7) !To use with 2D normal RS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    skip_fourth_idx = (/ 1, 2, 3, 5, 6, 7 /)
    skip_fifth_idx = (/ 1, 2, 3, 4, 6, 7 /)
    g=grav
    drytol=dry_tolerance

    !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
              s(mw,i)=0.d0
              do m=1,meqn
                 fwave(m,mw,i)=0.d0
              enddo
         enddo


         !skip problem if in a completely dry area
         if (qr(1,i-1).le.drytol.and.ql(1,i).le.drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(1,i-1)
         hR = ql(1,i)
         huL = qr(2,i-1)
         huR = ql(2,i)
         bL = auxr(1,i-1)
         bR = auxl(1,i)

         !Solve Riemann problem for BTEqs or SWEs (compute fwaves and speeds)
         if (useBouss(i).and. useBouss(i-1) .and. min(hL,hR)>drytol &
            .and. boussEquations==-1) then
               !EDC HypRef Boussinesq-type equations
               !If we get here, the solver should return 4 f-waves
               !of 4 components each, and 4 wave speeds.
               !Also, q has 4 components: h,hu,hw,hp
               uL = huL/hL
               uR = huR/hR
               wL = qr(3,i-1)/hL
               wR = ql(3,i)/hR
               pL = qr(4,i-1)/hL
               pR = ql(4,i)/hR
               call riemann_normal_HypRel_EDC(2,3,g,gamma,c_sq, &
                                    hL, uL, 0.d0, wL, pL, bL, c_sq, &
                                    hR, uR, 0.d0, wR, pR, bR, c_sq, &
                                    fwEDC, swEDC)
               
               
               !swEDC(3) and fwEDC(:,3) have no use in 1D
               do mw=1,mwaves
                  s(mw,i) = swEDC(skip_fourth_idx(mw))
                  fwave(1:2,mw,i) = fwEDC(1:2,skip_fourth_idx(mw)) 
                  !fwEDC(3,:) has no use in 1D
                  fwave(3:4,mw,i) = fwEDC(4:5,skip_fourth_idx(mw))
               end do
         
         else if (useBouss(i).and. useBouss(i-1) .and. min(hL,hR)>drytol &
            .and. boussEquations==-2) then
               !BBBD HypRef Boussinesq-type equations
               !If we get here, the solver should return 6 f-waves
               !of 6 components each, and 6 wave speeds.
               !Also, q has 6 components: h,hu,hw,h\sigma,hp,hp_b
               uL = huL/hL
               uR = huR/hR
               wL = qr(3,i-1)/hL
               wR = ql(3,i)/hR
               sgL = qr(4,i-1)/hL
               sgR = ql(4,i)/hR
               pL = qr(5,i-1)/hL
               pR = ql(5,i)/hR
               pbl = qr(6,i-1)/hL
               pbr = ql(6,i)/hR
               call riemann_normal_HypRel_BBBD(mu=2, nv=3, grav=g, c_sq=c_sq, &
                  hl=hl, ul=ul, vl=0.d0, wl=wl, sgl=sgl, pl=pl, pbl=pbl, bl=bl, csq_l=c_sq,&
                  hr=hr, ur=ur, vr=0.d0, wr=wr, sgr=sgr, pr=pr, pbr=pbr, br=br, csq_r=c_sq, &
                  fw=fwBBBD, sw=swBBBD)
               
               !swBBBD(3) and fwBBBD(:,3) have no use in 1D
               do mw=1,mwaves
                  s(mw,i) = swBBBD(skip_fifth_idx(mw))
                  fwave(1:2,mw,i) = fwBBBD(1:2,skip_fifth_idx(mw)) 
                  !fwBBBD(3,:) has no use in 1D
                  fwave(3:6,mw,i) = fwBBBD(4:7,skip_fifth_idx(mw)) 
               end do
         else
            !Solve shallow water equations
            hvL=0.d0
            hvR=0.d0

            !check for wet/dry boundary
            sE1= 1.d99
            sE2=-1.d99
            wall(2) = 1.d0
            wall(1) = 1.d0
            if (hR.le.drytol) then
               sLtest=min(-sqrt(g*hL),huL/hL-sqrt(g*hL)) !what would be the Einfeldt speed of wall problem
               hstartest=hL-(huL/sLtest) !what would be middle state in approx Riemann solution
               if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
                  wall(2)=0.d0
                  hR=hL
                  huR=-huL
                  bR=bL
                  phiR=phiL
               endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
               sRtest=max(sqrt(g*hR),huR/hR+sqrt(g*hR)) !what would be the Einfeldt speed of wall
               hstartest= hR-(huR/sRtest) !what would be middle state in approx Rimeann solution
               if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
                  wall(1)=0.d0
                  hL=hR
                  huL=-huR
                  bL=bR
                  phiL=phiR
               endif
            endif

            if (hR.gt.drytol) then
               uR=huR/hR
               vR=hvR/hR
               phiR = 0.5d0*g*hR**2 + huR**2/hR
            else
               hR = 0.d0
               huR = 0.d0
               uR = 0.d0
               vR = 0.d0
               phiR = 0.d0
               sE2 = max(sE2,huL/hL+2.d0*sqrt(g*hL))
            endif

            if (hL.gt.drytol) then
               uL=huL/hL
               vL=hvL/hL
               phiL = 0.5d0*g*hL**2 + huL**2/hL
            else
               hL=0.d0
               huL=0.d0
               uL=0.d0
               vL=0.d0
               phiL = 0.d0
               sE1 = min(sE1,huR/hR - 2.d0*sqrt(g*hR))
            endif

            !determine wave speeds
            sL=uL-sqrt(g*hL) ! 1 wave speed of left state
            sR=uR+sqrt(g*hR) ! 2 wave speed of right state

            uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
            chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
            sRoe1=uhat-chat ! Roe wave speed 1 wave
            sRoe2=uhat+chat ! Roe wave speed 2 wave

            sE1 = min(sE1,min(sL,sRoe1)) ! Eindfeldt speed 1 wave
            sE2 = max(sE2,max(sR,sRoe2)) ! Eindfeldt speed 2 wave

            !--------------------end initializing...finally----------
            !solve Riemann problem.

            maxiter = 1

            ! In case there is no pressure forcing
            pL = 0.d0
            pR = 0.d0

            call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
               huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR, &
               sE1,sE2,drytol,g,rho,sw,fw)

            !         call riemann_ssqfwave(maxiter,meqn+1,mwaves+1,hL,hR,huL,huR, &
            !         &  hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,rho,sw,fw)

            !         call riemann_fwave(meqn+1,mwaves+1,hL,hR,huL,huR,hvL,hvR, &
            !           &   bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,rho,sw,fw)

            s(1,i) = sw(1)*wall(1)
            s(2,i) = sw(3)*wall(2)

            do m=1,2 !meqn This is hard-coded for the moment
               fwave(m,1,i)=fw(m,1)*wall(1)
               fwave(m,2,i)=fw(m,3)*wall(2)
               if (sw(2)>0.0) then
                  fwave(m,2,i) = fwave(m,2,i) + fw(m,2)*wall(2)
               else
                  fwave(m,1,i) = fwave(m,1,i) + fw(m,2)*wall(1)
               endif
            enddo
         end if

 30      continue
      enddo


      do i=2-mbc,mx+mbc
         do m=1,meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do  mw=1,mwaves
               if (s(mw,i).lt.0.d0) then
                  amdq(m,i) = amdq(m,i) + fwave(m,mw,i)
               elseif (s(mw,i).gt.0.d0) then
                  apdq(m,i) = apdq(m,i) + fwave(m,mw,i)
               else
                  amdq(m,i) = amdq(m,i) + .5d0*fwave(m,mw,i)
                  apdq(m,i) = apdq(m,i) + .5d0*fwave(m,mw,i)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine rp1


