!======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,&
                      ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!======================================================================
!
! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!     with topography:
!     #        h_t + (hu)_x + (hv)_y = 0                           #
!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.

!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
!  From the basic clawpack routines, this routine is called with
!     ql = qr
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updated to also solve a hyperbolic-dispersive system, presented by
! Escalante, Dumbser, and Castro (EDC)

      use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
      use geoclaw_module, only: earth_radius, deg2rad,sea_level

      use amr_module, only: mcapa

      use storm_module, only: pressure_forcing, pressure_index

      !Parameters for Boussinesq equations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !For the moment, the only system that requires changes in the RS
      !is EDC HypRel (boussEquations = -1)
      use bouss_module, only: boussEquations
      use bouss_module, only: riemann_normal_HypRel_EDC
      use bouss_module, only: riemann_normal_HypRel_BBBD
      use bouss_module, only: c_sq => boussEDCcsq 
      use bouss_module, only: csq_index => bouss_csq_index
      use bouss_module, only: dr_index => boussDecayRate_index
      !Reference hyperbolic relaxation parameter (Just used for EDC)
      !Approximate SGN (3/2) or Sainte-Marie equations (2)
      use bouss_module, only: gamma => boussEDCgamma 
      use bouss_module, only: small_tol_bouss, large_tol_bouss
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !input
      !mwaves = meqn = 5 (7) for EDC (BBBD) HypRel
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local only
      integer m,i,mw,maxiter,mu,nv
      double precision wall(3)
      double precision fw(3,3)
      double precision sw(3)

      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL,pL,pR
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision s1m,s2m
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest
      double precision tw,dxdc

      logical rare1,rare2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Local HypRel Boussinesq equations
      double precision :: wL, wR !For EDC and BBBD HypRel
      double precision :: sgL, sgR, pbL, pbR !For BBBD HypRel
      double precision :: fwHypRel(meqn,mwaves), swHypRel(mwaves)
      double precision :: csq_l, csq_r, drL, drR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! In case there is no pressure forcing
      pL = 0.d0 !Used also for EDC HypRel
      pR = 0.d0 !Used also for EDC HypRel

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

         !set normal direction
         if (ixy.eq.1) then
            mu=2
            nv=3
         else
            mu=3
            nv=2
         endif

         !zero (small) negative values if they exist
         if (qr(1,i-1).lt.0.d0) then
               qr(:,i-1)=0.d0
         endif

         if (ql(1,i).lt.0.d0) then
               ql(:,i)=0.d0
         endif

         !skip problem if in a completely dry area
         if (qr(1,i-1) <= drytol .and. ql(1,i) <= drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(1,i-1) 
         hR = ql(1,i) 
         huL = qr(mu,i-1) 
         huR = ql(mu,i) 
         bL = auxr(1,i-1)
         bR = auxl(1,i)
         csq_l = auxr(csq_index,i-1)
         csq_r = auxl(csq_index,i)
         drL = auxr(dr_index,i-1)
         drR = auxl(dr_index,i)
         !No pressure forcing, I'll recycle these variables for 
         !non-hydrostatic pressure 
         !if (pressure_forcing) then
         !    pL = auxr(pressure_index, i-1)
         !    pR = auxl(pressure_index, i)
         !end if

         hvL=qr(nv,i-1) 
         hvR=ql(nv,i)
        

         if (min(hl,hr)<drytol .or. min(drL,drR)<=small_tol_bouss &
            .or. boussEquations==0) then
         !Solve SWEs
            !check for wet/dry boundary
            if (hR.gt.drytol) then
               uR=huR/hR
               vR=hvR/hR
               phiR = 0.5d0*g*hR**2 + huR**2/hR
            else
               hR = 0.d0
               huR = 0.d0
               hvR = 0.d0
               uR = 0.d0
               vR = 0.d0
               phiR = 0.d0
            endif

            if (hL.gt.drytol) then
               uL=huL/hL
               vL=hvL/hL
               phiL = 0.5d0*g*hL**2 + huL**2/hL
            else
               hL=0.d0
               huL=0.d0
               hvL=0.d0
               uL=0.d0
               vL=0.d0
               phiL = 0.d0
            endif

            wall(1) = 1.d0
            wall(2) = 1.d0
            wall(3) = 1.d0
            if (hR.le.drytol) then
               call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,&
                           rare1,rare2,1,drytol,g)
               hstartest=max(hL,hstar)
               if (hstartest+bL.lt.bR) then 
               !right state should become ghost values that mirror left for wall problem
                  wall(2)=0.d0
                  wall(3)=0.d0
                  hR=hL
                  huR=-huL
                  bR=bL
                  phiR=phiL
                  uR=-uL
                  vR=vL
               elseif (hL+bL.lt.bR) then
                  bR=hL+bL
               endif
            elseif (hL.le.drytol) then ! right surface is lower than left topo
               call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,&
                           rare1,rare2,1,drytol,g)
               hstartest=max(hR,hstar)
               if (hstartest+bR.lt.bL) then  
               !left state should become ghost values that mirror right
                  wall(1)=0.d0
                  wall(2)=0.d0
                  hL=hR
                  huL=-huR
                  bL=bR
                  phiL=phiR
                  uL=-uR
                  vL=vR
               elseif (hR+bR.lt.bL) then
                  bL=hR+bR
               endif
            endif

            !determine wave speeds
            sL=uL-sqrt(g*hL) ! 1 wave speed of left state
            sR=uR+sqrt(g*hR) ! 2 wave speed of right state

            uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
            chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
            sRoe1=uhat-chat ! Roe wave speed 1 wave
            sRoe2=uhat+chat ! Roe wave speed 2 wave

            sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
            sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

            !--------------------end initializing...finally----------
            !solve Riemann problem.

            maxiter = 1

            call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
                  huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,&
                  drytol,g,rho,sw,fw)

            !eliminate ghost fluxes for wall
            do mw=1,3
               sw(mw)=sw(mw)*wall(mw)

                  fw(1,mw)=fw(1,mw)*wall(mw) 
                  fw(2,mw)=fw(2,mw)*wall(mw)
                  fw(3,mw)=fw(3,mw)*wall(mw)
            enddo

            do mw=1,3
               s(mw,i)=sw(mw)
               fwave(1,mw,i)=fw(1,mw)
               fwave(mu,mw,i)=fw(2,mw)
               fwave(nv,mw,i)=fw(3,mw)
               !            write(51,515) sw(mw),fw(1,mw),fw(2,mw),fw(3,mw)
               !515         format("++sw",4e25.16)
            enddo
         else if (boussEquations==-1) then
         !EDC HypRel Boussinesq-type equations
            uR=huR/hR
            vR=hvR/hR
            uL=huL/hL
            vL=hvL/hL
            wL = qr(4,i-1)/hL
            wR = ql(4,i)/hR
            pL = qr(5,i-1)/hL
            pR = ql(5,i)/hR
   
            call riemann_normal_HypRel_EDC(mu,nv,g,gamma,c_sq,&
                                 hL, uL, vL, wL, pL, bL, csq_l, &
                                 hR, uR, vR, wR, pR, bR, csq_r, &
                                 fwHypRel, swHypRel)

            ! Computing f-waves and speeds
            do mw = 1, mwaves
                fwave(:, mw, i) = fwHypRel(:,mw)
                s(mw,i) = swHypRel(mw)
            END do

         else if (boussEquations==-2) then
         !BBBD HypRel Boussinesq-type equations
            uR=huR/hR
            vR=hvR/hR
            uL=huL/hL
            vL=hvL/hL
            wL = qr(4,i-1)/hL
            wR = ql(4,i)/hR
            sgL = qr(5,i-1)/hL
            sgR = ql(5,i)/hR
            pL = qr(6,i-1)/hL
            pR = ql(6,i)/hR
            pbL = qr(7,i-1)/hL
            pbR = ql(7,i)/hR
            
            call riemann_normal_HypRel_BBBD(mu, nv, g, c_sq, &
                     hl, ul, vl, wl, sgl, pl, pbl, bl, csq_l,&
                     hr, ur, vr, wr, sgr, pr, pbr, br, csq_r, &
                     fwHypRel, swHypRel)

            ! Computing f-waves and speeds
            do mw = 1, mwaves
                fwave(:, mw, i) = fwHypRel(:,mw)
                s(mw,i) = swHypRel(mw)
            END do
         
         else
            print *, "Not shallow water equations or Boussinesq equations"
            call abort
         end if

 30      continue
      enddo


!c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxl(3,i))*deg2rad
          endif

          do mw=1,mwaves
! c             if (s(mw,i) .gt. 316.d0) then
! c               # shouldn't happen unless h > 10 km!
! c                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
! c                endif
               s(mw,i)=dxdc*s(mw,i)
               fwave(:,mw,i)=dxdc*fwave(:,mw,i)
               !fwave(2,mw,i)=dxdc*fwave(2,mw,i)
               !fwave(3,mw,i)=dxdc*fwave(3,mw,i)
          enddo
         enddo
        endif

!===============================================================================


!============= compute fluctuations=============================================
         amdq(:,:) = 0.d0
         apdq(:,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
               if (s(mw,i) <-1.e-16) then
                     amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
               else if (s(mw,i)>1.e-16) then
                  apdq(:,i)  = apdq(:,i) + fwave(:,mw,i)
               else
                 amdq(:,i) = amdq(:,i) + 0.5d0 * fwave(:,mw,i)
                 apdq(:,i) = apdq(:,i) + 0.5d0 * fwave(:,mw,i)
               endif
            enddo
         enddo
!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo

      return
      end subroutine
