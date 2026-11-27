
      subroutine urqmd_cola_init_particles(parts)
         use cola
         implicit none
         
         type(EventParticles), intent(inout) :: parts
         type(Particle) :: one_part
         type(LorentzVector) :: pos, mom
         integer :: i
         integer, external :: getspin, fchg
         include 'coms.f'
         include 'options.f'
         
         

         do i = 1, npart
            one_part = parts%get(i)
            pos = one_part%get_position()
            mom = one_part%get_momentum()
            ! set coords
            r0(i) = acttime
            rx(i) = pos%get_x()
            ry(i) = pos%get_y()
            rz(i) = pos%get_z()
            ! set momentum
            p0(i) = mom%get_e() / 1000.0
            px(i) = (mom%get_x() / 1000.0 - ffermpx(i))
            py(i) = (mom%get_y() / 1000.0 - ffermpy(i))
            pz(i) = (mom%get_z() / 1000.0 - ffermpz(i))
            
            call pdg_to_ityp(one_part%get_pdg_code(), 
     &                    ityp(i), iso3(i))
            fmass(i) = sqrt(mom%Mag2()) / 1000.0
            !spin(i) = getspin(ityp(i), -1)
            charge(i)=fchg(iso3(i),ityp(i))
            
            
         enddo
         return
      end
      
      subroutine urqmd_cola_file14out_to_parts_conv(timestep,
     & parts, coll_to_r0)
      use cola
      use last_event_storage
      implicit none
      integer timestep
      type(EventParticles), intent(inout) :: parts
      type(Particle) :: p
      type(LorentzVector) :: mom, pos
      integer i, pdgid, pdg, pclass
      integer urqmd_cola_classify_pclass
      logical, external :: urqmd_cola_build_classification
      logical :: coll_to_r0
      real*8 p0v, pxv, pyv, pzv, last_coll_time
      include 'coms.f'
      include 'options.f'
      include 'freezeout.f'
      include 'comhis.f'
      if (.not.urqmd_cola_build_classification()) then
         write(6,*) '*** urqmd_cola: build_classification failed'
         return
      endif
      do i = 1, npart
         pdg = pdgid(ityp(i), iso3(i))
         p0v = 1000.0 * p0(i) ! to MeV
         pxv = 1000.0 * (px(i) + ffermpx(i)) ! to MeV
         pyv = 1000.0 * (py(i) + ffermpy(i)) ! to MeV
         pzv = 1000.0 * (pz(i) + ffermpz(i)) ! to MeV
         pclass = urqmd_cola_classify_pclass(uid(i), 
     &         pzv, ncoll(i), lstcoll(i))
         p = Particle()
         call p%set_pdg_code(pdg)
         call p%set_p_class(pclass)
         mom = LorentzVector(p0v, pxv, pyv, pzv)
         call p%set_momentum(mom)

         if (coll_to_r0) then
            if (lstcoll(i) .gt. 0) then
               last_coll_time = hisacttime(lstcoll(i))
            else
               last_coll_time = 0.d0
            endif
            pos = LorentzVector(last_coll_time, rx(i), ry(i), rz(i))
         else
            pos = LorentzVector(r0(i), rx(i), ry(i), rz(i))
         endif
         call p%set_position(pos)
         
         call parts%push_back(p)

      enddo
      call set_evolution_essentials
      return
      end


      subroutine pdg_to_az(pdg_in, A, Z)
      implicit none
      integer, intent(in)  :: pdg_in
      integer, intent(out) :: A, Z
      integer :: k
      A = 0;  Z = 0
      if (pdg_in == 2112) then; A=1; Z=0; return; end if
      if (pdg_in == 2212) then; A=1; Z=1; return; end if
      k = abs(pdg_in) / 10
      Z = mod(k,10) + mod(k/10,10)*10 + mod(k/100,10)*100
      k = k / 1000
      A = mod(k,10) + mod(k/10,10)*10 + mod(k/100,10)*100
      end subroutine

      subroutine urqmd_cola_run_propagation(
     & prop_time, ini, input_parts, 
     & output_parts, coll_to_r0, first_conv)
      use cola
      use last_event_storage
      implicit none
      real*8 prop_time
      type(EventParticles), intent(inout) :: input_parts
      type(EventParticles), intent(inout) :: output_parts
      type(EventIniState), intent(inout) :: ini
      type(AZ) :: AZp, AZt
      logical :: coll_to_r0, first_conv

      include 'coms.f'
      include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
      include 'freezeout.f'
      

      integer i,j,k,steps,ii,ocharge,ncharge, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx,CTOsave
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav
      real*8 thydro_start,thydro,nucrad
      logical lhydro

      npart = input_parts%size()
      call pdg_to_az(ini%get_pdg_code_a(), Ap, Zp)
      call pdg_to_az(ini%get_pdg_code_b(), At, Zt)

      
!     initialize random number generator (per event)
      if (first_conv) then
         if (.not.firstseed .and. (.not.fixedseed)) then
            ranseed = -(1*abs(ranseed))
            call sseed(ranseed)
         else
            firstseed = .false.
         endif

         call init_resrec
         call init_evolution_essentials
      endif
      call get_evolution_essentials
!     Enable collision history (comhis) for file15-like classification.
      CTOption(68) = 1
      bf13 = .false.

!     Disable coalescence
      CTOption(64) = 0

!     Disable final decay
      CTOption(18) = 1

!     Init resonance reconstruction (clears comhis for new event)
      if (.not.bf13.and.CTOption(68).eq.1) then
         call init_resrec
      endif


      if(prop_time.gt.0.0) then
         !dtimestep = prop_time
         dtimestep = prop_time
         nsteps = nint(0.01 + prop_time / dtimestep)
         outsteps = nint(0.01 + prop_time / dtimestep)
      else
         dtimestep = 0.2
         nsteps = 0
         outsteps = 0
      endif

      event = 1
      if (first_conv) then
         time = 0.0
      else 
         time = acttime
      endif
      lhydro = .true.
      

!     eccentricity (init sets nuclei; this computes Glauber observables)
      !call init
      !call init_eccentricity
      call urqmd_cola_init_particles(input_parts)

      if (CTOption(40).ne.0 .and. (.not.success)) return

!     hydro switch
      if (CTOption(45).eq.1) then
         thydro_start = CTParam(65)*2.d0*nucrad(Ap)*sqrt(2.d0*emnuc/ebeam)
         if (thydro_start.lt.CTParam(63)) then
            thydro_start = CTParam(63)
         endif
      endif

      if (CTOption(40).ne.0) time = acttime
!     output preparation

!      call output(13)
!      call output(14)
!      call output(15)
!      call output(16)
!      if (event.eq.1) then
!         call output(17)
!         call osc_header
!         call osc99_header
!      endif
!      call osc99_event(-1)

!     for CTOption(4)=1 : output of initialization configuration
!      if (CTOption(4).eq.1) call file14out(0)
      if (CTOption(4).eq.1) call urqmd_cola_file14out_to_parts_conv(0, output_parts, coll_to_r0)
!     participant/spectator model
      if (CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

      otime = outsteps*dtimestep
      steps = 0
      !if (nsteps.eq.0) then
      !   call urqmd_cola_file14out_to_parts_conv(0, output_parts, coll_to_r0)
      !endif
!     loop over all timesteps

      do steps = 1, nsteps
         if (eos.ne.0) then
            do j = 1, npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
            enddo
         endif
         acttime = time
         if (CTOption(16).eq.0) then
            call colload
            if (nct.gt.0) then
               coll_loop: do
                  k = 0
                  do
                     call getnext(k)
                     if (k.eq.0) exit coll_loop

                     if (CTOption(45).eq.1) then
                        if (cttime(k).gt.thydro_start .and. lhydro) then
                           if (CTOption(62).eq.1) then
                              call prepout
                              call urqmd_cola_file14out_to_parts_conv(0, output_parts, coll_to_r0)
                              call restore
                           endif
                           st = thydro_start - acttime
                           call cascstep(acttime,st)
                           acttime = thydro_start
                           lhydro = .false.
                           if (CTOption(50).eq.1) return
                           if (thydro.gt.1.d-8 .or. CTOption(48).eq.1) then
                              call colload
                              cycle coll_loop
                           endif
                        endif
                     endif

                     st = cttime(k) - acttime
                     call cascstep(acttime,st)
                     acttime = cttime(k)

                     if (cti2(k).gt.0 .and. abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3) then
                        write(6,*) ' ***(E) wrong collision update (col) ***'
                     else if (cti2(k).eq.0 .and. abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
                        write(6,*) ' *** main(W) wrong collision update (decay)'
                     endif

                     ocharge = charge(cti1(k))
                     if (cti2(k).gt.0) ocharge = ocharge + charge(cti2(k))
                     it1 = ityp(cti1(k))
                     if (cti2(k).gt.0) it2 = ityp(cti2(k))

                     cti1sav = cti1(k)
                     cti2sav = cti2(k)
                     call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k), ctcolfluc(k))

                     if (CTOption(17).eq.0) then
                        if (nexit.eq.0) then
                           if (cti1(k).ne.cti1sav .or. cti2(k).ne.cti2sav) then
                              cti1(k) = cti1sav
                              cti2(k) = cti2sav
                           endif
                           call collupd(cti1(k),1)
                           if (cti2(k).gt.0) call collupd(cti2(k),1)
                        else
                           ncharge = 0
                           do i = 1, nexit
                              ncharge = ncharge + charge(inew(i))
                              call collupd(inew(i),1)
                           enddo
                           do ii = 1, nsav
                              call collupd(ctsav(ii),1)
                           enddo
                           nsav = 0
                        endif
                     else
                        call colload
                     endif

                     if (CTOption(17).eq.0) cycle
                     exit
                  enddo
               enddo coll_loop
            endif
         endif

         time = time + dtimestep
         call cascstep(acttime,time-acttime)

         if (eos.ne.0) then
            do j = 1, npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
            enddo
            call proprk(time,dtimestep)
         endif

!     perform output if desired
         if (mod(steps,outsteps).eq.0 .and. steps.lt.nsteps) then
            if (CTOption(28).eq.2) call spectrans(otime)
            if (CTOption(62).eq.1) call prepout
!            call file14out(steps)
            call urqmd_cola_file14out_to_parts_conv(steps, output_parts, coll_to_r0)
!           if (CTOption(64).eq.1) call file13out(steps)
            if (CTOption(62).eq.1) then
               call restore
               call colload
            endif
!           if (CTOption(55).eq.1) call osc_vis(steps)
         endif
      enddo

!     final decay of unstable particles
      acttime = time
      if (CTOption(18).eq.0) then
         nct = 0
         actcol = 0
         CTOsave = CTOption(10)
         CTOption(10) = 1
         do i = 1, npart
            do while (dectime(i).lt.1.d30)
               isstable = .false.
               do stidx = 1, nstable
                  if (ityp(i).eq.stabvec(stidx)) isstable = .true.
               enddo
               if (isstable) exit
               call scatter(i,0,0.d0,fmass(i),xdummy)
            enddo
         enddo
         CTOption(10) = CTOsave
      endif
      if (CTOption(64).eq.1) call coalescence

!     Match pure UrQMD: file13out, file14out, file16out (same execution order)
!      call file13out(nsteps)
      if (CTOption(50).eq.0) then
!       call file14out(nsteps)
         call urqmd_cola_file14out_to_parts_conv(nsteps, output_parts, coll_to_r0)
      endif
!      call file16out
!      if (CTOption(50).eq.0.and.CTOption(55).eq.0) call osc_event
!      if (CTOption(50).eq.0.and.CTOption(55).eq.1) call osc_vis(nsteps)
!      call osc99_event(1)
!      call osc99_eoe

      return
      end
