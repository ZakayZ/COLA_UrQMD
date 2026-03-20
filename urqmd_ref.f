! ---------------------------------------------------------------------------
! Glue subroutines callable from cola_fortran_generator_impl.
! ---------------------------------------------------------------------------
      subroutine urqmd_cola_uinit
      implicit none

      call uinit(0)
      return
      end

      subroutine urqmd_cola_fill_inistate_particles(ini_parts, pza, pzb)
      use cola
      implicit none
      type(EventParticles), intent(inout) :: ini_parts
      real*8, intent(out) :: pza, pzb
      type(Particle) :: p
      type(LorentzVector) :: mom, pos
      integer i, pdg, pdgid
      include 'inputs.f'

      pza = 0.d0
      pzb = 0.d0
      do i = 1, PT_AA(1)
         pdg = pdgid(PT_ityp(i,1), PT_iso3(i,1))
         p = Particle()
         call p%set_pdgCode(pdg)
         call p%set_pClass(ParticleClass_SPECTATOR_A)
         mom = LorentzVector(PT_p0(i,1), PT_px(i,1), PT_py(i,1), PT_pz(i,1))
         call p%set_momentum(mom)
         pos = LorentzVector(PT_r0(i,1), PT_rx(i,1), PT_ry(i,1), PT_rz(i,1))
         call p%set_position(pos)
         call ini_parts%push_back(p)
         pza = pza + PT_pz(i,1)
      enddo
      do i = 1, PT_AA(2)
         pdg = pdgid(PT_ityp(i,2), PT_iso3(i,2))
         p = Particle()
         call p%set_pdgCode(pdg)
         call p%set_pClass(ParticleClass_SPECTATOR_B)
         mom = LorentzVector(PT_p0(i,2), PT_px(i,2), PT_py(i,2), PT_pz(i,2))
         call p%set_momentum(mom)
         pos = LorentzVector(PT_r0(i,2), PT_rx(i,2), PT_ry(i,2), PT_rz(i,2))
         call p%set_position(pos)
         call ini_parts%push_back(p)
         pzb = pzb + PT_pz(i,2)
      enddo
      return
      end

      subroutine urqmd_cola_disable_outputs
      implicit none
      include 'options.f'
!     bf13=false: f15outch populates comhis only when bf13
!     is false; needed for file15-like particle classification
      bf13 = .false.
      bf14 = .true.
      bf15 = .true.
      bf16 = .true.
      bf19 = .true.
      bf20 = .true.
      return
      end

      subroutine urqmd_cola_generate_tables(tabpath, ok)
      implicit none
      character*(*) tabpath
      logical ok
      integer ios
      character*77 tname

      tname = tabpath
      if (tname(1:4).eq.'    ') tname = 'tables.dat'

!     If an old file exists, remove it so uinit/loadwtab regenerates it.
      open (unit=75, iostat=ios, file=tname, form='unformatted', status='old')
      if (ios.eq.0) then
         close (unit=75, status='delete')
      endif

      call uinit(0)

      open (unit=75, iostat=ios, file=tname, form='unformatted', status='old')
      if (ios.eq.0) then
         close (unit=75, status='keep')
         ok = .true.
      else
         ok = .false.
      endif
      return
      end

      logical function is_elastic_proc(proc_id)
!     Elastic process IDs from UrQMD Table 12:
!     13=BB elastic (not pp/pn), 17=pn elastic, 19=pp elastic,
!     22=BB elastic, 26=MB elastic, 38=MM elastic
      implicit none
      integer proc_id
      is_elastic_proc = (proc_id.eq.13 .or. proc_id.eq.17 .or. proc_id.eq.19 .or.
     &     proc_id.eq.22 .or. proc_id.eq.26 .or. proc_id.eq.38)
      return
      end

! ---------------------------------------------------------------------------
! Particle classification into COLA ParticleClass
! ---------------------------------------------------------------------------
      logical function urqmd_cola_build_classification()
      use cola
      implicit none
      include 'coms.f'
      include 'comhis.f'
      include 'options.f'
      integer iww, i1, ii, iuid, nprocs
      integer uid_nucleus_a, uid_nucleus_b
      parameter (uid_nucleus_a = 1, uid_nucleus_b = 2)
      integer uid_initial(nmax), uid_produced(nmax)
      integer uid_nprocs(nmax), uid_procs(nmax, 32)
      integer proc_id
      common /urqmd_cola_f15/ uid_initial, uid_produced, uid_nprocs, uid_procs
      save /urqmd_cola_f15/

!     Parameter consistency checks (must match comhis.f)
      if (ctag.gt.wwmax) then
         write(6,*) '*** urqmd_cola: ctag.gt.wwmax, increase wwmax'
         urqmd_cola_build_classification = .false.
         return
      endif
      if (uid_cnt.gt.nmax) then
         write(6,*) '*** urqmd_cola: uid_cnt.gt.nmax, increase nmax'
         urqmd_cola_build_classification = .false.
         return
      endif

      do i1 = 1, min(uid_cnt, nmax)
         uid_initial(i1) = 0
         uid_produced(i1) = 0
         uid_nprocs(i1) = 0
      enddo

      do iww = 1, ctag
         if (hisiline(iww).eq.0) cycle
         if (hisnin(iww).gt.inmax .or. hisnexit(iww).gt.outmax) then
            write(6,*) '*** urqmd_cola: hisnin/hisnexit exceed limit'
            urqmd_cola_build_classification = .false.
            return
         endif
         proc_id = mod(iabs(hisiline(iww)), 100)

!        Process incoming: first-seen uids are from initial nucleus
         do i1 = 1, hisnin(iww)
            iuid = INind(iww, i1)
            if (iuid.le.0 .or. iuid.gt.nmax) cycle
            if (uid_produced(iuid).eq.0 .and. uid_initial(iuid).eq.0) then
               if (INpz(iww, i1).ge.0.d0) then
                  uid_initial(iuid) = uid_nucleus_a
               else
                  uid_initial(iuid) = uid_nucleus_b
               endif
            endif
            nprocs = uid_nprocs(iuid)
            if (nprocs.lt.32) then
               nprocs = nprocs + 1
               uid_procs(iuid, nprocs) = proc_id
               uid_nprocs(iuid) = nprocs
            endif
         enddo

!        Process outgoing: first-seen uids (not in initial) are produced
         do ii = 1, hisnexit(iww)
            iuid = OUTind(iww, ii)
            if (iuid.le.0 .or. iuid.gt.nmax) cycle
            if (uid_initial(iuid).eq.0) uid_produced(iuid) = 1
            nprocs = uid_nprocs(iuid)
            if (nprocs.lt.32) then
               nprocs = nprocs + 1
               uid_procs(iuid, nprocs) = proc_id
               uid_nprocs(iuid) = nprocs
            endif
         enddo
      enddo
      urqmd_cola_build_classification = .true.
      return
      end

      integer function urqmd_cola_classify_pclass(uidv, pzv, ncollv, lstcollv)
      use cola
      implicit none
      include 'coms.f'
      integer uidv, ncollv, lstcollv
      real*8 pzv
      integer uid_nucleus_a, uid_nucleus_b
      parameter (uid_nucleus_a = 1, uid_nucleus_b = 2)
      integer uid_initial(nmax), uid_produced(nmax)
      integer uid_nprocs(nmax), uid_procs(nmax, 32)
      common /urqmd_cola_f15/ uid_initial, uid_produced, uid_nprocs, uid_procs
      integer ip
      logical all_elastic
      logical, external :: is_elastic_proc

      ! coalesence
      if (lstcollv.lt.-1 .and. ncollv.eq.0) then
         urqmd_cola_classify_pclass = ParticleClass_PRODUCED
         return
      endif

      ! out of bounds
      if (uidv.le.0 .or. uidv.gt.nmax) then
         urqmd_cola_classify_pclass = ParticleClass_PRODUCED
         return
      endif

      ! marked produced
      if (uid_produced(uidv).eq.1) then
         urqmd_cola_classify_pclass = ParticleClass_PRODUCED
         return
      endif

      ! wasn't initiated for some reason
      if (uid_initial(uidv).eq.0) then
         urqmd_cola_classify_pclass = ParticleClass_PRODUCED
         return
      endif

      ! spectator (no interactions)
      if (uid_nprocs(uidv).eq.0) then
         if (uid_initial(uidv).eq.uid_nucleus_a) then
            urqmd_cola_classify_pclass = ParticleClass_SPECTATOR_A
         else
            urqmd_cola_classify_pclass = ParticleClass_SPECTATOR_B
         endif
         return
      endif

      ! elastic (if all interactions are elastic)
      all_elastic = .true.
      do ip = 1, uid_nprocs(uidv)
         if (.not.is_elastic_proc(uid_procs(uidv, ip))) all_elastic = .false.
      enddo

      if (all_elastic) then
         if (uid_initial(uidv).eq.uid_nucleus_a) then
            urqmd_cola_classify_pclass = ParticleClass_ELASTIC_A
         else
            urqmd_cola_classify_pclass = ParticleClass_ELASTIC_B
         endif
      else
         if (uid_initial(uidv).eq.uid_nucleus_a) then
            urqmd_cola_classify_pclass = ParticleClass_NONELASTIC_A
         else
            urqmd_cola_classify_pclass = ParticleClass_NONELASTIC_B
         endif
      endif
      return
      end

! ---------------------------------------------------------------------------
! Mimic file14out but write particle records to COLA EventParticles.
! ---------------------------------------------------------------------------
      subroutine urqmd_cola_file14out_to_parts(timestep, parts)
      use cola
      implicit none
      integer timestep
      type(EventParticles), intent(inout) :: parts
      type(Particle) :: p
      type(LorentzVector) :: mom, pos
      integer i, pdgid, pdg, pclass
      integer urqmd_cola_classify_pclass
      logical, external :: urqmd_cola_build_classification
      real*8 p0v, pxv, pyv, pzv
      include 'coms.f'
      include 'options.f'

      if (.not.urqmd_cola_build_classification()) then
         write(6,*) '*** urqmd_cola: build_classification failed'
         return
      endif

      do i = 1, npart
         pdg = pdgid(ityp(i), iso3(i))
         p0v = p0(i)
         pxv = px(i) + ffermpx(i)
         pyv = py(i) + ffermpy(i)
         pzv = pz(i) + ffermpz(i)
         pclass = urqmd_cola_classify_pclass(uid(i), pzv, ncoll(i), lstcoll(i))
         p = Particle()
         call p%set_pdgCode(pdg)
         call p%set_pClass(pclass)
         mom = LorentzVector(p0v, pxv, pyv, pzv)
         call p%set_momentum(mom)
         pos = LorentzVector(r0(i), rx(i), ry(i), rz(i))
         call p%set_position(pos)
         call parts%push_back(p)
      enddo
      return
      end

! ---------------------------------------------------------------------------
! One-event UrQMD loop copied from urqmd.f and adapted:
! file writes commented out
! all particles are going to cola particle vector
! ---------------------------------------------------------------------------
      subroutine urqmd_cola_run_one_event(ebeam_out, bimp_out, np, parts, ini_parts, pdga, pdgb, pza, pzb, eout, snnout,
     &     nco, npp, npn, nnn, npt, npa, nbt,
     &     phia, thaa, phib, thab)
      use cola
      implicit none
      real*8 ebeam_out, bimp_out
      integer np
      type(EventParticles), intent(inout) :: parts
      type(EventParticles), intent(inout) :: ini_parts
      integer, intent(out) :: pdga, pdgb
      real*8, intent(out) :: pza, pzb, eout, snnout
      integer, intent(out) :: nco, npp, npn, nnn
      integer, intent(out) :: npt, npa, nbt
      real*8, intent(out) :: phia, thaa, phib, thab
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
      real*8 sigtot

      parts = EventParticles()
      ini_parts = EventParticles()
      pdga = 0
      pdgb = 0
      pza = 0.d0
      pzb = 0.d0
      eout = 0.d0
      snnout = 0.d0
      nco = 0
      npp = 0
      npn = 0
      nnn = 0
      npt = 0
      npa = 0
      nbt = 0
      phia = 0.d0
      thaa = 0.d0
      phib = 0.d0
      thab = 0.d0

      event = 1
      time = 0.0
      lhydro = .true.

!     initialize random number generator (per event)
      if (.not.firstseed .and. (.not.fixedseed)) then
         ranseed = -(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed = .false.
      endif

!     Enable collision history (comhis) for file15-like classification.
      CTOption(68) = 1
      bf13 = .false.

!     Init resonance reconstruction (clears comhis for new event)
      if (.not.bf13.and.CTOption(68).eq.1) then
       call init_resrec
      endif

!     eccentricity (init sets nuclei; this computes Glauber observables)
      call init
      call init_eccentricity
      pdga = AZToPdg(MakeAZ(Ap, Zp))
      pdgb = AZToPdg(MakeAZ(At, Zt))
      call urqmd_cola_fill_inistate_particles(ini_parts, pza, pzb)
      if (abs(pzb).lt.1.d-12) then
         eout = ebeam
      else
         eout = ecm
      endif
      snnout = sigtot(1,2,ecm)

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
      if (CTOption(4).eq.1) call urqmd_cola_file14out_to_parts(0, parts)

!     participant/spectator model
      if (CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

      otime = outsteps*dtimestep
      steps = 0

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
                              call urqmd_cola_file14out_to_parts(0, parts)
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
            call urqmd_cola_file14out_to_parts(steps, parts)
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
         call urqmd_cola_file14out_to_parts(nsteps, parts)
      endif
!      call file16out
!      if (CTOption(50).eq.0.and.CTOption(55).eq.0) call osc_event
!      if (CTOption(50).eq.0.and.CTOption(55).eq.1) call osc_vis(nsteps)
!      call osc99_event(1)
!      call osc99_eoe

      ebeam_out = ebeam
      bimp_out = bimp
      np = npart
      nco = ctag - dectag
      do i = 1, Ap
         if (ncoll(i).gt.0) npa = npa + 1
      enddo
      do i = Ap + 1, Ap + At
         if (ncoll(i).gt.0) nbt = nbt + 1
      enddo
      npt = npa + nbt
      return
      end
