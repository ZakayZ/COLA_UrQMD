
        module last_event_storage
        implicit none
        include 'coms.f'
        include 'comres.f'
        include 'comhis.f'
        include 'colltab.f'
        include 'freezeout.f'
        include 'options.f'
        private

        public :: set_evolution_essentials,
     &          get_evolution_essentials,
     &          init_evolution_essentials

        type :: event_vars
            integer :: npart=0, nbar=0, nmes=0, uid_cnt=0,
     &              ctag=0, dectag=0, NHardRes=0, NSoftRes=0,
     &              NDecRes=0, NBlColl=0, NElColl=0,
     &              ranseed

            logical :: firstseed, fixedseed

            real*8 :: dectime(nmax), tform(nmax), xtotfac(nmax),
     &             ffermpx(nmax), ffermpy(nmax), ffermpz(nmax),
     &             hisacttime(wwmax), acttime

            integer :: ncoll(nmax), lstcoll(nmax), uid(nmax)

            integer :: eos=0
            integer :: nstable=0, stabvec(maxstables)

            !---- for build_classification ----
            integer :: hisnin(wwmax), hisnexit(wwmax),
     &              hisiline(wwmax)
            integer :: INind(wwmax,inmax)
            real*8  :: INpz(wwmax,inmax)
            integer :: OUTind(wwmax,outmax)
        end type event_vars

        type(event_vars) :: instance_
        !$omp threadprivate(instance_)

        contains

        subroutine set_evolution_essentials()
            instance_%npart     = npart
            instance_%nbar      = nbar
            instance_%nmes      = nmes
            instance_%uid_cnt   = uid_cnt
            instance_%ctag      = ctag
            instance_%dectag    = dectag
            instance_%NHardRes  = NHardRes
            instance_%NSoftRes  = NSoftRes
            instance_%NDecRes   = NDecRes
            instance_%NBlColl   = NBlColl
            instance_%NElColl   = NElColl
            instance_%acttime   = acttime

            instance_%ncoll(1:npart)   = ncoll(1:npart)
            instance_%lstcoll(1:npart) = lstcoll(1:npart)
            instance_%uid(1:npart)     = uid(1:npart)
            instance_%dectime(1:npart) = dectime(1:npart)
            instance_%tform(1:npart)   = tform(1:npart)
            instance_%xtotfac(1:npart) = xtotfac(1:npart)
            instance_%ffermpx(1:npart) = ffermpx(1:npart)
            instance_%ffermpy(1:npart) = ffermpy(1:npart)
            instance_%ffermpz(1:npart) = ffermpz(1:npart)

            instance_%firstseed = firstseed
            instance_%fixedseed = fixedseed
            instance_%ranseed   = ranseed

            instance_%eos        = eos
            instance_%nstable    = nstable
            if (nstable > 0) instance_%stabvec(1:nstable) = stabvec(1:nstable)

            instance_%hisiline(1:ctag) = hisiline(1:ctag)
            instance_%hisnin(1:ctag)   = hisnin(1:ctag)
            instance_%hisnexit(1:ctag) = hisnexit(1:ctag)
            instance_%INind(1:ctag,:)  = INind(1:ctag,:)
            instance_%INpz(1:ctag,:)   = INpz(1:ctag,:)
            instance_%OUTind(1:ctag,:) = OUTind(1:ctag,:)
            instance_%hisacttime = hisacttime

        end subroutine set_evolution_essentials

        subroutine get_evolution_essentials()
            npart     = instance_%npart
            nbar      = instance_%nbar
            nmes      = instance_%nmes
            uid_cnt   = instance_%uid_cnt
            ctag      = instance_%ctag
            dectag    = instance_%dectag
            NHardRes  = instance_%NHardRes
            NSoftRes  = instance_%NSoftRes
            NDecRes   = instance_%NDecRes
            NBlColl   = instance_%NBlColl
            NElColl   = instance_%NElColl
            acttime   = instance_%acttime

            ncoll(1:npart)   = instance_%ncoll(1:npart)
            lstcoll(1:npart) = instance_%lstcoll(1:npart)
            uid(1:npart)     = instance_%uid(1:npart)
            dectime(1:npart) = instance_%dectime(1:npart)
            tform(1:npart)   = instance_%tform(1:npart)
            xtotfac(1:npart) = instance_%xtotfac(1:npart)
            ffermpx(1:npart) = instance_%ffermpx(1:npart)
            ffermpy(1:npart) = instance_%ffermpy(1:npart)
            ffermpz(1:npart) = instance_%ffermpz(1:npart)

            firstseed = instance_%firstseed
            fixedseed = instance_%fixedseed
            ranseed   = instance_%ranseed

            eos        = instance_%eos
            nstable    = instance_%nstable
            if (nstable > 0) stabvec(1:nstable) = instance_%stabvec(1:nstable)

            hisiline(1:ctag) = instance_%hisiline(1:ctag)
            hisnin(1:ctag)   = instance_%hisnin(1:ctag)
            hisnexit(1:ctag) = instance_%hisnexit(1:ctag)
            INind(1:ctag,:)  = instance_%INind(1:ctag,:)
            INpz(1:ctag,:)   = instance_%INpz(1:ctag,:)
            OUTind(1:ctag,:) = instance_%OUTind(1:ctag,:)
            hisacttime = instance_%hisacttime

        end subroutine get_evolution_essentials

        subroutine init_evolution_essentials()
            instance_%npart     = 0
            instance_%nbar      = 0
            instance_%nmes      = 0
            instance_%uid_cnt   = 0
            instance_%ctag      = 0
            instance_%dectag    = 0
            instance_%NHardRes  = 0
            instance_%NSoftRes  = 0
            instance_%NDecRes   = 0
            instance_%NBlColl   = 0
            instance_%NElColl   = 0
            instance_%acttime   = 0.0d0

            instance_%ncoll     = 0
            instance_%lstcoll   = 0
            instance_%uid       = 0

            instance_%firstseed = .true.
            instance_%fixedseed = .true.
            instance_%ranseed   = 0

            instance_%dectime   = 0.0d0
            instance_%tform     = 0.0d0
            instance_%xtotfac   = 0.0d0
            instance_%ffermpx   = 0.0d0
            instance_%ffermpy   = 0.0d0
            instance_%ffermpz   = 0.0d0
            instance_%hisacttime = 0.0d0

            instance_%eos        = 0
            instance_%nstable    = 0
            instance_%stabvec    = 0
            instance_%hisiline = 0
            instance_%hisnin   = 0
            instance_%hisnexit = 0
            instance_%INind    = 0
            instance_%INpz     = 0.0d0
            instance_%OUTind   = 0

        end subroutine init_evolution_essentials

        end module last_event_storage