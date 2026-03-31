module cola_fortran_generator_impl
  use cola
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(AbstractFortranGenerator) :: URQMDGenerator
    character(len=512) :: input_file = ''
    character(len=512) :: generated_config_file = ''
    character(len=512) :: tables_file = ''
  contains
    procedure :: init => generator_init
    procedure :: run => generator_run
    final :: generator_final
  end type URQMDGenerator

  interface
    function setenv(name, value, overwrite) bind(c)
      import :: c_int, c_char
      integer(c_int) :: setenv
      character(kind=c_char), intent(in) :: name(*), value(*)
      integer(c_int), value :: overwrite
    end function
  end interface

  interface
    subroutine urqmd_cola_uinit()
    end subroutine

    subroutine urqmd_cola_disable_outputs()
    end subroutine

    subroutine urqmd_cola_generate_tables(tabpath, ok)
      character(len=*), intent(in) :: tabpath
      logical, intent(out) :: ok
    end subroutine

    subroutine urqmd_cola_run_one_event(ebeam_out, bimp_out, np, parts, ini_parts, &
                                        pdga, pdgb, pza, pzb, eout, snnout, nco, npp, &
                                        npn, nnn, npt, npa, nbt, phia, thaa, phib, thab)
      import :: EventParticles
      real(8), intent(out) :: ebeam_out, bimp_out
      integer, intent(out) :: np
      type(EventParticles), intent(inout) :: parts
      type(EventParticles), intent(inout) :: ini_parts
      integer, intent(out) :: pdga, pdgb
      real(8), intent(out) :: pza, pzb, eout, snnout
      integer, intent(out) :: nco, npp, npn, nnn, npt, npa, nbt
      real(8), intent(out) :: phia, thaa, phib, thab
    end subroutine
  end interface

contains

  subroutine urqmd_cola_set_env(name_no_nul, path)
    character(len=*), intent(in) :: name_no_nul, path
    character(len=64) :: name
    character(len=1024) :: path_c
    integer(c_int) :: ierr
    name = trim(name_no_nul) // c_null_char
    path_c = trim(path) // c_null_char
    ierr = setenv(name, path_c, 1_c_int)
  end subroutine urqmd_cola_set_env

  logical function is_control_key(key)
    character(len=*), intent(in) :: key
    is_control_key = (key == 'config_file' .or. key == 'generated_config_file' .or. key == 'tables_file')
  end function is_control_key

  subroutine get_tmpdir(tmpdir)
    character(len=512), intent(out) :: tmpdir
    integer :: len, status
    tmpdir = ''
    call get_environment_variable('TMPDIR', tmpdir, len, status)
    if (len_trim(tmpdir) == 0) tmpdir = '/tmp'
  end subroutine get_tmpdir

  subroutine resolve_tables_path(self, path, exists)
    class(URQMDGenerator), intent(in) :: self
    character(len=512), intent(out) :: path
    logical, intent(out) :: exists
    path = ''
    exists = .false.
    if (len_trim(self%tables_file) > 0) then
      path = trim(self%tables_file)
      inquire(file=trim(path), exist=exists)
      if (exists) return
    end if
    path = 'tables.dat'
    inquire(file=trim(path), exist=exists)
  end subroutine resolve_tables_path

  subroutine generate_tables_external(table_path, ok)
    character(len=*), intent(in) :: table_path
    logical, intent(out) :: ok
    call urqmd_cola_set_env('URQMD_TAB', trim(table_path))
    call urqmd_cola_generate_tables(trim(table_path), ok)
  end subroutine generate_tables_external

  subroutine generator_init(self, pmap, err)
    class(URQMDGenerator), intent(inout) :: self
    type(ParametersMap), intent(in) :: pmap
    character(len=:), allocatable, intent(out) :: err
    type(ParametersMapItem) :: kv
    character(len=:), allocatable :: key, val
    character(len=512) :: tmpdir
    character(len=512) :: tables_path
    integer :: i, n, iostat, u
    logical :: tables_ok

    self%input_file = ''
    self%generated_config_file = ''
    self%tables_file = ''
    err = ''
    n = pmap%size()

    do i = 1, n
      kv = pmap%get(i)
      key = kv%get_first()
      val = kv%get_second()

      if (key == 'config_file') then
        self%input_file = trim(val)
      else if (key == 'generated_config_file') then
        self%generated_config_file = trim(val)
      else if (key == 'tables_file') then
        self%tables_file = trim(val)
      end if
    end do

    if (len_trim(self%input_file) == 0) then
      if (len_trim(self%generated_config_file) == 0) then
        call get_tmpdir(tmpdir)
        tmpdir = trim(tmpdir)
        self%generated_config_file = trim(tmpdir) // '/urqmd_cola_config.txt'
      end if

      ! Generate UrQMD config file
      open(newunit=u, file=trim(self%generated_config_file), status='replace', action='write', iostat=iostat)
      if (iostat /= 0) then
        err = 'URQMD init failed: failed to create generated config file.'
        return
      end if
      do i = 1, n
        kv = pmap%get(i)
        key = kv%get_first()
        val = kv%get_second()
        if (is_control_key(key)) cycle
        write(u, '(a)') trim(key) // ' ' // trim(val)
      end do
      write(u, '(a)') 'xxx'
      close(u)
      self%input_file = trim(self%generated_config_file)
    end if

    call urqmd_cola_set_env('ftn09', trim(self%input_file))

    ! Stop UrQMD from writing to files
    call urqmd_cola_set_env('ftn13', ' ')
    call urqmd_cola_set_env('ftn14', ' ')
    call urqmd_cola_set_env('ftn15', ' ')
    call urqmd_cola_set_env('ftn16', ' ')
    call urqmd_cola_set_env('ftn19', ' ')
    call urqmd_cola_set_env('ftn20', ' ')

    call resolve_tables_path(self, tables_path, tables_ok)
    if (.not. tables_ok) then
      call generate_tables_external(trim(tables_path), tables_ok)
      if (tables_ok) then
        call resolve_tables_path(self, tables_path, tables_ok)
      end if
    end if
    if (.not. tables_ok) then
      err = 'URQMD init failed: failed to generate tables.dat.'
      return
    end if
    call urqmd_cola_set_env('URQMD_TAB', trim(tables_path))
    call urqmd_cola_set_env('ftn09', trim(self%input_file))

    call urqmd_cola_uinit()
    call urqmd_cola_disable_outputs()
  end subroutine generator_init

  function generator_run(self, err) result(ed)
    class(URQMDGenerator), intent(in) :: self
    character(len=:), allocatable, intent(out) :: err
    type(EventData) :: ed
    type(EventIniState) :: ini
    type(EventParticles) :: parts, ini_parts
    integer :: np, pdga, pdgb
    integer :: nco, npp, npn, nnn
    integer :: npt, npa, nbt
    real(8) :: ebeam_val, bimp_val, pza, pzb, eout, snnout
    real(8) :: phia, thaa, phib, thab

    ed = EventData()
    err = ''
    parts = EventParticles()
    ini_parts = EventParticles()

    call urqmd_cola_run_one_event(ebeam_val, bimp_val, np, parts, ini_parts, &
                                  pdga, pdgb, pza, pzb, eout, snnout, nco, npp, &
                                  npn, nnn, npt, npa, nbt, phia, thaa, phib, thab)

    ini = ed%get_iniState()
    call ini%set_pdgCodeA(pdga)
    call ini%set_pdgCodeB(pdgb)
    call ini%set_pZA(pza)
    call ini%set_pZB(pzb)
    call ini%set_energy(eout)
    call ini%set_sectNN(real(snnout, kind(0.0)))
    call ini%set_b(real(bimp_val, kind(0.0)))
    call ini%set_nColl(nco)
    call ini%set_nCollPP(npp)
    call ini%set_nCollPN(npn)
    call ini%set_nCollNN(nnn)
    call ini%set_nPart(npt)
    call ini%set_nPartA(npa)
    call ini%set_nPartB(nbt)
    call ini%set_phiRotA(real(phia, kind(0.0)))
    call ini%set_thetaRotA(real(thaa, kind(0.0)))
    call ini%set_phiRotB(real(phib, kind(0.0)))
    call ini%set_thetaRotB(real(thab, kind(0.0)))
    call ini%set_iniStateParticles(ini_parts)
    if (np < 0) np = 0
    call ed%set_iniState(ini)
    call ed%set_particles(parts)

  end function generator_run

  subroutine generator_final(self)
    type(URQMDGenerator), intent(inout) :: self
    integer :: u, ios
    logical :: fexists

    if (len_trim(self%generated_config_file) == 0) return

    inquire(file=trim(self%generated_config_file), exist=fexists)
    if (.not. fexists) return

    open(newunit=u, file=trim(self%generated_config_file), status='old', iostat=ios)
    if (ios /= 0) return
    close(u, status='delete', iostat=ios)
  end subroutine generator_final
end module cola_fortran_generator_impl
