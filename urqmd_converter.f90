
module cola_fortran_converter_impl
  use cola
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public, extends(AbstractFortranConverter) :: URQMDConverter
    real(8)            :: propagate_time = 0.0 
    logical            :: pass_coll = .true.  
    logical            :: first_conv = .false.  
    character(len=512) :: input_file = ''
    character(len=512) :: generated_config_file = ''
    character(len=512) :: tables_file = ''
  contains
    procedure :: init => converter_init
    procedure :: run => converter_run
    final :: converter_final
  end type URQMDConverter

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

    subroutine urqmd_cola_run_propagation(prop_time, ini, parts_in, parts_out, coll, first)
      import :: EventParticles
      import :: EventIniState
      real(8), intent(in) :: prop_time
      logical, intent(in) :: coll, first
      
      type(EventParticles), intent(inout) :: parts_in
      type(EventParticles), intent(inout) :: parts_out
      type(EventIniState), intent(inout) :: ini
       
    end subroutine
  end interface

  contains

  !=============SETUP================
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
    class(URQMDConverter), intent(in) :: self
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
  !=============SETUP================


  subroutine converter_init(self, pmap, err)
    class(URQMDConverter), intent(inout) :: self
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
      else if (key == 'prop_time') then
        read(val, *, iostat=iostat) self%propagate_time
        if (iostat /= 0) then
          print *, "Error in prop_time reading"
        endif
      else if (key == 'first_conv') then
        read(val, *, iostat=iostat) self%first_conv
        if (iostat /= 0) then
          print *, "Error in first_conv reading"
        endif
      else if (key == 'coll_to_r0') then
        read(val, *, iostat=iostat) self%pass_coll
        if (iostat /= 0) then
          print *, "Error in coll_to_r0 reading"
        endif
      endif
    end do

    
    if (.not.self%first_conv) return
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
  end subroutine converter_init


  
  subroutine converter_run(self, ed, err)
    class(URQMDConverter), intent(in) :: self
    type(EventData), intent(inout) :: ed
    character(len=:), allocatable, intent(out) :: err
    type(EventParticles) :: parts_in
    type(EventParticles) :: parts_out
    type(EventIniState) :: ini

    parts_in = ed%get_particles()
    ini = ed%get_ini_state()
    err = ''    
    parts_out = EventParticles()
    call urqmd_cola_run_propagation(self%propagate_time, ini, parts_in, parts_out, self%pass_coll, self%first_conv)
    
    call ed%set_particles(parts_out)


  end subroutine converter_run

  subroutine converter_final(self)
    type(URQMDConverter), intent(inout) :: self
  end subroutine converter_final

end module cola_fortran_converter_impl