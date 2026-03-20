! Symlink installed (or source-tree) eosfiles into the process cwd so legacy
! OPEN(..., file='eosfiles/...') paths resolve.
subroutine cola_ensure_eosfiles_link()
  implicit none
  include 'eosfiles_path.inc'
  character(len=2048) :: cmd
  integer :: istat
  logical :: exists

  if (.not. cola_eosfiles_auto) return

  inquire (file='eosfiles', exist=exists)
  if (exists) return

  cmd = 'ln -sf "' // trim(cola_eosfiles_install_dir) // '" eosfiles'
  call execute_command_line(trim(cmd), exitstat=istat)
end subroutine cola_ensure_eosfiles_link
