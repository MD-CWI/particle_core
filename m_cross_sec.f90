! Copyright 2005-2012, Chao Li, Anbang Sun, Jannis Teunissen
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that contains routines to read in cross section data from textfiles.
module m_cross_sec

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: tiny_len = 20
  integer, parameter :: name_len = 40
  integer, parameter :: line_len = 200

  type CS_coll_t
     integer  :: type       = -1
     real(dp) :: part_mass  = 0
     real(dp) :: rel_mass   = 0
     real(dp) :: en_loss    = 0
     ! Flag to optionally indicate the neutral molecule
     integer  :: gas_index  = 0
  end type CS_coll_t

  !> The type of cross section table
  type CS_t
     type(CS_coll_t) :: coll
     real(dp), allocatable   :: en_cs(:,:)    ! Stores the energy vs cross sec table
     integer                 :: n_rows        ! Number of rows in the table
     real(dp)                :: min_energy    ! Minimum energy for the collision
     real(dp)                :: max_energy    ! Maximum energy for the collision
     character(LEN=tiny_len) :: gas_name      ! Name of the colliding neutral molecule
     character(LEN=line_len) :: description   ! Description of the collision
     character(LEN=line_len) :: comment       ! Additional comments
  end type CS_t

  integer, parameter, public :: CS_attach_t  = 1
  integer, parameter, public :: CS_elastic_t = 2
  integer, parameter, public :: CS_excite_t  = 3
  integer, parameter, public :: CS_ionize_t  = 4
  integer, parameter, public :: CS_effective_t = 5
  integer, parameter, public :: CS_num_types = 5

  ! Index indiciting a collision with neutral molecule CH4
  integer, parameter, public :: CS_is_CH4 = 1

  !> Maximum number of cross sections per gas
  integer, parameter :: max_processes_per_gas = 200

  !> Maximum number of rows per cross section in input files
  integer, parameter :: max_num_rows = 2000

  ! Public variables
  public :: CS_t
  public :: CS_coll_t

  ! Methods
  public :: CS_add_from_file
  public :: CS_write_summary
  public :: CS_create_ledger
  public :: CS_write_ledger

  ! Types of behaviours for dealing with data outside of the input data range
  integer, parameter, public :: CS_extrapolate_error = 0
  integer, parameter, public :: CS_extrapolate_constant = 1
  integer, parameter, public :: CS_extrapolate_zero = 2
  integer, parameter, public :: CS_extrapolate_linear = 3

  ! Types of behaviours for dealing with EFFECTIVE cross sections
  integer, parameter, public :: CS_effective_error = 0
  integer, parameter, public :: CS_effective_auto = 1
  integer, parameter, public :: CS_effective_keep = 2

contains

  !> Read cross section data for a gas from a file
  subroutine CS_add_from_file(filename, gas_name, number_dens, &
       req_energy, cross_secs, opt_out_of_bounds_lower, &
       opt_out_of_bounds_upper, handle_effective)
    use m_units_constants
    use iso_fortran_env, only: error_unit
    !> Name of the file with cross sections
    character(len=*), intent(in)           :: filename
    !> Name of the gas (e.g. N2, O2)
    character(len=*), intent(in)           :: gas_name
    !> Number density of the gas (1/m3)
    real(dp), intent(in)                   :: number_dens
    !> Up to which energy the cross sections are required (eV)
    real(dp), intent(in)                   :: req_energy
    !> The found cross sections are added to this array
    type(CS_t), intent(inout), allocatable :: cross_secs(:)
    !> How to handle data below the first data point, if the first cross section
    !> is non-zero (default: error)
    integer, intent(in), optional          :: opt_out_of_bounds_lower
    !> How to handle data after the last data point, if the last cross section
    !> is non-zero (default: error)
    integer, intent(in), optional          :: opt_out_of_bounds_upper
    !> How to convert EFFECTIVE cross sections (default: do not)
    integer, intent(in), optional          :: handle_effective

    type(CS_t), allocatable :: cs_cpy(:)
    type(CS_t)              :: cs_buf(max_processes_per_gas)
    integer                 :: n, cIx, nL, n_rows, col_type
    integer                 :: my_unit, io_state, len_gas_name
    character(LEN=name_len) :: lineFMT, unit_string
    character(LEN=line_len) :: line, prev_line
    ! Use max_num_rows+1 in case we need to add a data point
    real(dp)                :: cs(2, max_num_rows+1)
    real(dp)                :: x_scaling, y_scaling, tmp_value
    real(dp)                :: two_reals(2), temp
    integer                 :: out_of_bounds_lower, out_of_bounds_upper
    integer                 :: ahandle_effective

    ! Default behaviour for when the input data does not go up to the upper
    ! bound energy (req_energy)
    if (.not. present(opt_out_of_bounds_upper)) then
       out_of_bounds_upper = CS_extrapolate_error
    else
       out_of_bounds_upper = opt_out_of_bounds_upper
    end if

    ! Default behavior for data at lower energies than the input data
    if (.not. present(opt_out_of_bounds_lower)) then
       out_of_bounds_lower = CS_extrapolate_constant
    else
       out_of_bounds_lower = opt_out_of_bounds_lower
    end if

    ahandle_effective = CS_effective_auto
    if (present(handle_effective)) ahandle_effective = handle_effective

    my_unit      = 333
    nL           = 0 ! Set the number of lines to 0
    cIx          = 0
    len_gas_name = len(trim(gas_name))

    ! Set the line format to read, only depends on line_len currently
    write(lineFMT, FMT = "(A,I0,A)") "(A", line_len, ")"

    ! Open 'filename' (with error checking)
    open(my_unit, FILE = trim(filename), ACTION = "READ", &
         ERR = 999, IOSTAT = io_state)

    ! Look for collision processes with the correct gas name in the file,
    ! which should look for example like:

    !     ATTACHMENT                    [description of the type of process, always in CAPS]
    !     H2O -> H2O^-                  [the gas name possibly followed by the result of the process]
    !     COMMENT: total attachment     [possibly comments]
    !     UPDATED: 2010-06-24 15:04:36
    !     SCALING: 1.0 1.0              [optionally scale factors for the columns]
    !     TIMES_N: CM3 (or M3)          [optional; multiply with gas number dens, for 3-body processes]
    !     ------------------            [at least 5 dashes]
    !     xxx   xxx                     [cross section data in two column format]
    !     ...   ...
    !     xxx   xxx
    !     ------------------

    ! So if we find the gas name the previous line holds the type of collision, then
    ! there is possibly some extra information and between the dashes the actual cross
    ! sections are found.

    ! The outer DO loop, running until the end of the file is reached
    do
       ! Search for 'gas_name' in the file
       line = ' '
       do
          prev_line = line
          read(my_unit, FMT = lineFMT, ERR = 999, end = 666) line; nL = nL+1
          line = adjustl(line)
          if (line(1:len_gas_name) == gas_name) exit
       end do

       ! Check prev_line for the type of collision
       select case (prev_line)
       case ("ATTACHMENT")
          col_type = CS_attach_t
       case ("ELASTIC")
          col_type = CS_elastic_t
       case ("EFFECTIVE", "MOMENTUM")
          if (prev_line == "MOMENTUM") print *, "Warning: MOMENTUM deprecated"
          col_type = CS_effective_t
       case ("EXCITATION")
          col_type = CS_excite_t
       case ("IONIZATION")
          col_type = CS_ionize_t
       case ("COMMENT")
          cycle
       case DEFAULT
          print *, "CS_read_file warning: ignoring unknown process type for ", &
               trim(gas_name), " in ", trim(filename), " at line ", nL
          cycle
       end select

       ! Update the number of processes and set the gas name and collision type
       cIx = cIx + 1

       ! Add the reaction description to the table
       cs_buf(cIx)%description = adjustl(line)

       ! For all collisions except attachment, there is a value on the next line
       if (col_type /= CS_attach_t) then
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          read(line, FMT = *, ERR = 999, end = 555) tmp_value
       else
          tmp_value = 0
       end if

       cs_buf(cIx)%coll%type = col_type
       cs_buf(cIx)%coll%part_mass = UC_elec_mass
       if (gas_name == "CH4") cs_buf(cIx)%coll%gas_index = CS_is_CH4

       select case(col_type)
       case (CS_elastic_t, CS_effective_t)
          cs_buf(cIx)%coll%rel_mass = tmp_value ! Store relative mass e/M
       case (CS_excite_t, CS_ionize_t)
          ! Energy loss in Joule
          cs_buf(cIx)%coll%en_loss  = tmp_value * UC_elec_volt
       end select

       cs_buf(cIx)%gas_name = gas_name
       cs_buf(cIx)%comment = "COMMENT: (empty)"
       x_scaling = 1.0_dp
       y_scaling = 1.0_dp

       ! Now we can check whether there is a ZDPLASKIN statement and a
       ! reaction description, while scanning lines until dashes are
       ! found, which indicate the start of the cross section data
       do
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:9) == "ZDPLASKIN" ) then
             cs_buf(cIx)%description = trim(gas_name) // " [" // &
                  trim(adjustl(line(11:))) // "]"
          else if ( line(1:7) == "COMMENT") then
             cs_buf(cIx)%comment = line
          else if ( line(1:7) == "TIMES_N") then
             read(line(9:), *) unit_string
             unit_string = adjustl(unit_string)
             if (unit_string == "CM3") then
                y_scaling = y_scaling * number_dens * 1.0e-6_dp
             else if (unit_string == "M3") then
                y_scaling = y_scaling * number_dens
             else
                write(error_unit, *) "TIMES_N: statement followed by wrong symbol"
                write(error_unit, *) "Allowed are CM3 and M3"
                go to 999
             end if
          else if ( line(1:7) == "SCALING") then
             read(line(9:), *) two_reals
             x_scaling = x_scaling * two_reals(1)
             y_scaling = y_scaling * two_reals(2)
          else if ( line(1:5) == "-----" ) then
             exit
          end if
       end do

       ! Read the cross section data into a temporary array
       n_rows = 0
       do
          read(my_unit, FMT = lineFMT, ERR = 999, end = 555) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit  ! Dashes mark the end of the data
          else if (trim(line) == "" .or. line(1:1) == "#") then
             cycle ! Ignore whitespace or comments
          else if (n_rows < max_num_rows) then
             n_rows = n_rows + 1
             read(line, FMT = *, ERR = 999, end = 555) cs(:, n_rows)
          else
             write(error_unit, *) "Too many rows in ", &
                  trim(filename), " at line ", nL
             error stop
          end if
       end do

       if (n_rows < 2) then
          write(error_unit, *) "Need >= 2 rows in ", &
               trim(filename), " at line number ", nL
          error stop
       end if

       ! Apply scaling
       cs(1, 1:n_rows) = cs(1, 1:n_rows) * x_scaling
       cs(2, 1:n_rows) = cs(2, 1:n_rows) * number_dens * y_scaling

       ! Check whether the tables that have been read in go up to high enough
       ! energies for our simulation. They can also have 0.0 as their highest
       ! listed cross section, in which case we assume the cross section is 0.0
       ! for all higher energies. Cross sections will be taken as constant for
       ! energies above the ones found inside the tables.
       if (cs(1, n_rows) < req_energy .and. cs(2, n_rows) > 0.0D0) then
          select case (out_of_bounds_upper)
          case (CS_extrapolate_error)
             write(error_unit, *) "Cross section data at line ", &
                  nL, " does not go up to high enough x-values (energy)."
             write(error_unit, *) "Required: ", req_energy, "found: ", cs(1, n_rows)
             error stop
          case (CS_extrapolate_zero)
             ! Add an extra line to the end of the input data with a cross
             ! section of 0 and an energy slightly higher than the last value
             n_rows = n_rows + 1
             cs(1, n_rows) = cs(1, n_rows - 1) * (1 + epsilon(1.0_dp))
             cs(2, n_rows) = 0
          case (CS_extrapolate_constant)
             continue
             ! Do nothing
          case (CS_extrapolate_linear)
             ! Linearly extrapolate from the last value in the table to
             ! req_energy. If this cross section value becomes negative we set
             ! it to 0 and warn the user
             cs(1, n_rows + 1) = req_energy
             temp = (req_energy - cs(1, n_rows - 1)) / &
                  (cs(1, n_rows) - cs(1, n_rows - 1))
             cs(2, n_rows + 1) = (1 - temp) * cs(2, n_rows - 1) + &
                  temp * cs(2, n_rows)
             n_rows = n_rows + 1

             if (cs(2, n_rows) < 0) then
                cs(2, n_rows) = 0
                write(error_unit, *) "CS_read_file warning: Linearly ", &
                     " extrapolating cross section to ", req_energy, &
                     " eV resulted in a negative cross section. ", &
                     "Setting this cross section to 0."
             end if
          case default
             write(error_unit, *) "Method with enum ", out_of_bounds_upper, &
                  " for handling input data which does not go to high enough "&
                  "x-values (energy) was not implemented."
             error stop
          end select
       end if

       ! Adjust the tables at energies below the lowest value in the input data
       ! to ensure wanted behavior. Cross sections will be taken as constant for
       ! energies below the ones found inside the tables.
       if (cs(2, 1) > 0.0D0 .and. cs(1, 1) > 0.0D0) then
          select case (out_of_bounds_lower)
          case (CS_extrapolate_constant)
             continue
             ! Do nothing
          case (CS_extrapolate_zero)
             ! Add an extra line at the beginning of the input data with a cross
            ! section of 0 and an energy slightly lower than the first value
            n_rows = n_rows + 1
            ! Shift all array elements (energy and cs) by 1 further in the array
            cs = cshift(cs, shift=-1, dim=2)

            cs(1, 1) = cs(1, 2) * (1 - epsilon(1.0_dp))
            cs(2, 1) = 0

            ! If the collision has a threshold energy we want the extra 0 cs point to be at this threshold
            if (col_type == CS_ionize_t .or. col_type == CS_excite_t) then
               ! tmp_value stores the threshold energy in eV for ionization and excitation collision types
               if (cs(1, 1) > tmp_value) then
                  cs(1, 1) = tmp_value
               end if
            end if
          case default
             error stop "Invalid value for opt_out_of_bounds_upper"
          end select
       end if

       ! Store the data in the actual table
       allocate(cs_buf(cIx)%en_cs(2, n_rows))
       cs_buf(cIx)%n_rows     = n_rows
       cs_buf(cIx)%en_cs(:, :) = cs(:, 1:n_rows)

       ! Locate minimum energy (first value followed by non-zero cross sec)
       cs_buf(cIx)%min_energy = 0.0_dp
       do n = 1, n_rows-1
          if (cs_buf(cIx)%en_cs(2, n+1) > 0.0_dp) then
             cs_buf(cIx)%min_energy = cs_buf(cIx)%en_cs(1, n)
             exit
          end if
       end do

       ! Locate maximum energy (last value preceded by non-zero)
       cs_buf(cIx)%max_energy = 0.0_dp
       do n = n_rows, 2, -1
          if (cs_buf(cIx)%en_cs(2, n-1) > 0.0_dp) then
             cs_buf(cIx)%max_energy = cs_buf(cIx)%en_cs(1, n)
             exit
          end if
       end do

    end do

555 continue ! Routine ends here if the end of "filename" is reached erroneously
    close(my_unit, ERR = 999, IOSTAT = io_state)
    write(error_unit, *) "iostat = ", io_state, " while reading from [", &
         trim(filename), "] at line ", nL
    error stop "End of file while searching"
    return

666 continue ! Routine ends here if the end of "filename" is reached correctly
    close(my_unit, ERR = 999, IOSTAT = io_state)

    if (cIx == 0) then
       write(error_unit, *) "While searching [", &
            trim(gas_name), "] in [", trim(filename), "]"
       error stop "No cross sections found"
    end if

    ! Handle EFFECTIVE cross sections
    n = count(cs_buf(1:cIx)%coll%type == CS_effective_t)

    if (n > 1) then
       write(error_unit, *) "For [", &
            trim(gas_name), "] in [", trim(filename), "]"
       error stop "Multiple EFFECTIVE cross sections"
    else if (n == 1) then
       select case (ahandle_effective)
       case (CS_effective_error)
          print *, "Found EFFECTIVE cross section for ", &
               trim(gas_name), " in ", trim(filename), " at line ", nL
          error stop "handle_effective=CS_effective_error"
       case (CS_effective_auto)
          if (any(cs_buf(1:cIx)%coll%type == CS_elastic_t)) then
             ! Also have ELASTIC cross section, remove EFFECTIVE
             cs_buf(1:cIx-1) = pack(cs_buf(1:cIx), &
                  mask=(cs_buf(1:cIx)%coll%type /= CS_effective_t))
          else
             write(error_unit, *) "Only found EFFECTIVE cross section for [", &
                  trim(gas_name), "] in [", trim(filename), "]"
             error stop
          end if
       case (CS_effective_keep)
          continue
       case default
          error stop "Invalid option for handle_effective"
       end select
    end if

    ! Set the output data
    if (allocated(cross_secs)) then
       n = size(cross_secs)
       call move_alloc(cross_secs, cs_cpy)
       allocate(cross_secs(n + cIx))
       cross_secs(1:n) = cs_cpy
       cross_secs(n+1:) = cs_buf(1:cIx)
    else
       allocate(cross_secs(cIx))
       cross_secs(:) = cs_buf(1:cIx)
    end if

    return

999 continue ! If there was an error, the routine will end here
    write(error_unit, *) "CS_read_file error at line ", nL, &
         " io_state = ", io_state, " while searching [", trim(gas_name), &
         "] in [", trim(filename), "]"
    error stop

  end subroutine CS_add_from_file

  subroutine CS_write_summary(cross_secs, filename)
    use iso_fortran_env, only: error_unit
    type(CS_t), intent(in)       :: cross_secs(:)
    character(LEN=*), intent(in) :: filename
    character(LEN=name_len)      :: col_name
    integer                      :: n, io_state, my_unit
    my_unit = 333

    open(my_unit, FILE = trim(filename), ACTION = "WRITE", &
         ERR = 999, IOSTAT = io_state)

    write(my_unit, ERR = 999, FMT = "(A)") "# List of collision processes"
    write(my_unit, ERR = 999, FMT = "(A)") &
         "Index      Gasname            Coltype     Description           Activation Energy (J)"
    write(my_unit, ERR = 999, FMT = "(A)") &
         "---------------------------------------------------------------------------------------"

    do n = 1, size(cross_secs)
       select case (cross_secs(n)%coll%type)
       case (CS_elastic_t)
          col_name = "Elastic"
       case (CS_effective_t)
          col_name = "Effective"
       case (CS_excite_t)
          col_name = "Excitation"
       case (CS_attach_t)
          col_name = "Attachment"
       case (CS_ionize_t)
          col_name = "Ionization"
       case default
          error stop "Unknown collision type"
       end select

       write(my_unit, ERR = 999, FMT = "((I4),(A),(A12),(A),(A15),(A),(A30),(A),(E13.6))") &
            n, "    ", trim(cross_secs(n)%gas_name), "  ", trim(col_name), &
            "     ", cross_secs(n)%description, " ", cross_secs(n)%coll%en_loss
    end do

    close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)
    return

999 continue ! If there was an error, the routine will end here
    write(error_unit, *) "CS_write_summary error, io_state = ", io_state, &
         " while writing to ", trim(filename)
    error stop

  end subroutine CS_write_summary

  subroutine CS_create_ledger(cross_secs, filename)
    use iso_fortran_env, only: error_unit
    type(CS_t), intent(in)       :: cross_secs(:)
    character(len=*), intent(in) :: filename
    character(len=name_len)      :: col_name
    integer                      :: my_unit, ii, io_state

    my_unit = 333

    open(my_unit, FILE = trim(filename), ACTION = "WRITE", &
         ERR = 999, IOSTAT = io_state)

    write(my_unit, '(A)') "This document stores that the total (=space and time &
      &integrated) number of collisions that have occurred. The output is given  &
      &per collision index (cIx, c.f. *_cs_summary.txt.). The plotting-script &
      &associated with this data can be found in: &
      &<path-to-afivo-pic>/tools/plot_ledger.py"
    write(my_unit, '(A)') " "
    write(my_unit, '(A)') "The next three rows are: 1. cIx , 2. gas molecule &
      &involved in collision, 3. collision type. This legend is followed by the &
      &timesteps the program has completed. Note that a timestamp is appended &
      &to each output line (i.e. the last column of each row corresponds to the &
      &simulated time of that output)"
    write(my_unit, '(A)') " ==================================================== "

    ! Write cIx (given here by loop-counter) to first line
    do ii = 1, (size(cross_secs(:))-1)
      write(my_unit, '(I3, A)', advance = "no") ii, " "
    end do
    write(my_unit, '(I3)') size(cross_secs(:))

    ! Write neutral molecule type to second line
    do ii = 1, (size(cross_secs(:))-1)
      write(my_unit, '(A, A)', advance = "no") trim(cross_secs(ii)%gas_name), " "
    end do
    write(my_unit, '(A)') trim(cross_secs(size(cross_secs(:)))%gas_name)

    ! Write collision type to third line
    do ii = 1, size(cross_secs(:))
      select case (cross_secs(ii)%coll%type)
      case (CS_elastic_t)
        col_name = "Elastic"
      case (CS_excite_t)
        col_name = "Excitation"
      case (CS_attach_t)
        col_name = "Attachment"
      case (CS_ionize_t)
        col_name = "Ionization"
      end select

      if (ii < size(cross_secs(:))) then
        write(my_unit, '(A, A)', advance = "no") trim(col_name), " "
      else
        write(my_unit, '(A)') trim(col_name)
      end if
    end do

    close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)
    return

999 continue ! If there was an error, the routine will end here
    write(error_unit, *) "CS_create_ledger error, io_state = ", io_state, &
         " while writing to ", trim(filename)
    error stop

  end subroutine CS_create_ledger

  subroutine CS_write_ledger(coll_ledger, filename, timestamp)
    use iso_fortran_env, only: error_unit
    real(dp), intent(in)         :: coll_ledger(:)
    character(LEN=*), intent(in) :: filename
    real(dp), intent(in)         :: timestamp
    integer                      :: io_state, my_unit, ii

    my_unit = 333

    open(my_unit, FILE = trim(filename), ACTION = "WRITE", position = "append", &
          ERR = 998, IOSTAT = io_state)

    do ii = 1, size(coll_ledger)
      write(my_unit, '(E13.6)', advance = 'no') coll_ledger(ii)
    end do
    write(my_unit, *) timestamp
    close(my_unit, STATUS = "KEEP", ERR = 999, IOSTAT = io_state)

    return

998 continue ! If there was an error, the routine will end here
    write(error_unit, *) " Error while opening ", trim(filename), ", io_state = ", io_state
    error stop
999 continue ! If there was an error, the routine will end here
    write(error_unit, *) "Error while writing to ", trim(filename), ", io_state = ", io_state
    error stop

  end subroutine CS_write_ledger


end module m_cross_sec
