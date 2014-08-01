module m_particle_par
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: PP_share
  public :: PP_share_mpi
  public :: PP_reorder_by_bins
  public :: PP_reorder_by_bins_mpi

contains

  ! Share particles between PC_t objects
  subroutine PP_share(pcs)
    type(PC_t), intent(inout) :: pcs(:)

    integer :: i_temp(1), n_pc
    integer :: n_avg, i_min, i_max, n_min, n_max, n_send

    n_pc = size(pcs)
    n_avg = ceiling(sum(pcs(:)%n_part) / real(n_pc, dp))

    do
       i_temp = maxloc(pcs(:)%n_part)
       i_max  = i_temp(1)
       n_max  = pcs(i_max)%n_part
       i_temp = minloc(pcs(:)%n_part)
       i_min  = i_temp(1)
       n_min  = pcs(i_min)%n_part

       ! Difference it at most n_pc - 1, if all lists get one more particle
       ! than the last list
       if (n_max - n_min < n_pc) exit

       ! Send particles from i_max to i_min
       n_send = min(n_max - n_avg, n_avg - n_min)
       pcs(i_min)%particles(n_min+1:n_min+n_send) = &
            pcs(i_max)%particles(n_max-n_send+1:n_max)

       ! Always at the end of a list, so do not need to clean up later
       pcs(i_min)%n_part = pcs(i_min)%n_part + n_send
       pcs(i_max)%n_part = pcs(i_max)%n_part - n_send
    end do
  end subroutine PP_share

  subroutine PP_share_mpi(pc, myrank, ntasks)
      use mpi
      type(PC_t), intent(inout) :: pc
      integer, intent(IN) :: myrank, ntasks

      integer :: i_more, i_less, n_send
      integer :: ierr, mean_n_part, i_min, i_max, tag
      integer, dimension(0:ntasks-1) :: n_part_task

      print *, myrank, "Before sharing: pc%n_part", pc%n_part
      ! Get the number of particles each task has
      call MPI_ALLGATHER(pc%n_part, 1, MPI_integer, n_part_task, 1, &
           MPI_integer, MPI_COMM_WORLD, ierr)
      mean_n_part = (sum(n_part_task) + ntasks - 1) / ntasks

      i_less = 0
      i_more = 0

      do while (i_more < ntasks)
         if (n_part_task(i_more) > mean_n_part) then

            ! Find first task with less than mean particles
            do while (n_part_task(i_less) >= mean_n_part)
               i_less = i_less + 1
            end do

            ! Determine amount of particles to send
            n_send = min(n_part_task(i_more) - mean_n_part, &
                 mean_n_part - n_part_task(i_less))

            if (myrank == i_more) then ! Send
               i_min = n_part_task(i_more) - n_send + 1
               i_max = i_min + n_send - 1
               tag   = i_more
               call send_parts_mpi(pc%particles(i_min:i_max), i_less, tag)
            else if (myrank == i_less) then ! Receive
               i_min = n_part_task(i_less) + 1
               i_max = i_min + n_send - 1
               call pc%check_space(i_max)
               tag   = i_more
               call recv_parts_mpi(pc%particles(i_min:i_max), i_more, tag)
            end if

            ! Update the number of particles for each task
            n_part_task(i_less) = n_part_task(i_less) + n_send
            n_part_task(i_more) = n_part_task(i_more) - n_send
         end if

         if (n_part_task(i_more) <= mean_n_part) i_more = i_more + 1
      end do

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      pc%n_part = n_part_task(myrank)
      print *, myrank, "After sharing: pc%n_part", pc%n_part
  end subroutine PP_share_mpi

  subroutine PP_reorder_by_bins(pcs, bin_func, n_bins, bin_func_args)
    use m_mrgrnk
    use m_linked_list
    type(PC_t), intent(inout) :: pcs(:)
    interface
       integer function bin_func(my_part, r_args)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(in) :: r_args(:)
       end function bin_func
    end interface
    integer, intent(in) :: n_bins
    real(dp), intent(in) :: bin_func_args(:)
    integer, allocatable :: bin_counts(:,:)
    integer, allocatable :: bin_counts_sum(:)
    integer, allocatable :: bin_owner(:)
    integer :: prev_num_part(size(pcs))
    integer :: n_pc, n_avg, tmp_sum
    integer :: ll, io, ib, ip, new_loc

    type tmp_t
       integer, allocatable :: ib(:)
    end type tmp_t
    type(tmp_t), allocatable :: pcs_bins(:)

    n_pc = size(pcs)
    n_avg = ceiling(sum(pcs(:)%n_part) / real(n_pc, dp))

    allocate(bin_counts(n_bins, n_pc))
    allocate(bin_counts_sum(n_bins))
    allocate(bin_owner(n_bins))
    allocate(pcs_bins(n_pc))
    bin_counts(:, :) = 0

    ! Get the counts in the bins for each pcs(ip)
    do ip = 1, n_pc
       allocate(pcs_bins(ip)%ib(pcs(ip)%n_part))
       do ll = 1, pcs(ip)%n_part
          ib                  = bin_func(pcs(ip)%particles(ll), bin_func_args)
          pcs_bins(ip)%ib(ll) = ib
          bin_counts(ib, ip)  = bin_counts(ib, ip) + 1
       end do
    end do

    bin_counts_sum = sum(bin_counts, dim=2)
    tmp_sum        = 0
    ip             = 1

    ! Set the owners of the bins
    do ib = 1, n_bins
       tmp_sum = tmp_sum + bin_counts_sum(ib)
       bin_owner(ib) = ip
       if (tmp_sum >= n_avg) then
          ip = ip + 1
          tmp_sum = 0
       end if
    end do

    prev_num_part(:) = pcs(:)%n_part

    do ip = 1, n_pc
       do ll = 1, prev_num_part(ip)
          ib = pcs_bins(ip)%ib(ll)
          io = bin_owner(ib)
          if (io /= ip) then
             ! Insert at owner
             new_loc = pcs(io)%n_part + 1
             pcs(io)%particles(new_loc) = pcs(ip)%particles(ll)
             pcs(io)%n_part = new_loc
             call LL_add(pcs(ip)%clean_list, ll)
          end if
       end do
    end do

    do ip = 1, n_pc
       call pcs(ip)%clean_up()
       print *, "REMOVECHECK", ip, pcs(ip)%n_part
    end do
  end subroutine PP_reorder_by_bins

  ! Divide the particles over the tasks based on their bin index
  subroutine PP_reorder_by_bins_mpi(pc, bin_func, n_bins, bin_func_args, &
       myrank, ntasks)
    use mpi
    use m_lookup_table
    use m_mrgrnk
    use m_linked_list
    type(PC_t), intent(inout) :: pc
    real(dp), intent(in) :: bin_func_args(:)
    integer, intent(in) :: n_bins, myrank, ntasks

    interface
       integer function bin_func(my_part, r_args)
         import
         type(PC_part_t), intent(in) :: my_part
         real(dp), intent(in) :: r_args(:)
       end function bin_func
    end interface

    integer :: ll, ix, rank, ierr, n_behind, n_part_interval, n_send
    integer :: i_min, i_max
    integer :: tag, part_count
    integer :: n_sends, n_recvs, sender, recver, n_part_mean
    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: n_part_task(0:ntasks-1)
    integer :: n_part_interval_task(0:ntasks-1, 0:ntasks-1)
    integer :: split_ix(-1:ntasks)
    integer, allocatable :: n_part_per_bin(:), ix_list(:)
    real(dp), allocatable :: bin_ixs(:)

    print *, "Before divide ", myrank, " has ", pc%n_part, " particles"

    ! Get the number of particles each task has
    call MPI_ALLGATHER(pc%n_part, 1, MPI_integer, n_part_task, 1, &
         MPI_integer, MPI_COMM_WORLD, ierr)

    n_part_mean = (sum(n_part_task) + ntasks - 1) / ntasks
    allocate(n_part_per_bin(n_bins))
    allocate(bin_ixs(pc%n_part))
    allocate(ix_list(pc%n_part))
    n_part_per_bin = 0

    ! Set the bin indexes for the particles
    do ll = 1, pc%n_part
       ix                   = bin_func(pc%particles(ll), bin_func_args)
       bin_ixs(ll)          = ix
       n_part_per_bin(ix)   = n_part_per_bin(ix) + 1
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, n_part_per_bin, n_bins, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Sort by bin_ixs
    call mrgrnk(bin_ixs, ix_list)
    pc%particles(1:pc%n_part) = pc%particles(ix_list)

    ! Find the split_ix(:) s.t. task n will hold particles with ix >
    ! split_ix(n-1) and ix <= split_ix(n)
    split_ix(-1)              = 0
    split_ix(ntasks)          = n_bins
    rank                      = 0
    n_behind                  = 0
    n_part_interval_task(:,:) = 0
    n_part_interval           = 0

    do ix = 1, n_bins
       n_part_interval = n_part_interval + n_part_per_bin(ix)

       if (n_part_interval >= n_part_mean .or. ix == n_bins) then
          split_ix(rank)  = ix
          n_part_interval = 0
          part_count      = 0

          do ll = n_behind+1, pc%n_part
             if (bin_ixs(ll) > ix) exit
             part_count = part_count + 1
          end do

          n_part_interval_task(rank, myrank) = part_count
          n_behind                           = n_behind + part_count
          rank                               = rank + 1

          if (rank == ntasks) exit
       end if

    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, n_part_interval_task, ntasks*ntasks, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Send particles to task that holds a region
    n_sends = 0
    n_recvs = 0

    do recver = 0, ntasks - 1
       if (myrank == recver) then ! Receive
          do sender = 0, ntasks - 1
             if (sender == myrank) cycle ! Skip ourselves
             n_send = n_part_interval_task(recver, sender)

             if (n_send > 0) then
                n_recvs   = n_recvs + 1
                tag       = sender
                i_min     = pc%n_part + 1
                i_max     = pc%n_part + n_send
                pc%n_part = i_max
                call pc%check_space(i_max)
                print *, myrank, ": receives", n_send, "from", sender, "i_min", i_min
                call recv_parts_mpi(pc%particles(i_min:i_max), sender, tag)
             end if
          end do
       else ! Send
          n_send = n_part_interval_task(recver, myrank)

          if (n_send > 0) then
             n_sends = n_sends + 1
             tag     = myrank

             ! Find index of first particle to send
             i_min    = sum(n_part_interval_task(0:recver-1, myrank)) + 1
             i_max = i_min + n_send - 1

             print *, myrank, ": sends", n_send, "to", recver, "i_min", i_min
             call send_parts_mpi(pc%particles(i_min:i_max), recver, tag)

             ! Mark the particles that we have send away as inactive
             do ll = i_min, i_max
                call LL_add(pc%clean_list, ll)
             end do

          end if
       end if
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, "After divide ", myrank, " has ", pc%n_part, "pc%n_part (some dead)"
    call pc%clean_up()
    print *, "After divide ", myrank, " has ", pc%n_part, " particles"

  end subroutine PP_reorder_by_bins_mpi

  subroutine send_parts_mpi(parts, recver, tag)
    use mpi
    type(PC_part_t), intent(in) :: parts(:)
    integer, intent(in) :: recver, tag

    integer :: n_send, req_size, rel_size, ierr
    real(dp), allocatable :: r_buf(:)

    ! We use a buffer because we cannot send the part type *neatly* with mpi
    rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
    n_send   = size(parts)
    req_size = n_send * rel_size
    allocate(r_buf(req_size))

    r_buf = transfer(parts, r_buf, req_size)
    call MPI_SEND(r_buf, req_size, MPI_DOUBLE_PRECISION, recver, tag, &
         & MPI_COMM_WORLD, ierr)
  end subroutine send_parts_mpi

  subroutine recv_parts_mpi(parts, sender, tag)
    use mpi
    type(PC_part_t), intent(inout) :: parts(:)
    integer, intent(in) :: sender, tag

    integer :: n_send, req_size, rel_size, ierr
    integer :: mpi_status(MPI_STATUS_SIZE)
    real(dp), allocatable :: r_buf(:)

    ! We use a buffer because we cannot send the part type *neatly* with mpi
    n_send   = size(parts)
    rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
    req_size = n_send * rel_size
    allocate(r_buf(req_size))

    call MPI_RECV(r_buf, req_size, MPI_DOUBLE_PRECISION, sender, tag, &
         & MPI_COMM_WORLD, mpi_status, ierr)
    parts = transfer(r_buf, parts, n_send)
  end subroutine recv_parts_mpi

end module m_particle_par