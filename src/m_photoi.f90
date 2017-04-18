! Jannis Teunissen, Casper Rutjes
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

!> Module that adds photo ionization functionallity to particle core
module m_photoi
  use m_random
  use m_units_constants

  implicit none
  private

  integer, parameter                     :: dp = kind(0.0d0)
  
  !> The settings for the photo ionization
  type PI_t
    logical                                :: enable = .false.
    real(dp)                               :: quench_fac
    real(dp)                               :: min_inv_abs_len
    real(dp)                               :: max_inv_abs_len
    integer                                :: size_photo_eff_table
    real(dp), dimension(:, :), allocatable :: photo_eff_table
    type(RNG_t) :: rng
    integer     :: separator(100) ! Separate rng data
  
  contains
    procedure, non_overridable  :: get_photoi_eff
    procedure, non_overridable  :: get_photoi_lambda
  end type PI_t

  public PI_t
  public PI_initialize
  public PI_from_ionization
  
contains
  
  subroutine PI_initialize(photoi, enabled, quench_fac, min_inv_abs_len, &
       max_inv_abs_len, size_photo_eff_table, photo_eff_table, rng_seed)
    use m_units_constants
    type(PI_t), intent(inout)         :: photoi
    logical, intent(in)               :: enabled
    real(dp), intent(in)              :: quench_fac
    real(dp), intent(in)              :: min_inv_abs_len
    real(dp), intent(in)              :: max_inv_abs_len
    integer, intent(in)               :: size_photo_eff_table
    real(dp), intent(in), allocatable :: photo_eff_table(:,:)
    integer, intent(in), optional :: rng_seed(4)

    integer, parameter            :: i8 = selected_int_kind(18)
    integer(i8)                   :: rng_seed_8byte(2)
    
    photoi%enabled = enabled
    
    photoi%quench_fac = quench_fac
    photoi%min_inv_abs_len = min_inv_abs_len
    photoi%max_inv_abs_len = max_inv_abs_len
    
    photoi%size_photo_eff_table = size_photo_eff_table
    allocate(photoi%photo_eff_table(2,photoi%size_photo_eff_table))
    photoi%photo_eff_table = photo_eff_table
        
    if (present(rng_seed)) then
       rng_seed_8byte = transfer(rng_seed, rng_seed_8byte)
       call photoi%rng%set_seed(rng_seed_8byte)
    else
       call photoi%rng%set_seed([8972134_i8, 21384823409_i8])
    end if    
  end subroutine PI_initialize

  subroutine PI_from_ionization(pc,ll)
    use m_particle_core
    use m_units_constants
    type(PC_t), intent(inout)            :: pc
    integer, intent(in)                  :: ll
    real(dp)                             :: mean_gammas, en_frac, fly_len
    real(dp)                             :: fld, psi, chi, x_end(3)
    integer                              :: n, n_photons

    fld         = norm2(my_part%a / UC_elec_q_over_m)
    mean_gammas = pc%photoi%get_photoi_eff(fld) * my_part%w * pc%photoi%quench_fac
    n_photons   = pc%photoi%rng%poisson(mean_gammas)

    do n = 1, n_photons
       ! Select random direction and absorption length
       en_frac  = pc%photoi%rng%unif_01()
       fly_len  = -log(1.0_dp - pc%photoi%rng%unif_01()) / pc%photoi%get_photoi_lambda(en_frac)
       psi      = 2 * UC_photoi * pc%photoi%rng%unif_01()
       chi      = acos(1.0_dp - 2 * pc%photoi%rng%unif_01())

       x_end(1) = pc%particles(ll)%x(1) + fly_len * sin(chi) * cos(psi)
       x_end(2) = pc%particles(ll)%x(2) + fly_len * sin(chi) * sin(psi)
       x_end(3) = pc%particles(ll)%x(3) + fly_len * cos(chi)
       
       call pc%create_part_and_outside_check(x_end, & ! position (m)
            &                   (/ 0.0D0, 0.0D0, 0.0D0/), & ! velocity (m/s)
            &                   (/0.0D0, 0.0D0, 0.0D0/),  & ! accelaration (m/s2) 
            &                                      1.0D0, & ! weight 
            &                                     0.0D0)   ! time left for move_and_collide (s)
    end do
  end subroutine PI_from_ionization

  ! Returns the photo-efficiency coefficient corresponding to an electric
  ! field of strength fld
  real(dp) function get_photoi_eff(self,fld)
    use m_lookup_table
    class(PI_t), intent(in) :: self
    real(dp), intent(in) :: fld
    call LT_lin_interp_list(self%photo_eff_table(1,:), &
         self%photo_eff_table(2,:), fld, get_photoi_eff)
  end function get_photoi_eff

   ! Returns the inverse mean free path for a photon.
   real(dp) function get_photoi_lambda(self,en_frac)
    class(PI_t), intent(in) :: self
    real(dp), intent(in) :: en_frac
      get_photoi_lambda = self%min_inv_abs_len * &
           (self%max_inv_abs_len/self%min_inv_abs_len)**en_frac
   end function
end module m_photoi