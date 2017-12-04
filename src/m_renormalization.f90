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

!> Program to renormalize cross sections according anisotropic scattering
!> following Chao Li et al. (2012).
program m_renormalization



  use m_cross_sec
  use m_units_constants
  implicit none

  integer, parameter  :: dp        = kind(0.0d0)



  type(cs_t), allocatable :: cross_secs(:)
  real(dp) :: tmp_rn
  real(dp) :: max_energy_eV = 1000.0
  character(len=80), allocatable :: cross_files(:)
  character(len=80), allocatable :: molecule_names(:)
  integer                        :: cIx, j,n               ! coll index
  integer                        :: n_gas_comp        = 2        ! number of gas components
  integer :: file_unit = 87
  character(len=80) :: filename
  real(dp) :: dens = 1.0



    allocate(cross_files(n_gas_comp))
    allocate(molecule_names(n_gas_comp))
    cross_files = ["cs_example.txt","cs_example.txt"]   
    molecule_names = ["N2","O2"]
    filename = "renormalized_x-sections.txt"
    filename = trim(filename)
  do n = 1, n_gas_comp
     call CS_add_from_file(cross_files(n), molecule_names(n), &
          dens, max_energy_eV, cross_secs)
  end do

 do cIx = 1,size(cross_secs)
       if ((trim(cross_secs(cIx)%gas_name) == 'N2').and.(cross_secs(cIx)%coll%type == CS_elastic_t))then
          do j=1,cross_secs(cIx)%n_rows
             if(cross_secs(cIx)%en_cs(1,j) > 0 ) then
                tmp_rn = (0.065*cross_secs(cIx)%en_cs(1,j)+0.26* sqrt(cross_secs(cIx)%en_cs(1,j)))/&
                (1.0 + 0.05*cross_secs(cIx)%en_cs(1,j)+0.2* sqrt(cross_secs(cIx)%en_cs(1,j))) -&
                12.0*sqrt(cross_secs(cIx)%en_cs(1,j))/(1.0+40.0*sqrt(cross_secs(cIx)%en_cs(1,j)))
                ! To avoid division by zero
                if ( abs(tmp_rn) < 64.0*UC_tiny ) then
                    tmp_rn = 1.0 ! This is the limit of the formula below for vec -> 0
                else
                    tmp_rn = (1.0-tmp_rn)/(2.0 * tmp_rn**2)* &
                    ((1.0+tmp_rn)*log((1.0+tmp_rn)/(1.0-tmp_rn))-2.0*tmp_rn)
                end if
                if(abs(tmp_rn) > 0) then
                    cross_secs(cIx)%en_cs(2,j) = cross_secs(cIx)%en_cs(2,j)/tmp_rn
                end if
              end if                     
            end do
         end if
         if ((trim(cross_secs(cIx)%gas_name) == 'O2').and.(cross_secs(cIx)%coll%type == CS_elastic_t))then
             do j=1,cross_secs(cIx)%n_rows
                if(cross_secs(cIx)%en_cs(1,j) > 0 ) then
                   tmp_rn =(2.0/log(1.0+cross_secs(cIx)%en_cs(1,j)))*&
                   (1.0 - log(1.0+cross_secs(cIx)%en_cs(1,j))/cross_secs(cIx)%en_cs(1,j))
                   if(abs(tmp_rn) >0) then 
                     cross_secs(cIx)%en_cs(2,j) = cross_secs(cIx)%en_cs(2,j)/tmp_rn
                   end if
                end if
             end do
         end if
         if ((trim(cross_secs(cIx)%gas_name) == 'Ar').and.(cross_secs(cIx)%coll%type == CS_elastic_t))then
            do j=1,cross_secs(cIx)%n_rows
                if(cross_secs(cIx)%en_cs(1,j) > 0 ) then
                  tmp_rn =(2.0/log(1.0+cross_secs(cIx)%en_cs(1,j)))*&
                  (1.0 - log(1.0+cross_secs(cIx)%en_cs(1,j))/cross_secs(cIx)%en_cs(1,j))
                  if(abs(tmp_rn) >0) then 
                    cross_secs(cIx)%en_cs(2,j) = cross_secs(cIx)%en_cs(2,j)/tmp_rn
                  end if
                end if
             end do
         end if
  end do


 open(unit=file_unit, file=filename,Access = 'append',STATUS="UNKNOWN")

do cIx = 1,size(cross_secs)

write(file_unit,*) cross_secs(cIx)%gas_name, 'collision type:', cross_secs(cIx)%coll%type
write(file_unit,*) cross_secs(cIx)%description
write(file_unit,*) "energy in eV // cross section in m2"
write(file_unit,*) '--------------------------------------------------'
do j=1,cross_secs(cIx)%n_rows
write(file_unit,*) cross_secs(cIx)%en_cs(1,j) , cross_secs(cIx)%en_cs(2,j)
end do
write(file_unit,*) '--------------------------------------------------'

end do

 close(file_unit)



end program m_renormalization
