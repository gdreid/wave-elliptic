      subroutine res_psi(n_psi,np1_psi,NZ,NR,ht,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer NZ
      integer NR
      real*8 ht
      real*8 n_psi(NZ,NR)
      real*8 np1_psi(NZ,NR)
      integer phys_bdy(4)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=1, NR, 1
      qb = (-0.1D1 * n_psi(i, j) + np1_psi(i, j)) / ht + 0.5000000000000
     #000D0 * np1_psi(i, j) + 0.5000000000000000D0 * n_psi(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, NR, 1
      qb = (-0.1D1 * n_psi(i, j) + np1_psi(i, j)) / ht + 0.5000000000000
     #000D0 * np1_psi(i, j) + 0.5000000000000000D0 * n_psi(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      qb = (-0.1D1 * n_psi(i, j) + np1_psi(i, j)) / ht + 0.5000000000000
     #000D0 * np1_psi(i, j) + 0.5000000000000000D0 * n_psi(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=2, NZ-1, 1
      do j=1, 1, 1
      qb = (-0.1D1 * n_psi(i, j) + np1_psi(i, j)) / ht + 0.5000000000000
     #000D0 * np1_psi(i, j) + 0.5000000000000000D0 * n_psi(i, j)
      res = res + qb**2
      end do
      end do
      endif
      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      qb = (-0.1D1 * n_psi(i, j) + np1_psi(i, j)) / ht + 0.5000000000000
     #000D0 * np1_psi(i, j) + 0.5000000000000000D0 * n_psi(i, j)
      res = res + qb**2
      end do
      end do
      res = sqrt(res/(1*NZ*NR))
      END
