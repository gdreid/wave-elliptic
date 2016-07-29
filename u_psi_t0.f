      subroutine u_psi_t0(R,Z,n_f,n_psi_t0,NZ,NR,hR,hZ,myzero,phys_bdy,r
     &es)
      implicit none
      integer i
      integer j
      integer NZ
      integer NR
      real*8 hR
      real*8 hZ
      real*8 myzero
      real*8 R(NR)
      real*8 Z(NZ)
      real*8 n_f(NZ,NR)
      real*8 n_psi_t0(NZ,NR)
      integer phys_bdy(4)
      real*8 res(NZ,NR)
      real*8 qb
      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=1, NR, 1
      qb = myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, NR, 1
      qb = myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      qb = myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=2, NZ-1, 1
      do j=1, 1, 1
      qb = myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      qb = n_psi_t0(i, j) + 0.5000000000000000D0 / (exp(Z(i)) ** 2 * exp
     #(R(j)) ** 2 * hR ** 2 + exp(Z(i)) ** 4 * hZ ** 2 + 0.2D1 * exp(Z(i
     #)) ** 2 * hZ ** 2 + hZ ** 2) * exp(R(j)) ** 2 * hR ** 2 * (exp(Z(i
     #)) ** 2 + 0.1D1) ** 2 * hZ ** 2 * (-0.5000000000000000D0 / (exp(R(
     #j)) - 0.1D1) / exp(R(j)) * (n_psi_t0(i, j - 1) - 0.1D1 * n_psi_t0(
     #i, j + 1)) / hR + 0.5000000000000000D0 / exp(R(j)) ** 2 * (n_psi_t
     #0(i, j - 1) - 0.1D1 * n_psi_t0(i, j + 1)) / hR + 0.1D1 / exp(R(j))
     # ** 2 * (n_psi_t0(i, j - 1) - 0.2D1 * n_psi_t0(i, j) + n_psi_t0(i,
     # j + 1)) / hR ** 2 - 0.5000000000000000D0 / (exp(Z(i)) + 0.1D1 / e
     #xp(Z(i))) ** 3 * (-0.1D1 * n_psi_t0(i - 1, j) + n_psi_t0(i + 1, j)
     #) / hZ * exp(Z(i)) + 0.5000000000000000D0 / (exp(Z(i)) + 0.1D1 / e
     #xp(Z(i))) ** 3 * (-0.1D1 * n_psi_t0(i - 1, j) + n_psi_t0(i + 1, j)
     #) / hZ / exp(Z(i)) + 0.1D1 / (exp(Z(i)) + 0.1D1 / exp(Z(i))) ** 2 
     #* (n_psi_t0(i - 1, j) - 0.2D1 * n_psi_t0(i, j) + n_psi_t0(i + 1, j
     #)) / hZ ** 2 - 0.1D1 * n_f(i, j) ** 2)
      res(i,j)=qb
      end do
      end do
      END
