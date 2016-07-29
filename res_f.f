      subroutine res_f(R,Z,n_f,n_f_t,np1_f,np1_f_t,NZ,NR,ht,myzero,phys_
     &bdy,res)
      implicit none
      integer i
      integer j
      integer NZ
      integer NR
      real*8 ht
      real*8 myzero
      real*8 R(NR)
      real*8 Z(NZ)
      real*8 n_f(NZ,NR)
      real*8 n_f_t(NZ,NR)
      real*8 np1_f(NZ,NR)
      real*8 np1_f_t(NZ,NR)
      integer phys_bdy(4)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=1, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res = res + qb**2
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=2, NZ-1, 1
      do j=1, 1, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.8181818181818182D0 * n_f
     #(i, j + 1) + 0.4090909090909091D0 * n_f(i, j + 2) - 0.909090909090
     #9091D-1 * n_f(i, j + 3) + 0.5000000000000000D0 * np1_f(i, j) - 0.8
     #181818181818182D0 * np1_f(i, j + 1) + 0.4090909090909091D0 * np1_f
     #(i, j + 2) - 0.9090909090909091D-1 * np1_f(i, j + 3)
      res = res + qb**2
      end do
      end do
      endif
      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      qb = (-0.1D1 * n_f(i, j) + np1_f(i, j)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, j) - 0.5000000000000000D0 * n_f_t(i, j)
      res = res + qb**2
      end do
      end do
      res = sqrt(res/(1*NZ*NR))
      END
