      subroutine eval_res_f(R,Z,n_f,n_f_t,np1_f,np1_f_t,NZ,NR,ht,
     &myzero,phys_bdy,res)
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
      real*8 res(NZ,NR)
      real*8 qb
      res = 0.0D0
      do i=NZ, NZ, 1
      do j=1, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res(i,j)=qb
      end do
      end do
      do i=1, 1, 1
      do j=1, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res(i,j)=qb
      end do
      end do
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.1D1 * myzero * Z(i) * R(
     #j) + 0.5000000000000000D0 * np1_f(i, j)
      res(i,j)=qb
      end do
      end do
      do i=2, NZ-1, 1
      do j=1, 1, 1
      qb = 0.5000000000000000D0 * n_f(i, j) - 0.8181818181818182D0 * n_f
     #(i, j + 1) + 0.4090909090909091D0 * n_f(i, j + 2) - 0.909090909090
     #9091D-1 * n_f(i, j + 3) + 0.5000000000000000D0 * np1_f(i, j) - 0.8
     #181818181818182D0 * np1_f(i, j + 1) + 0.4090909090909091D0 * np1_f
     #(i, j + 2) - 0.9090909090909091D-1 * np1_f(i, j + 3)
      res(i,j)=qb
      end do
      end do
      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      qb = (-0.1D1 * n_f(i, j) + np1_f(i, j)) / ht - 0.5000000000000000D
     #0 * np1_f_t(i, j) - 0.5000000000000000D0 * n_f_t(i, j)
      res(i,j)=qb
      end do
      end do
      END
