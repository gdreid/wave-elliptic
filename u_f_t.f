      subroutine u_f_t(R,Z,n_f,n_f_t,np1_f,np1_f_t,NZ,NR,hR,hZ,ht,myzero
     &,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer NZ
      integer NR
      real*8 hR
      real*8 hZ
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
      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=1, NR, 1
      qb = -0.1D1 * n_f_t(i, j) + 0.2D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, NR, 1
      qb = -0.1D1 * n_f_t(i, j) + 0.2D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      qb = -0.1D1 * n_f_t(i, j) + 0.2D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=2, NZ-1, 1
      do j=1, 1, 1
      qb = -0.1D1 * n_f_t(i, j) + 0.1636363636363636D1 * n_f_t(i, j + 1)
     # - 0.8181818181818182D0 * n_f_t(i, j + 2) + 0.1818181818181818D0 *
     # n_f_t(i, j + 3) + 0.1636363636363636D1 * np1_f_t(i, j + 1) - 0.81
     #81818181818182D0 * np1_f_t(i, j + 2) + 0.1818181818181818D0 * np1_
     #f_t(i, j + 3)
      res(i,j)=qb
      end do
      end do
      endif
      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      qb = np1_f_t(i, j) - 0.1D1 * ht * ((-0.1D1 * n_f_t(i, j) + np1_f_t
     #(i, j)) / ht + 0.2500000000000000D0 * (0.1D1 / (exp(R(j)) - 0.1D1)
     # / exp(R(j)) - 0.1D1 / exp(R(j)) ** 2) * (np1_f(i, j - 1) - 0.1D1 
     #* np1_f(i, j + 1)) / hR - 0.5000000000000000D0 / exp(R(j)) ** 2 * 
     #(np1_f(i, j - 1) - 0.2D1 * np1_f(i, j) + np1_f(i, j + 1)) / hR ** 
     #2 - 0.2500000000000000D0 / (exp(Z(i)) + 0.1D1 / exp(Z(i))) ** 3 * 
     #(np1_f(i - 1, j) - 0.1D1 * np1_f(i + 1, j)) / hZ * exp(Z(i)) + 0.2
     #500000000000000D0 / (exp(Z(i)) + 0.1D1 / exp(Z(i))) ** 3 * (np1_f(
     #i - 1, j) - 0.1D1 * np1_f(i + 1, j)) / hZ / exp(Z(i)) - 0.50000000
     #00000000D0 / (exp(Z(i)) + 0.1D1 / exp(Z(i))) ** 2 * (np1_f(i - 1, 
     #j) - 0.2D1 * np1_f(i, j) + np1_f(i + 1, j)) / hZ ** 2 + 0.25000000
     #00000000D0 * (0.1D1 / (exp(R(j)) - 0.1D1) / exp(R(j)) - 0.1D1 / ex
     #p(R(j)) ** 2) * (n_f(i, j - 1) - 0.1D1 * n_f(i, j + 1)) / hR - 0.5
     #000000000000000D0 / exp(R(j)) ** 2 * (n_f(i, j - 1) - 0.2D1 * n_f(
     #i, j) + n_f(i, j + 1)) / hR ** 2 - 0.2500000000000000D0 / (exp(Z(i
     #)) + 0.1D1 / exp(Z(i))) ** 3 * (n_f(i - 1, j) - 0.1D1 * n_f(i + 1,
     # j)) / hZ * exp(Z(i)) + 0.2500000000000000D0 / (exp(Z(i)) + 0.1D1 
     #/ exp(Z(i))) ** 3 * (n_f(i - 1, j) - 0.1D1 * n_f(i + 1, j)) / hZ /
     # exp(Z(i)) - 0.5000000000000000D0 / (exp(Z(i)) + 0.1D1 / exp(Z(i))
     #) ** 2 * (n_f(i - 1, j) - 0.2D1 * n_f(i, j) + n_f(i + 1, j)) / hZ 
     #** 2)
      res(i,j)=qb
      end do
      end do
      END
