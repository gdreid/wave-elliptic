      subroutine lop_psi_t0(R,Z,n_f,n_psi_t0,NZ,NR,hR,hZ,myzero,phys_bdy
     &,cmask,res)
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
      
      real*8 cmask(NZ,NR)
      include 'cmask.inc'

      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=1, NR, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = n_psi_t0(i, j) - 0.1D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, NR, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = n_psi_t0(i, j) - 0.1D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
      if (phys_bdy(4) .eq. 1) then
      do i=2, NZ-1, 1
      do j=NR, NR, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = n_psi_t0(i, j) - 0.1D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
c======================================================================
c Neumann boundary conditions at R=0
c======================================================================
      
c      if (phys_bdy(3) .eq. 1) then
c      do i=2, NZ-1, 1
c      do j=1, 1, 1
c      if (cmask(i,j).eq.CMASK_ON) then
c      qb = n_psi_t0(i, j) - 0.1333333333333333D1 * n_psi_t0(i, j + 1) + 
c     #0.3333333333333333D0 * n_psi_t0(i, j + 2)
c      res(i,j)=qb
c      endif
c      end do
c      end do
c      endif
      
c======================================================================
c Dirichlet boundary conditions at R=0
c======================================================================

      if (phys_bdy(3) .eq. 1) then
      do i=2, NZ-1, 1
      do j=1, 1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = n_psi_t0(i, j) - 0.1D1 * myzero * Z(i) * R(j)
      res(i,j)=qb
      endif
      end do
      end do
      endif

      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = (n_psi_t0(i, j - 1) - 0.2D1 * n_psi_t0(i, j) + n_psi_t0(i, j 
     #+ 1)) / hR ** 2 + (n_psi_t0(i - 1, j) - 0.2D1 * n_psi_t0(i, j) + n
     #_psi_t0(i + 1, j)) / hZ ** 2 - 0.1D1 * n_f(i, j) ** 2
      res(i,j)=qb
      endif
      end do
      end do
      END
