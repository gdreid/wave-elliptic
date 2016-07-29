      subroutine relax_psi_t0(R,Z,n_f,n_psi_t0,n_psi_t0_rhs,NZ,NR,hR,hZ,
     &myzero,phys_bdy,cmask,res)
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
      real*8 n_psi_t0_rhs(NZ,NR)
      integer phys_bdy(4)
      real*8 res(NZ,NR)
      real*8 qb
      
      real*8 cmask(NZ,NR)
      include 'cmask.inc'

      if (phys_bdy(2) .eq. 1) then
      do i=NZ, NZ, 1
      do j=2, NR-1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = 0.1D0*n_psi_t0(i,j) + 0.9D0*(n_psi_t0_rhs(i, j) 
     #+ myzero * Z(i) * R(j))
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=2, NR-1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = 0.1D0*n_psi_t0(i,j) + 0.9D0*(n_psi_t0_rhs(i, j) 
     #+ myzero * Z(i) * R(j))
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
      if (phys_bdy(4) .eq. 1) then
      do i=1, NZ, 1
      do j=NR, NR, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = 0.1D0*n_psi_t0(i,j) + 0.9D0*(n_psi_t0_rhs(i, j) 
     #+ myzero * Z(i) * R(j))
      res(i,j)=qb
      endif
      end do
      end do
      endif
      
c======================================================================
c Neumann boundary conditions at R=0
c======================================================================
      
c      if (phys_bdy(3) .eq. 1) then
c      do i=1, NZ, 1
c      do j=1, 1, 1
c      if (cmask(i,j).eq.CMASK_ON) then
c      qb = 0.1D0*n_psi_t0(i,j) + 0.9D0*(0.1333333333333333D1 
c     #* n_psi_t0(i, j + 1) - 0.333333333333333
c     #3D0 * n_psi_t0(i, j + 2) + n_psi_t0_rhs(i, j))
c      res(i,j)=qb
c      endif
c      end do
c      end do
c      endif
      
c======================================================================
c Dirichlet boundary conditions at R=0
c======================================================================

      if (phys_bdy(3) .eq. 1) then
      do i=1, NZ, 1
      do j=1, 1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = 0.1D0*n_psi_t0(i,j) + 0.9D0*(n_psi_t0_rhs(i, j) 
     #+ myzero * Z(i) * R(j))
      res(i,j)=qb
      endif
      end do
      end do
      endif

      do i=2, NZ-1, 1
      do j=2, NR-1, 1
      if (cmask(i,j).eq.CMASK_ON) then
      qb = n_psi_t0(i, j) + 0.4500000000000000D0 / (hR ** 2 + hZ ** 2) *
     # hR ** 2 * hZ ** 2 * ((n_psi_t0(i, j - 1) - 0.2D1 * n_psi_t0(i, j)
     # + n_psi_t0(i, j + 1)) / hR ** 2 + (n_psi_t0(i - 1, j) - 0.2D1 * n
     #_psi_t0(i, j) + n_psi_t0(i + 1, j)) / hZ ** 2 - 0.1D1 * n_f(i, j) 
     #** 2 - 0.1D1 * n_psi_t0_rhs(i, j))
      res(i,j)=qb
      endif
      end do
      end do
      
c      if (NZ .eq. 9 .and. NR .eq. 5) then
c      write(*,*) 'bdry: ', phys_bdy(1), phys_bdy(2), 
c     #phys_bdy(3), phys_bdy(4)
c      write(*,*) 'relax diff:'
c      do i=1, NZ,1
c      write(*,"(30e15.5)") ( (n_psi_t0(i, j) - 0.1333333333333333D1 
c     #* n_psi_t0(i, j + 1) + 
c     #0.3333333333333333D0 * n_psi_t0(i, j + 2) - 0.1D1 * n_psi_t0_rhs(i
c     #, j)), j=1,NR )
c      end do
c      endif
      
      END
