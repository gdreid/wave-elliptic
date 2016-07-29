      subroutine init_f(R,Z,NZ,NR,delr,delz,famp,r0,z0,res)
      implicit none
      integer i
      integer j
      integer NZ
      integer NR
      real*8 delr
      real*8 delz
      real*8 famp
      real*8 r0
      real*8 z0
      real*8 R(NR)
      real*8 Z(NZ)
      real*8 res(NZ,NR)
      real*8 qb
      do i=1, NZ, 1
      do j=1, NR, 1
      qb = famp * exp(-0.1D1 * (0.2D1 * exp(-0.1D1 * Z(i)) * z0 * delr *
     #* 2 - 0.2D1 * exp(Z(i)) * z0 * delr ** 2 - 0.2D1 * exp(R(j)) * r0 
     #* delz ** 2 + z0 ** 2 * delr ** 2 + r0 ** 2 * delz ** 2 + exp(-0.2
     #D1 * Z(i)) * delr ** 2 + exp(0.2D1 * Z(i)) * delr ** 2 + exp(0.2D1
     # * R(j)) * delz ** 2 - 0.2D1 * exp(R(j)) * delz ** 2 + 0.2D1 * r0 
     #* delz ** 2 - 0.2D1 * delr ** 2 + delz ** 2) / delr ** 2 / delz **
     # 2)
      res(i,j)=qb
      end do
      end do
      END
