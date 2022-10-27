
      subroutine stackima(pix,flag,spix,NPX,NPY,NIMS)
      implicit none
      
      integer NPX,NPY,NIMs, NIM
      real pix(NIMs, NPY, NPX)
      real spix(NPY, NPX)

      real mbar_sky, ss(NIMs)
      integer i, ix, jy, j, Ls
      !real pixs(NIMs, NPY, NPX)
      integer flag(NIMs)
      real tomed(NIMs), bar, sig


      !f2py intent(in) :: pix, flag
      !f2py intent(hide), depend(pix) :: NIMs = shape(pix,0)
      !f2py intent(hide), depend(pix) :: NPX  = shape(pix,2)
      !f2py intent(hide), depend(pix) :: NPY  = shape(pix,1)
      !f2py intent(out) spix


      ix = 16
      jy = 16
      do NIM = 1, NIMs
         ss(NIM) = mbar_sky(ix,jy,7,15,pix(NIM,:,:),NPX,NPY)
      enddo

      do i = 1, NPX
         do j = 1, NPY
            Ls = 0
            do NIM = 1, NIMs
               if (flag(NIM).gt.0) cycle
               Ls = Ls + 1
               tomed(Ls) = pix(NIM,i,j)  - ss(NIM)
            enddo
            call barsig(tomed,ls,bar,sig)
            spix(i,j) = bar
         enddo
      enddo


      end subroutine

      subroutine centroid(pix,mag,xc,yc,sxc,syc,NPX,NPY)
      implicit none

      integer NPX, NPY
      real pix(NPX,NPY)
      real pixx(NPX,NPY)
      real mag, lato
      real xc, yc, xorin, yorin
      real sxc, syc
      integer i, j, ii, jj
      integer xdiff, ydiff
      integer LLs
      real, allocatable:: tomedx(:),tomedy(:),tomedz(:)

      !f2py intent(in) :: pix, mag
      !f2py intent(hide), depend(pix) :: NPX  = shape(pix,1)
      !f2py intent(hide), depend(pix) :: NPY  = shape(pix,0)
      !f2py intent(out) xc, yc, sxc, syc

      if (mag>12.75) lato = 1
      if (mag<12.75.and.mag>9.5) lato = 2
      if (mag<9.5.and.mag>6) lato = 3
      if (mag<6) lato = 4

      do i = 1, NPX
         do j = 1, NPY
            pixx(i,j) = pix(j,i)
         enddo
      enddo

      xorin = 16.
      yorin = 16.
      xc = xorin
      yc = yorin
      xdiff = nint(xorin)
      ydiff = nint(yorin)
      do i = -3, 3
         do j = -3, 3
            ii = nint(xorin-i)
            jj = nint(yorin-j)
            if (pixx(ii,jj).ge.pixx(xdiff,ydiff)) then
               xdiff = ii
               ydiff = jj
            endif
         enddo
      enddo
      LLs = 0
      do i = -5, 5
         do j = -5, 5
            ii = xdiff-i
            jj = ydiff-j
            if (sqrt((1.*ii-xdiff)**2+(1.*jj-ydiff)**2)
     .                       .le.1.0*lato) then 
               LLs = LLs + 1
            endif
         enddo
      enddo
      allocate(tomedx(LLs), tomedy(LLs), tomedz(LLs))
      LLs=0
      do i = -5, 5
         do j = -5, 5
            ii = xdiff-i
            jj = ydiff-j
            if (sqrt((1.*ii-xdiff)**2+(1.*jj-ydiff)**2)
     .                       .le.1.0*lato) then 
               LLs = LLs + 1
               tomedx(LLs) = ii
               tomedy(LLs) = jj
               tomedz(LLs) = (pixx(ii,jj))
            endif
         enddo
      enddo
      
      call barsigw(tomedx, tomedz, LLs, xc, sxc)
      call barsigw(tomedy, tomedz, LLs, yc, syc)

      deallocate(tomedx,tomedy,tomedz)
      return

      end subroutine





      subroutine radec2dy(r,d,r0,d0,dy)
      implicit none

      real*8 r, d
      real*8 r0,d0

      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 yrad
      real*8   x, y, z
      real*8   xx,yy,zz
      real*8 dy, rd2y
      !f2py intent(in) :: r,d,r0,d0
      !f2py intent(out) dy

      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)
      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)

      rrrr = sind0*sinde + cosd0*cosde*cosra

      yrad = (cosd0*sinde-sind0*cosde*cosra)/rrrr
      rd2y = yrad*180.0d0/3.141592654d0

      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0

      if (x*xx + y*yy + z*zz.lt.0) rd2y = 90

      dy = rd2y


      return
      end



      subroutine radec2dx(r,d,r0,d0,dx)
      implicit none

      real*8 r, d
      real*8 r0,d0

      real*8 cosra, sinra
      real*8 cosde, sinde
      real*8 cosd0, sind0
      real*8 rrrr
      real*8 xrad
      real*8   x, y, z
      real*8   xx,yy,zz
      real*8 dx, rd2x
      !f2py intent(in) :: r,d,r0,d0
      !f2py intent(out) dx

      cosra = cos((r-r0)*3.141592654d0/180.0d0)
      sinra = sin((r-r0)*3.141592654d0/180.0d0)

      cosde = cos(d *3.141592654d0/180.0d0)
      sinde = sin(d *3.141592654d0/180.0d0)
      cosd0 = cos(d0*3.141592654d0/180.0d0)
      sind0 = sin(d0*3.141592654d0/180.0d0)
  
      rrrr = sind0*sinde + cosd0*cosde*cosra


      xrad = cosde*sinra/rrrr
      rd2x = xrad*180.0d0/3.141592654d0
  
      x  = cosde*cos(r *3.14159/180)
      y  = cosde*sin(r *3.14159/180)
      z  = sinde
      xx = cosd0*cos(r0*3.14159/180)
      yy = cosd0*sin(r0*3.14159/180)
      zz = sind0

      if (x*xx + y*yy + z*zz.lt.0) rd2x = 90

      dx = rd2x
      
      return
      end






















      subroutine barsigw(xl,wl,NTOT,bar,sig)
      implicit none

      integer NTOT
      real*4 xl(NTOT)
      real*4 wl(NTOT)
      real*4 bar
      real*4 sig

      integer n
      real*4  bsum, ssum
      real*4  wsum
      integer NIT
      real*4 sig0


      bar = 0.e0
      sig = 9e9
      NIT = 0
  1   continue
      NIT = NIT + 1
      sig0 = sig
      bsum = 0.
      ssum = 0.
      wsum = 0.
      do n = 1, NTOT
         if (abs(xl(n)-bar).le.2.50*sig) then
            bsum = bsum + wl(n)*xl(n)
            ssum = ssum + wl(n)*abs(xl(n)-bar)
            wsum = wsum + wl(n)
            endif
         enddo
      if (wsum.gt.0) then
         bar = bsum/wsum
         sig = ssum/wsum
         endif
      if (abs(sig-sig0).gt.0.10*sig0.and.NIT.lt.20) goto 1

      return
      end


      subroutine barsig(xlist,NTOT,bar,sig)
      implicit none
 
      integer NTOT,NUSE
      real xlist(NTOT)
      real bar
      real sig
 
      integer n
      real*4  bsum, ssum
      integer NIT
      integer NSUM
      real sig0, bar0
 
      bar = 0.e0
      sig = 9e9
      NIT = 0
  1   continue
      NIT = NIT + 1
      sig0 = sig
      bar0 = bar
      bsum = 0.
      ssum = 0.
      nsum = 0
      do n = 1, NTOT
         if (abs(xlist(n)-bar).le.3.00*sig) then
            bsum = bsum + xlist(n)
            ssum = ssum + abs(xlist(n)-bar)
            nsum = nsum + 1
            endif
         enddo
      if (nsum.gt.0) bar = bsum/nsum
      if (nsum.gt.1) sig = ssum/nsum*1.10
      if (NIT.le.2) goto 1
      if (abs(sig-sig0).gt.0.05*sig0.and.NIT.lt.6) goto 1
      if (abs(bar-bar0).gt.0.50*sig0.and.NIT.lt.6) goto 1
      if (nsum.le.1) sig = 0.999
      NUSE = NSUM
 
      return
      end
 

      

      real function mbar_sky(ixc,iyc,inr,ior,pixarr,NX,NY)
      implicit none

      integer NX, NY
      integer ixc, iyc
      integer inr, ior
      real*4 pixarr(NX,NY)
c      integer*2 pixarr(_NPX_,_NPY_)

      integer i  , j
      integer ixh, iyh
      real    rij
      integer nuse

      real sklist(999)
      integer nsk
      real bar, sig    

      integer HIFLAG
!      common / HIFLAG_ / HIFLAG

      integer LOFLAG
!      common / LOFLAG_ / LOFLAG


      HIFLAG = 80000
      LOFLAG = -100

      nsk = 0
      do i = -ior+1, ior-1
      do j = -ior+1, ior-1
         rij = i**2+j**2
         if (rij.ge.ior**2) goto 333
         if (rij.lt.inr**2) goto 333
         ixh = ixc + i
         iyh = iyc + j
         if (ixh.lt.1.or.ixh.gt.(NX)) goto 333
         if (iyh.lt.1.or.iyh.gt.(NY)) goto 333

         if (pixarr(ixh,iyh).le.LOFLAG) goto 333
         if (pixarr(ixh,iyh).ge.HIFLAG) goto 333
         if (pixarr(ixh,iyh).le.  -250) goto 333
         if (pixarr(ixh,iyh).ge.001000) goto 333
         nsk = nsk + 1
         if (nsk.gt.99999) goto 333
         sklist(nsk) = pixarr(ixh,iyh)
 333     continue
         enddo 
         enddo 

      call barsigg2(sklist,nsk,bar,sig,nuse)      

      mbar_sky = bar
 
      return
      end



      subroutine barsigg2(xlist,NTOT,bar,sig,NUSE)
      implicit none
    
      integer NTOT
      real xlist(NTOT)
      real bar 
      real sig
      integer NUSE

      integer n
      real*4    bsum, ssum
      integer nsum
      integer NIT


      bar = 0.e0
      sig = 9e9
      do NIT = 1, 20
         bsum = 0.
         ssum = 0.
         nsum = 0.
         do n = 1, NTOT
            if (abs(xlist(n)-bar).le.1.75*sig) then
               bsum = bsum + xlist(n)
               ssum = ssum + abs(xlist(n)-bar)
               nsum = nsum + 1
               endif
            enddo
         if (nsum.gt.0) bar = bsum / nsum 
         if (nsum.gt.1) sig = ssum/(nsum-1)
         if (nsum.lt.0.35*NTOT) return
         enddo
      NUSE = nsum
      if (nsum.le.1) sig = 0.999 

      return
      end
 

