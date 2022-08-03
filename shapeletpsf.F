

c read in a PSF map file. Coefficients go into array psf. First index of
c psf() gives spatial variation; second index is shapelet
c coefficient. At position x,y in the image the shapelet model for the
c PSF is
c
c sum_{ik} PSF(i,k) SPACE(i) SHAPE(k)
c
c SPACE(i) are terms of the form X^aY^b in the order
c (ab=00, 10,01, 20,11,02, 30,21,12,03, 40,...) up to a+b=nfit, 
c where X =(2*x/xmax-1), Y=(2*y/ymax-1)
c
c SHAPE(k) are normalized shapelets S_{ab}(dx,dy,beta) in the same order
c as above (ab=00,10,01,20,11,02,30,21,12,03,40,...) up to a+b=nsh.
c
c Format of psfmap file:
c
c  nsh, beta   ! shapelet order, scale radius
c  nfit,xmax,ymax  ! spatial variation polynomial order, scaling for x and y
c  <blank line>
c and then for each shapelet coefficient k, i.e. (nsh+1)*(nsh+2)/2 times:
c  - a line with the (nfit+1)*(nfit+2)/2 polynomial components psf(i,k)  
c  - followed by a blank line

      logical function readpsfmap(filename,psf,nsh,beta,nfit,xmax,ymax)
c read in a shapelet psfmap called 'filename'.
c argument is TRUE on success, FALSE on failure
      character *(*) filename
      real psf(mmpoly2,msh-2)
      readpsfmap=.false.
      open(23,file=filename,status='old',err=2)
c      write(*,*) filename,'opened.'
 1    read(23,*,err=1) nsh,beta
      read(23,*) nfit,xmax,ymax
      read(23,*) 
      do k=1,(nsh+1)*(nsh+2)/2
         read(23,*,err=2) (psf(kfit,k),kfit=1,(nfit+1)*(nfit+2)/2)
         read(23,*)
      enddo
      readpsfmap=.true.
 2    return
      end

c---------------------------------------------------------------
c interpolate psf map coefs coef() to spatial position x,y
c use scale factors xm,ym for field position
c n,beta are shapelet order and scale radius
c nfit is spatial polynomial order of coef()
c sh is the resulting shapelet array

      subroutine interpolpsf(x,y,xm,ym,coef,n,beta,nfit,sh)
      real coef(mmpoly2,msh-2),sh(msh)
      sh(1)=beta
      sh(2)=n
      do k=1,(n+1)*(n+2)/2
          sum=0
          kfit=0
          do i=0,nfit
              do m=0,i
                  l=i-m
                  kfit=kfit+1
                  sum=sum+coef(kfit,k)*polymodel(x,y,xm,ym,l,m)
              enddo
          enddo
          sh(k+2)=sum
      enddo
      return
      end

      function polymodel(x,y,xm,ym,l,m)
c simple power of x and y, with range 0..xm remapped to -1..1
      if (l.eq.0) then
          px=1
      else
          px=(2*x/xm-1)**l
      endif
      if (m.eq.0) then
          py=1
      else
          py=(2*y/ym-1)**m
      endif
      polymodel=px * py
      return
      end


c---------------------------------------------------------
c evaluate shapelet model sh at position x,y with respect to the center
c of the parent gaussian

      function shapeletpix(x,y,sh)
      real sh(msh),xtab(mm+1),ytab(mm+1)
      beta=sh(1)
      n=sh(2)
      k=2
      xnorm=x/beta
      call htab(xnorm,n,xtab)
      ex=exp(-0.5*xnorm**2)
      ynorm=y/beta
      call htab(ynorm,n,ytab)
      ey=exp(-0.5*ynorm**2)
      zz=0
      fac=ex*ey/sqrt(3.1415926535)/beta
      do nn=0,n
         do j=0,nn
            i=nn-j
            k=k+1
            if (j.ne.0) then
               fac=fac/sqrt(j/real(i+1))
            else 
               if (i.ne.0) fac=fac/sqrt(real(i))
            endif
            zz=zz+fac*sh(k)*xtab(i+1)*ytab(j+1)
         enddo
         fac=fac/sqrt(2.)
      enddo
      shapeletpix=zz
      return
      end

      subroutine htab(x,nn,t)
c make table t of Hermite polynomial evaluated at x, for order 0..n
      real t(nn+1)
      t(1)=1
      if (nn.gt.0) t(2)=2*x
c Hn is stored in t(n+1)
c recurrence Hn+1 = 2x Hn - 2n Hn-1
      do n=1,nn-1
          i=n+1
          t(n+2)=2*x*t(n+1) - 2*n*t(n)
      enddo
c      write(0,*) 'HTAB: ',x,n
c      write(0,*) (t(i),i=1,n+1)
      return
      end



