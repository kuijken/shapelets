
      parameter (mstamp=20)
      real psf(mmpoly2,msh-2),sh(msh),stamp(-mstamp:mstamp,-mstamp:mstamp)
      character*132 mapfile
      logical readpsfmap
      
      read (*,'(a132)') mapfile
      if (.not. readpsfmap(mapfile,psf,nsh,beta,nfit,xmax,ymax)) then
        write(0,*) 'Error reading map ',mapfile
        stop
      endif

c read x,y from stdin, compute PSF stamp
1     read(*,*,err=1,end=2) x,y

c     get shapelet model at image position x,y
      call interpolpsf(x,y,xmax,ymax,psf,nsh,beta,nfit,sh)
      do ix=-mstamp,mstamp
       do iy=-mstamp,mstamp
        xx=real(ix)
        yy=real(iy)
        stamp(ix,iy)=shapeletpix(xx,yy,sh)
       enddo
      enddo
      write(*,*) (sh(i),i=1,(nsh+1)*(nsh+2)/2+2)
      write(*,*) x,y,mstamp
      write(*,*) ((stamp(i,j),i=-mstamp,mstamp),j=-mstamp,mstamp)
 2    end

