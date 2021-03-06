      epsimin=-1.d0
      epsimax=1.d0
      depsi=(epsimax - epsimin)/dfloat(Lepsi)
      do i=0,Lepsi
         epsi(i)=epsimin+depsi*dfloat(i)
      enddo

      open(40,file='ddArpes.analytic')
      open(41,file='ppArpes.analytic')
      open(50,file='nkdVSepsi.analytic')
      open(51,file='nkpVSepsi.analytic')
      open(52,file='nktotVSepsi.analytic')
      kplus=0.d0
      
      omestep=(omemax-omemin)/dfloat(Lepsi)
      do k=0,Lepsi
         kk=mod(k,20)
         if(kk.eq.0.d0)kplus=kplus+0.1d0
         znd=0.d0
         znp=0.d0
         zntot=0.d0
         do i=0,Lepsi
            omega=omemin+omestep*dfloat(i) !nu(i)
            zr=omega+xi*0.02d0
            alpha=zr+xmu-ed0
            sigmapp=V**2/(alpha)
            Gminus=one/(zr+xmu-ep0-epsi(k)-sigmapp)
            Gd=one/(zr+xmu-ed0-V**2/(zr+xmu-ep0-epsi(k)))
            dosd=-dimag(Gd)/pi
            dosminus=-dimag(Gminus)/pi
            if(kk.eq.0.d0.and.abs(omega).le.omemax)then
               write(40,*)omega,dosd+kplus
               write(41,*)omega,dosminus+kplus
            endif
            xfermi=fermi(omega,0.d0,beta)
            znp=znp+1.d0*dosminus*omestep*xfermi
            znd=znd+1.d0*dosd*omestep*xfermi
            zntot=znd+znp
         enddo
         write(40,*)''
         write(41,*)''
         write(50,*)epsi(k),znd
         write(51,*)epsi(k),znp
         write(52,*)epsi(k),zntot
      enddo

