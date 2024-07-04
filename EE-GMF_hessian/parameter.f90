       program parameter
       implicit none
      
       integer(kind=4)::i,j,k,natm,k1,k2,k3,atomnum(999),resnum(999),&
                        sub,markC(999),markN(999),markCA(999),ini,nc,&
                        nd,ii,jj,nii,step,na,nai,naj,ist,mark1,mark2,&
                        n1,n2,nxm,nym,nzm,nxl,nyl,nzl,naij,ic,jc,kc,&
                        nres,iii,jjj,kkk,groupn(0:50)
       real(kind=8)::f(500,500),ftemp(1000000),&
                     dtemp(1000000),ptemp(1000000),&
                     d(3,9999),p(6,9999),f1(999,999),dist,charge(999),&
                     coord(999,3),fmm(999,999),clata,clatb,clatc,xa,xb,&
                     ya,yb,za,zb,temp(3)
       character(len=4)::atomname(999)
       character(len=3)::resname(999)
       character(len=200)::fname
       character(len=200)::fline
       character(len=200)::cmdstr
       logical::connect(50,50),existence

       i=0
       groupn=0
       open(5,file='cluster.pdb')
       open(6,file='mass')
       do while (.true.)
       read(5,'(A80)',end=100)fline
       if(fline(1:4).eq.'ATOM')then
        i=i+1
        read(fline,1000)atomnum(i),atomname(i),resname(i),resnum(i)
        groupn(resnum(i))=groupn(resnum(i))+1
       endif
       enddo
100    close(5)
1000   format(7x,I4,1x,a4,1x,a3,2x,I4)
         natm=i
         write(6,*)natm
       do i=1,natm
         if(atomname(i)(2:2).eq.'N')then
          write(6,'(A10)') '  14.00307'
         endif
         if(atomname(i)(2:2).eq.'H')then
          write(6,'(A10)') '   1.00783'
         endif
         if(atomname(i)(2:2).eq.'C')then
          write(6,'(A10)') '  12.00000'
         endif
         if(atomname(i)(2:2).eq.'O')then
          write(6,'(A10)') '  15.99491'
         endif
         if(atomname(i)(2:2).eq.'S')then
          write(6,'(A10)') '  31.97207'
         endif
       enddo
       close(6)
       
      open(5,file='min.rst')
      read(5,*)
      read(5,*)
      read(5,'(6f15.8)')((coord(i,j),j=1,3),i=1,natm)
      close(5)

         nres=resnum(natm)

      do i=1,nres
       groupn(i)=groupn(i)+groupn(i-1)
      enddo

1001  format(5E16.8)

      do i=1,natm
       if((resname(i).eq.'H3O').and.(atomname(i).eq.' O  '))then
        charge(i)=-0.960222
       endif
       if((resname(i).eq.'H3O').and.(atomname(i)(2:2).eq.'H'))then
        charge(i)=0.653407
       endif
       if((resname(i).eq.'H2O').and.(atomname(i).eq.' O  '))then
        charge(i)=-0.820d0
       endif
       if((resname(i).eq.'H2O').and.(atomname(i)(2:2).eq.'H'))then
        charge(i)=0.410d0
       endif
      enddo

!======================chain fragmentation=================================
       f(:,:)=0.0d0
       d(:,:)=0.0d0
       p(:,:)=0.0d0
       do i=1,nres  
        fname='000cell'//char(48+(i)/100)//&
              char(48+((i)-(i)/100*100)/10)//&
              char(48+((i)-(i)/100*100-((i)-(i)/100*100)/10*10))
        cmdstr='formchk '//trim(adjustl(fname))//'.chk '&
               //trim(adjustl(fname))//'.fchk'
!       call system(cmdstr)
        existence=.false.
        open(8,file=trim(adjustl(fname))//'.fchk')
        do while (.true.)
         read(8,'(A80)',iostat=ist)fline
         if(fline(1:25).eq.'Cartesian Force Constants')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ftemp(j),j=1,sub)
          existence=.true.
         endif
         if(fline(1:18).eq.'Dipole Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(dtemp(j),j=1,sub)
          na=sub/9
         endif
         if(fline(1:26).eq.'Polarizability Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ptemp(j),j=1,sub)
         endif
         if(ist/=0)exit
        enddo
        close(8)
        if(.not.existence)then
         write(*,*)"Cartesian Force Constants Not Found"
         stop
        endif
           do j=1,na*3
            do k=1,j
             f1(j,k)=ftemp(j*(j-1)/2+k)
            enddo
           enddo
           step=groupn(i-1)*3
           do j=1,na*3
            do k=1,j
              f(j+step,k+step)=f(j+step,k+step)+f1(j,k)
              f(k+step,j+step)=f(j+step,k+step)
            enddo
           enddo
           do j=1,3
            do k=1,na*3
             d(j,k+step)=d(j,k+step)+dtemp((k-1)*3+j)
            enddo
           enddo
           do j=1,6
            do k=1,na*3 
             p(j,k+step)=p(j,k+step)+ptemp((k-1)*6+j)
            enddo
           enddo
       enddo
!=========================two body parts===============================
      connect=.false.
      do i=1,nres-1
       do j=i+1,nres
        do ii=1,natm
        if((resnum(ii)==i).and.(atomname(ii)(2:2).eq.'O'))then
         do jj=1,natm
         if((resnum(jj)==j).and.(atomname(jj)(2:2).eq.'O'))then
       dist=sqrt((coord(ii,1)-coord(jj,1))**2+&
                 (coord(ii,2)-coord(jj,2))**2+&
                 (coord(ii,3)-coord(jj,3))**2)
          if(dist.le.99.0d0)then
           connect(i,j)=.true.
          endif
         endif
         enddo
        endif
        enddo
       enddo
      enddo

      do ii=1,nres-1
       do jj=ii+1,nres
        if(connect(ii,jj))then

        fname='000cell'//char(48+ii/100)//&
               char(48+(ii-ii/100*100)/10)//&
               char(48+(ii-ii/10*10))
!       cmdstr='formchk '//trim(adjustl(fname))//'.chk '&
!               //trim(adjustl(fname))//'.fchk'
!       call system(cmdstr)
        open(8,file=trim(adjustl(fname))//'.fchk')
        do while (.true.)
         read(8,'(A80)',iostat=ist)fline
         if(fline(1:25).eq.'Cartesian Force Constants')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ftemp(j),j=1,sub)
         endif
         if(fline(1:18).eq.'Dipole Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(dtemp(j),j=1,sub)
          nai=sub/9
         endif
         if(fline(1:26).eq.'Polarizability Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ptemp(j),j=1,sub)
         endif
         if(ist/=0)exit
        enddo
        close(8)
           do j=1,nai*3
            do k=1,j
             f1(j,k)=ftemp(j*(j-1)/2+k)
            enddo
           enddo
           step=groupn(ii-1)*3
           do j=1,nai*3
            do k=1,j
              f(j+step,k+step)=f(j+step,k+step)-f1(j,k)
              f(k+step,j+step)=f(j+step,k+step)
            enddo
           enddo
           do j=1,3
            do k=1,nai*3
             d(j,k+step)=d(j,k+step)-dtemp((k-1)*3+j)
            enddo
           enddo
           do j=1,6
            do k=1,nai*3
             p(j,k+step)=p(j,k+step)-ptemp((k-1)*6+j)
            enddo
           enddo
         
        fname='000cell'//char(48+jj/100)//&
               char(48+(jj-jj/100*100)/10)//&
               char(48+(jj-jj/10*10))
!       cmdstr='formchk '//trim(adjustl(fname))//'.chk '&
!              //trim(adjustl(fname))//'.fchk'
!       call system(cmdstr)
        open(8,file=trim(adjustl(fname))//'.fchk')
        do while (.true.)
         read(8,'(A80)',iostat=ist)fline
         if(fline(1:25).eq.'Cartesian Force Constants')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ftemp(j),j=1,sub)
         endif
         if(fline(1:18).eq.'Dipole Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(dtemp(j),j=1,sub)
          naj=sub/9
         endif
         if(fline(1:26).eq.'Polarizability Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ptemp(j),j=1,sub)
         endif
         if(ist/=0)exit
        enddo
        close(8)
           do j=1,naj*3
            do k=1,j
             f1(j,k)=ftemp(j*(j-1)/2+k)
            enddo
           enddo
           step=groupn(jj-1)*3
           do j=1,naj*3
            do k=1,j
              f(j+step,k+step)=f(j+step,k+step)-f1(j,k)
              f(k+step,j+step)=f(j+step,k+step)
            enddo
           enddo
           do j=1,3
            do k=1,naj*3
             d(j,k+step)=d(j,k+step)-dtemp((k-1)*3+j)
            enddo
           enddo
           do j=1,6
            do k=1,naj*3
             p(j,k+step)=p(j,k+step)-ptemp((k-1)*6+j)
            enddo
           enddo   

       fname='000cell'//char(48+ii/100)//&
              char(48+(ii-ii/100*100)/10)//&
              char(48+(ii-ii/10*10))//'-000cell'//char(48+jj/100)//&
              char(48+(jj-jj/100*100)/10)//char(48+(jj-jj/10*10))
       cmdstr='formchk '//trim(adjustl(fname))//'.chk '&
              //trim(adjustl(fname))//'.fchk'
!      call system(cmdstr)
       existence=.false.
        open(8,file=trim(adjustl(fname))//'.fchk')
        do while (.true.)
         read(8,'(A80)',iostat=ist)fline
         if(fline(1:25).eq.'Cartesian Force Constants')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ftemp(j),j=1,sub)
          existence=.true.
         endif
         if(fline(1:18).eq.'Dipole Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(dtemp(j),j=1,sub)
         endif
         if(fline(1:26).eq.'Polarizability Derivatives')then
          read(fline,'(50x,I11)') sub
          read(8,'(5E16.8)')(ptemp(j),j=1,sub)
         endif
         if(ist/=0)exit
        enddo
        close(8)
        if(.not.existence)then
         write(*,*)"Cartesian Force Constants Not Found"
         stop
        endif
           do j=1,nai*3
            do k=1,j
             f1(j,k)=ftemp(j*(j-1)/2+k)
            enddo
           enddo
           step=groupn(ii-1)*3
           do j=1,nai*3
            do k=1,j
              f(j+step,k+step)=f(j+step,k+step)+f1(j,k)
              f(k+step,j+step)=f(j+step,k+step)
            enddo
           enddo
           do j=1,3
            do k=1,nai*3
             d(j,k+step)=d(j,k+step)+dtemp((k-1)*3+j)
            enddo
           enddo
           do j=1,6
            do k=1,nai*3
             p(j,k+step)=p(j,k+step)+ptemp((k-1)*6+j)
            enddo
           enddo

           do j=nai*3+1,(nai+naj)*3
            do k=nai*3+1,j
             f1(j-nai*3,k-nai*3)=ftemp(j*(j-1)/2+k)
            enddo
           enddo
           step=groupn(jj-1)*3
           do j=1,naj*3
            do k=1,j
              f(j+step,k+step)=f(j+step,k+step)+f1(j,k)
              f(k+step,j+step)=f(j+step,k+step)
            enddo
           enddo
           do j=nai*3+1,(nai+naj)*3
            do k=1,nai*3
             f1(j-nai*3,k)=ftemp(j*(j-1)/2+k)
            enddo
           enddo 
           do j=1,naj*3
            do k=1,nai*3
             f(j+step,k+groupn(ii-1)*3)=f(j+step,k+groupn(ii-1)*3)+&
                                          f1(j,k)
             f(k+groupn(ii-1)*3,j+step)=f(j+step,k+groupn(ii-1)*3)
            enddo
           enddo
           do j=1,3
            do k=1,naj*3
             d(j,k+step)=d(j,k+step)+dtemp((k-1)*3+j+3*nai*3)
            enddo
           enddo
           do j=1,6
            do k=1,naj*3
             p(j,k+step)=p(j,k+step)+ptemp((k-1)*6+j+6*nai*3)
            enddo
           enddo

        else

         do i=1,natm
         if(resnum(i)==ii)then
          do j=1,natm
          if(resnum(j)==jj)then
           dist=sqrt((coord(i,1)-coord(j,1))**2+&
                     (coord(i,2)-coord(j,2))**2+&
                     (coord(i,3)-coord(j,3))**2)
           do k1=1,3
            do k2=1,3
             if(k1==k2)then
      f(k1+j*3-3,k2+i*3-3)=f(k1+j*3-3,k2+i*3-3)+&
      ((0.5291772083)**3)*(dist**2-3.0d0*(coord(i,k1)-coord(j,k2))**2)*&
                     charge(i)*charge(j)/(dist**5)
      f(k2+i*3-3,k1+j*3-3)=f(k1+j*3-3,k2+i*3-3)
             else
      f(k1+j*3-3,k2+i*3-3)=f(k1+j*3-3,k2+i*3-3)+&
                ((0.5291772083)**3)*(-3.0d0)*(coord(i,k1)-coord(j,k1))*&
                                (coord(i,k2)-coord(j,k2))*&
                  charge(i)*charge(j)/(dist**5)
      f(k2+i*3-3,k1+j*3-3)=f(k1+j*3-3,k2+i*3-3)
             endif
            enddo
           enddo
          endif 
          enddo
         endif
         enddo

        endif
       enddo
      enddo

       open(222,file='dipole_deriv.info')
       open(333,file='raman.txt')
       open(444,file='hessian.info')
       write(222,'(A)')'Dipole Derivative'
       write(222,'(I6)')natm*3*3
       write(222,'(5ES16.8)')((d(i,j),i=1,3),j=1,natm*3)
       write(333,'(5ES16.8)')((p(i,j),i=1,6),j=1,natm*3)
       write(444,'(A)')'Hessian'
       write(444,'(I6)')(natm*3+1)*natm*3/2 
       write(444,'(5ES16.8)')((f(i,j),j=1,i),i=1,natm*3)
       close(222)
       close(333)
       close(444)

       end program
















