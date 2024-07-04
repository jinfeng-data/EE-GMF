              subroutine energy_cal
              use variables
              implicit none
              include 'parallel.h'
              include 'mpif.h'

              integer(kind=4)::niter
              real(kind=8)::maxx,maxy,maxz,enthalpy
              real(kind=8)::gr(natom*3+3)
              logical::ccsd

      energy2=0.0
      tenergy2=0.0
      efragment=0.0
      tefragment=0.0
      etotal=0.0
      force=0.0
      tforce=0.0
      pforce=0.0
      pforce_global=0.0
      field=0.0
      tfield=0.0
      pfield=0.0
      pfield_global=0.0

      charge=charge*18.2223
      njob=0

       do i=1,nres
        njob=njob+1
        cmd(njob)='g16 000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo

       do i=1,nres-1
        do j=i+1,nres
         if(connect(i,j))then
          njob=njob+1
          cmd(njob)='g16 000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.gjf'
         endif
        enddo
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
       do i=mytid+1,njob,numprocs
        call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       njob=0

       do i=1,nres
        njob=njob+1
        cmd(njob)='g16 field_000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo

       do i=1,nres-1
        do j=i+1,nres
         if(connect(i,j))then
          njob=njob+1
          cmd(njob)='g16 field_000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.gjf'
         endif
        enddo
       enddo
      
       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       do i=mytid+1,njob,numprocs
        call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       njob=0

       do i=1,nres
        njob=njob+1
        cmd(njob)='formchk 000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.chk'
       enddo

       do i=1,nres-1
        do j=i+1,nres
         if(connect(i,j))then
          njob=njob+1
          cmd(njob)='formchk 000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.chk'
         endif
        enddo
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       do i=mytid+1,njob,numprocs
        call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

!====================== QM energy calculation=====================
      njob=0
      do i=1,nres
       njob=njob+1
       cmd(njob)='000cell'//char(48+i/100)//&
                   char(48+(i-i/100*100)/10)//&
                   char(48+(i-i/10*10))//'.log'
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).eq.'000')then

2017  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3)')l
      normal=.false.
      ccsd=.false.
      do while(.true.)
      read(8,'(a80)',iostat=ist)fline
      if(ist/=0)exit
      if((fline(2:9).eq.'CCSD(T)=').and.(.not.ccsd))then
!     indexx=index(fline,'=')+2
!     read(fline(indexx:indexx+18),*)efragment(i)
      read(fline(11:30),*)tefragment(l)
      tefragment(l)=tefragment(l)-selfeng
      ccsd=.true.
      endif
      if(fline(28:35).eq.'E(CORR)=')then
      read(fline(37:53),*)tefragment(l)
      tefragment(l)=tefragment(l)-selfeng
      endif
      if(fline(2:27).eq.'Self energy of the charges')then
      read(fline(30:50),*)selfeng
      endif
!     if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
!        read(8,*)
!        read(8,*)
!        ic=0
!        do m=1,natom
!        if(resnum(m)==l)then
!        read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!        tforce(m,1)=tforce(m,1)+temp(1)
!        tforce(m,2)=tforce(m,2)+temp(2)
!        tforce(m,3)=tforce(m,3)+temp(3)
!        ic=ic+1
!        pforce(l,ic,1)=temp(1)
!        pforce(l,ic,2)=temp(2)
!        pforce(l,ic,3)=temp(3)
!        endif
!        enddo
!        normal=.true.
!        exit
!     endif
      enddo
      close(8)
      open(9,file=cmd(i)(1:10)//'.fchk')
       do while(.true.)
       read(9,'(a80)')fline
       if(fline(1:18).eq.'Cartesian Gradient')then
        read(fline(54:61),'(I8)')n
        read(9,'(5E16.8)')(temp(k),k=1,n)
        temp=-temp
        ic=0
        do m=1,natom       
        if(resnum(m)==l)then
        ic=ic+1
        tforce(m,1)=tforce(m,1)+temp(ic*3-2)
        tforce(m,2)=tforce(m,2)+temp(ic*3-1)
        tforce(m,3)=tforce(m,3)+temp(ic*3)
        pforce(l,ic,1)=temp(ic*3-2)
        pforce(l,ic,2)=temp(ic*3-1)
        pforce(l,ic,3)=temp(ic*3)                                                                
        endif
        enddo
        exit
       endif
       enddo
      close(9)
!     if(.not.normal)then
!      k=len(trim(adjustl(cmd(i))))
!      call system('g16 '//cmd(i)(1:k-3)//'gjf')
!      goto 2017
!     endif

1990   open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
       normal=.false.
       do while(.true.)
        read(6,'(a70)',iostat=ist)fline                                                 
        if(ist/=0)exit
        if(fline(42:55).eq.'Electric Field')then
        read(6,*)
        read(6,*)
        do m=1,natom
         if(resnum(m)==l)then
          read(6,*)
         endif
        enddo
        do m=1,natom
         if(resnum(m)/=l)then
        read(6,'(28x,3f18.10)')(temp(k),k=1,3)
        tfield(m,1)=tfield(m,1)+temp(1)
        tfield(m,2)=tfield(m,2)+temp(2)
        tfield(m,3)=tfield(m,3)+temp(3)
        pfield(l,m,1)=temp(1)
        pfield(l,m,2)=temp(2)
        pfield(l,m,3)=temp(3)
         endif
        enddo
        normal=.true.
        exit
        endif
       enddo
       close(6)
!      if(.not.normal)then
!       k=len(trim(adjustl(cmd(i))))
!       call system('g16 field_'//cmd(i)(1:k-3)//'gjf')
!       goto 1990
!      endif

      endif
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(pforce,pforce_global,nres*50*3,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(pfield,pfield_global,nres*natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      
      njob=0
      do i=1,nres-1
       do j=i+1,nres
        if(connect(i,j))then
         njob=njob+1
         cmd(njob)='000cell'//char(48+i/100)//&
                char(48+(i-i/100*100)/10)//&
                char(48+(i-i/10*10))//&
                '-000cell'//char(48+j/100)//&
                char(48+(j-j/100*100)/10)//&
                char(48+(j-j/10*10))//'.log'
        endif
       enddo
      enddo

      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).eq.'000')then

2019  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3,8x,I3)')j,l
      normal=.false.
      ccsd=.false.
          do while(.true.)
          read(8,'(a80)',iostat=ist)fline
          if(ist/=0)exit 
          if((fline(2:9).eq.'CCSD(T)=').and.(.not.ccsd))then
!          indexx=index(fline,'=')+2                                            
!          read(fline(indexx:indexx+18),*)energy2(i,j)    
           read(fline(11:30),*)tenergy2(j,l)
           tenergy2(j,l)=tenergy2(j,l)-selfeng
           ccsd=.true.
          endif
          if(fline(28:35).eq.'E(CORR)=')then
           read(fline(37:53),*)tenergy2(j,l) 
           tenergy2(j,l)=tenergy2(j,l)-selfeng
          endif
          if(fline(2:27).eq.'Self energy of the charges')then  
           read(fline(30:50),*)selfeng                          
          endif
!         if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
!          read(8,*)  
!          read(8,*)  
!          ic=0
!          do m=1,natom
!          if(resnum(m)==j)then
!          read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!          ic=ic+1
!          tforce(m,1)=tforce(m,1)+temp(1)-pforce_global(j,ic,1)
!          tforce(m,2)=tforce(m,2)+temp(2)-pforce_global(j,ic,2)
!          tforce(m,3)=tforce(m,3)+temp(3)-pforce_global(j,ic,3)
!          endif
!          enddo
!          ic=0
!          do m=1,natom
!          if(resnum(m)==l)then
!          ic=ic+1
!          read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!          tforce(m,1)=tforce(m,1)+temp(1)-pforce_global(l,ic,1)
!          tforce(m,2)=tforce(m,2)+temp(2)-pforce_global(l,ic,2)
!          tforce(m,3)=tforce(m,3)+temp(3)-pforce_global(l,ic,3)
!          endif
!          enddo
!          normal=.true.
!          exit
!         endif
          enddo
         close(8)
      open(9,file=cmd(i)(1:21)//'.fchk')
       do while(.true.)
       read(9,'(a80)')fline
       if(fline(1:18).eq.'Cartesian Gradient')then                                             
        read(fline(54:61),'(I8)')n
        read(9,'(5E16.8)')(temp(k),k=1,n)                                                      
        temp=-temp                                                                             
        ic=0
        ii=0
        jj=0
        do m=1,natom       
        if(resnum(m)==j)then
        ic=ic+1
        ii=ii+1
        tforce(m,1)=tforce(m,1)+temp(ic*3-2)-pforce_global(j,ii,1)
        tforce(m,2)=tforce(m,2)+temp(ic*3-1)-pforce_global(j,ii,2)
        tforce(m,3)=tforce(m,3)+temp(ic*3)-pforce_global(j,ii,3)
        endif
        enddo
        do m=1,natom       
        if(resnum(m)==l)then
        ic=ic+1
        jj=jj+1
        tforce(m,1)=tforce(m,1)+temp(ic*3-2)-pforce_global(l,jj,1)
        tforce(m,2)=tforce(m,2)+temp(ic*3-1)-pforce_global(l,jj,2)
        tforce(m,3)=tforce(m,3)+temp(ic*3)-pforce_global(l,jj,3)
        endif
        enddo
        exit
       endif
       enddo
      close(9)
!     if(.not.normal)then
!      k=len(trim(adjustl(cmd(i))))
!      call system('g16 '//cmd(i)(1:k-3)//'gjf')
!      goto 2019
!     endif

1991     open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
         normal=.false.
         do while(.true.)
           read(6,'(a70)',iostat=ist)fline
           if(ist/=0)exit
           if(fline(42:55).eq.'Electric Field')then
          read(6,*) 
          read(6,*)
           do m=1,natom
            if((resnum(m)==j).or.(resnum(m)==l))then
             read(6,*)
            endif
           enddo
           do m=1,natom
           if((resnum(m)/=j).and.&
              (resnum(m)/=l))then
           read(6,'(28x,3f18.10)')(temp(k),k=1,3)
           tfield(m,1)=tfield(m,1)+temp(1)
           tfield(m,2)=tfield(m,2)+temp(2)
           tfield(m,3)=tfield(m,3)+temp(3)
           endif
           enddo
           do m=1,natom
           if(resnum(m)/=j)then
          tfield(m,1)=tfield(m,1)-pfield_global(j,m,1)
          tfield(m,2)=tfield(m,2)-pfield_global(j,m,2)
          tfield(m,3)=tfield(m,3)-pfield_global(j,m,3)
           endif
           enddo
           do m=1,natom
           if(resnum(m)/=l)then
          tfield(m,1)=tfield(m,1)-pfield_global(l,m,1)
          tfield(m,2)=tfield(m,2)-pfield_global(l,m,2)
          tfield(m,3)=tfield(m,3)-pfield_global(l,m,3)
           endif
           enddo
           normal=.true.
           exit
           endif
         enddo
         close(6)
!        if(.not.normal)then
!         k=len(trim(adjustl(cmd(i))))
!         call system('g16 field_'//cmd(i)(1:k-3)//'gjf')
!         goto 1991
!        endif

      endif
      enddo
      
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_Allreduce(tforce,force,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tfield,field,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tefragment,efragment,nres,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tenergy2,energy2,nres*nres,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      
      if(master)then
        
          fele=0.0
          do i=1,nres-1
           do j=i+1,nres
            if(.not.connect(i,j))then
             do k=1,natom
             if(resnum(k)==i)then
              do l=1,natom
              if(resnum(l)==j)then
            xa=coord(k,1)
            ya=coord(k,2)
            za=coord(k,3)
            xc=coord(l,1)
            yc=coord(l,2)
            zc=coord(l,3)
            dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
            fele(k,1)=fele(k,1)+charge(k)*charge(l)/(dist**3)*(xa-xc)
            fele(k,2)=fele(k,2)+charge(k)*charge(l)/(dist**3)*(ya-yc)
            fele(k,3)=fele(k,3)+charge(k)*charge(l)/(dist**3)*(za-zc)
            fele(l,1)=fele(l,1)-charge(k)*charge(l)/(dist**3)*(xa-xc)
            fele(l,2)=fele(l,2)-charge(k)*charge(l)/(dist**3)*(ya-yc)
            fele(l,3)=fele(l,3)-charge(k)*charge(l)/(dist**3)*(za-zc)
              endif
              enddo
             endif
             enddo
            endif
           enddo
          enddo  

          fele=fele*0.5291771d0/627.51

      do i=1,natom
      force(i,1)=force(i,1)+field(i,1)*charge(i)/18.2223-fele(i,1)
      force(i,2)=force(i,2)+field(i,2)*charge(i)/18.2223-fele(i,2)
      force(i,3)=force(i,3)+field(i,3)*charge(i)/18.2223-fele(i,3)
      enddo

!==================end QM energy calculation=====================
!==============two body QM energy calculation====================

       open(99,access='append',file='out.data') 
       do i=1,nres
        etotal=etotal+efragment(i)
       enddo
       
2000   format(a11,I4,1x,a12,f16.8)
2001   format(3I2,1x,a4,I4,1x,a12,f16.8)

      energy2body=0.0
      
      do i=1,nres-1
       do j=i+1,nres
        if(connect(i,j))then
        energy2body=energy2body+energy2(i,j)-efragment(i)-&
                    efragment(j)  
        else
         ele2=0.0
         do ii=1,natom
         if(resnum(ii)==i)then
          do jj=1,natom
          if(resnum(jj)==j)then
           dist=sqrt((coord(ii,1)-coord(jj,1))**2+&
                     (coord(ii,2)-coord(jj,2))**2+&
                     (coord(ii,3)-coord(jj,3))**2)
           ele2=ele2+charge(ii)*charge(jj)/dist
          endif
          enddo
         endif 
         enddo
          energy2body=energy2body-ele2/627.51
        endif
       enddo
      enddo

      write(99,*)'The self-energy of the unit cell (au) '//&
                 '(without correction):',etotal+energy2body

!=========================Long range interaction===================================     

2011  format(a11,1x,2I4,a20,2f16.10)
2012  format(3I2,1x,a4,1x,2I4,a20,2f16.10)
      
      etotal=etotal+energy2body

!==============================================================

      maxx=0.0
      maxy=0.0
      maxz=0.0
     open(001,file='force.dat')
      write(001,'(f16.9)')etotal
      do i=1,natom
       write(001,'(I5,3x,3E18.9)')i,(force(i,k),k=1,3)
       if(abs(force(i,1))>maxx)then
         maxx=abs(force(i,1))
       endif
       if(abs(force(i,2))>maxy)then
         maxy=abs(force(i,2))
       endif
       if(abs(force(i,3))>maxz)then
         maxz=abs(force(i,3))
       endif
      enddo
     close(001)

     write(99,*)'Total electrical energy of the unit cell (au):',etotal
     write(99,'(a19,2x,3f10.6)')'The maximal force:',maxx,maxy,maxz
     write(99,*)'======================================================'

       close(99)

      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

       end subroutine
	
