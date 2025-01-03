              subroutine energy_cal
              use variables
              implicit none
              include 'parallel.h'
              include 'mpif.h'

      real(kind=8)::enthalpy

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
      da=0.0
      db=0.0
      dc=0.0
      xd=0.0
      txd=0.0
      yd=0.0
      tyd=0.0
      zd=0.0
      tzd=0.0
      dip=0.0
      tdip=0.0

      njob=0

       do i=1,ngroup
        njob=njob+1
        cmd(njob)='g16 000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo

        do ic=-nxm,nxm
         do jc=-nym,nym
          do kc=-nzm,nzm
           if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
            do i=1,ngroup
             if(onebody(ic,jc,kc,i))then
              njob=njob+1
              cmd(njob)='g16 '//char(48+ic+nxm*2)//&
        char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
             endif
            enddo
           endif
          enddo
         enddo
        enddo

       do i=1,ngroup-1
        do j=i+1,ngroup
         if(connect(0,0,0,i,j))then
          njob=njob+1
          cmd(njob)='g16 000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.gjf'
         endif
        enddo
       enddo

       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do i=1,ngroup
            do j=1,ngroup
             if(connect(ic,jc,kc,i,j))then
              njob=njob+1
              cmd(njob)='g16 '//char(48+ic+nxm*2)//&
         char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
         char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))//'-000cell'//&
         char(48+j/100)//char(48+(j-j/100*100)/10)//&
                         char(48+(j-j/10*10))//'.gjf'
             endif
            enddo
           enddo
          endif
         enddo
        enddo
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
       do i=mytid+1,njob,numprocs
        call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       njob=0

       do i=1,ngroup
        njob=njob+1
        cmd(njob)='g16 field_000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo
!       do ic=-nxm,nxm
!        do jc=-nym,nym
!         do kc=-nzm,nzm
!          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!           do i=1,nres
!            if(onebody(ic,jc,kc,i))then
!             njob=njob+1
!             cmd(njob)='g16 field_'//char(48+ic+nxm*2)//&
!        char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
!                        char(48+(i-i/100*100)/10)//&
!                        char(48+(i-i/10*10))//'.gjf'
!            endif
!           enddo
!          endif
!         enddo
!        enddo
!       enddo
       do i=1,ngroup-1
        do j=i+1,ngroup
         if(connect(0,0,0,i,j))then
          njob=njob+1
          cmd(njob)='g16 field_000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.gjf'
         endif
        enddo
       enddo
!      do ic=-nxm,nxm
!       do jc=-nym,nym
!        do kc=-nzm,nzm
!         if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!          do i=1,nres
!           do j=1,nres
!            if(connect(ic,jc,kc,i,j))then
!             njob=njob+1
!             cmd(njob)='g16 field_'//char(48+ic+nxm*2)//&
!        char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
!        char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))//'-000cell'//&
!        char(48+j/100)//char(48+(j-j/100*100)/10)//&
!                        char(48+(j-j/10*10))//'.gjf'
!            endif
!           enddo
!          enddo
!         endif
!        enddo
!       enddo
!      enddo
      
       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       do i=mytid+1,njob,numprocs
        call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

!====================== QM energy calculation=====================
      njob=0
      do i=1,ngroup
       njob=njob+1
       cmd(njob)='000cell'//char(48+i/100)//&
                   char(48+(i-i/100*100)/10)//&
                   char(48+(i-i/10*10))//'.log'
      enddo
      do ic=-nxm,nxm
       do jc=-nym,nym
        do kc=-nzm,nzm
         if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
          do i=1,ngroup
           if(onebody(ic,jc,kc,i))then
            njob=njob+1 
            cmd(njob)=char(48+ic+nxm*2)//char(48+jc+nym*2)//&
                     char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
                   char(48+(i-i/100*100)/10)//&
                   char(48+(i-i/10*10))//'.log'
           endif
          enddo
         endif
        enddo
       enddo
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      
      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).eq.'000')then

      bgc=.true.
2017  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3)')l
      normal=.false.
      do while(.true.)
      read(8,'(a80)',iostat=ist)fline
      if(ist/=0)exit
      if(fline(28:32).eq.'EUMP2')then
!     indexx=index(fline,'=')+2
!     read(fline(indexx:indexx+18),*)efragment(i)
      read(fline(36:59),*)tefragment(0,0,0,l)
      tefragment(0,0,0,l)=tefragment(0,0,0,l)-selfeng
      endif
      if(fline(2:27).eq.'Self energy of the charges')then
      read(fline(30:50),*)selfeng
      endif
      if(fline(2:14).eq.'Dipole moment')then
       if(group(l,4).eq.0)then
       read(8,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')txd(l),&
           tyd(l),tzd(l),tdip(l)
       endif
      endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
         read(8,*)
         read(8,*)
         ic=0
         do m=1,20
         if(group(l,m).ne.0)then
         read(8,'(23x,3f15.9)')(temp(k),k=1,3)
         bgc(group(l,m))=.false.
         tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)
         tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)
         tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)
         ic=ic+1
         pforce(0,0,0,l,ic,1)=temp(1)
         pforce(0,0,0,l,ic,2)=temp(2)
         pforce(0,0,0,l,ic,3)=temp(3)
         else
          exit
         endif
         enddo
         normal=.true.
         exit
      endif
      enddo
      close(8)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g09 '//cmd(i)(1:k-3)//'gjf')
       goto 2017
      endif

1990   open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
       normal=.false.
       do while(.true.)
        read(6,'(a70)',iostat=ist)fline                                                 
        if(ist/=0)exit
        if(fline(42:55).eq.'Electric Field')then
        read(6,*)
        read(6,*)
        do m=1,20
         if(group(l,m).ne.0)then
          read(6,*)
         else
          exit
         endif
        enddo
        do m=1,natom
         if(bgc(m))then
        read(6,'(28x,3f18.10)')(temp(k),k=1,3)
        tfield(0,0,0,m,1)=tfield(0,0,0,m,1)+temp(1)
        tfield(0,0,0,m,2)=tfield(0,0,0,m,2)+temp(2)
        tfield(0,0,0,m,3)=tfield(0,0,0,m,3)+temp(3)
        pfield(l,0,0,0,m,1)=temp(1)
        pfield(l,0,0,0,m,2)=temp(2)
        pfield(l,0,0,0,m,3)=temp(3)
         endif
        enddo
!       do ic=-nxm,nxm
!        do jc=-nym,nym
!         do kc=-nzm,nzm
!          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!           do m=1,natom
!       read(6,'(28x,3f18.10)')(temp(k),k=1,3)
!       tfield(ic,jc,kc,m,1)=tfield(ic,jc,kc,m,1)+temp(1)
!       tfield(ic,jc,kc,m,2)=tfield(ic,jc,kc,m,2)+temp(2)
!       tfield(ic,jc,kc,m,3)=tfield(ic,jc,kc,m,3)+temp(3)
!       pfield(l,ic,jc,kc,m,1)=temp(1)
!       pfield(l,ic,jc,kc,m,2)=temp(2)
!       pfield(l,ic,jc,kc,m,3)=temp(3)
!           enddo
!          endif
!         enddo
!        enddo
!       enddo
        normal=.true.
        exit
        endif
       enddo
       close(6)
       if(.not.normal)then
!       call sleep(999999)
        k=len(trim(adjustl(cmd(i))))
        call system('g09 field_'//cmd(i)(1:k-3)//'gjf')
        goto 1990
       endif

      endif
      enddo

      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).ne.'000')then
          
2018  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(3I1,4x,I3)')ic,jc,kc,l
      ic=ic-nxm*2
      jc=jc-nym*2
      kc=kc-nzm*2
      normal=.false.
       do while(.true.)
       read(8,'(a80)',iostat=ist)fline
       if(ist/=0)exit
       if(fline(28:32).eq.'EUMP2')then
!!      indexx=index(fline,'=')+2
!!      read(fline(indexx:indexx+18),*)efragment(i)
        read(fline(36:59),*)tefragment(ic,jc,kc,l)
        tefragment(ic,jc,kc,l)=tefragment(ic,jc,kc,l)-selfeng
        normal=.true.
        exit
       endif
       if(fline(2:27).eq.'Self energy of the charges')then
        read(fline(30:50),*)selfeng
       endif
!!     if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
!!      read(8,*)
!!      read(8,*)
!!      ic1=0
!!      do m=1,7
!!       if(group(l,m).ne.0)then
!!       read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!!       tforce(m,1)=tforce(m,1)+temp(1)
!!       tforce(m,2)=tforce(m,2)+temp(2)
!!       tforce(m,3)=tforce(m,3)+temp(3)
!!       ic1=ic1+1
!!       pforce(ic,jc,kc,l,ic1,1)=temp(1)
!!       pforce(ic,jc,kc,l,ic1,2)=temp(2)
!!       pforce(ic,jc,kc,l,ic1,3)=temp(3)
!!       endif
!!      enddo
!!      normal=.true.
!!      exit
!!     endif
       enddo
      close(8)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g09 '//cmd(i)(1:k-3)//'gjf')
       goto 2018
      endif
!      open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
!      do while(.true.)
!       read(6,'(a70)',iostat=ist)fline                                                 
!       if(ist/=0)exit
!       if(fline(42:55).eq.'Electric Field')then
!       read(6,*)
!       read(6,*)
!       do m=1,natom
!        if(resnum(m)==l)then
!         read(6,*)
!        endif
!       enddo
!       do ic1=-nxm,nxm
!        do jc1=-nym,nym
!         do kc1=-nzm,nzm
!          do m=1,natom
!        if((ic1.ne.ic).or.(jc1.ne.jc).or.(kc1.ne.kc).or.(resnum(m)/=l))then
!           read(6,'(28x,3f18.10)')(temp(k),k=1,3)
!           tfield(m,1)=tfield(m,1)+temp(1)
!           tfield(m,2)=tfield(m,2)+temp(2)
!           tfield(m,3)=tfield(m,3)+temp(3)
!           pfield(ic1,jc1,kc1,ic,jc,kc,l,m,1)=temp(1)
!           pfield(ic1,jc1,kc1,ic,jc,kc,l,m,2)=temp(2)
!           pfield(ic1,jc1,kc1,ic,jc,kc,l,m,3)=temp(3)
!        endif
!          enddo
!         enddo
!        enddo
!       enddo 
!       endif
!      enddo
!      close(6)

      endif        
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(pforce,pforce_global,nwat*3*10*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(pfield,pfield_global,nwat*27*natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      njob=0
      do i=1,ngroup-1
       do j=i+1,ngroup
        if(connect(0,0,0,i,j))then
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
      do ic=-nxm,nxm
       do jc=-nym,nym
        do kc=-nzm,nzm
         if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
          do i=1,ngroup
           do j=1,ngroup
            if(connect(ic,jc,kc,i,j))then
            njob=njob+1
            cmd(njob)=char(48+ic+nxm*2)//char(48+jc+nym*2)//&
                     char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
                char(48+(i-i/100*100)/10)//&
                char(48+(i-i/10*10))//&
                '-000cell'//char(48+j/100)//&
                char(48+(j-j/100*100)/10)//&
                char(48+(j-j/10*10))//'.log'
            endif
           enddo
          enddo
         endif
        enddo
       enddo
      enddo

      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).eq.'000')then

2019  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3,8x,I3)')j,l
      normal=.false.
      bgc=.true.
      bgc2=.true.
      do while(.true.)
      read(8,'(a80)',iostat=ist)fline
      if(ist/=0)exit 
      if(fline(28:32).eq.'EUMP2')then
!      indexx=index(fline,'=')+2                                            
!      read(fline(indexx:indexx+18),*)energy2(i,j)    
       read(fline(36:59),*)tenergy2(0,0,0,j,l)
       tenergy2(0,0,0,j,l)=tenergy2(0,0,0,j,l)-selfeng
      endif
      if(fline(2:27).eq.'Self energy of the charges')then  
       read(fline(30:50),*)selfeng                          
      endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(8,*)  
       read(8,*)  
       ic=0
       do m=1,20
       if(group(j,m).ne.0)then
       read(8,'(23x,3f15.9)')(temp(k),k=1,3)
       bgc(group(j,m))=.false.
       ic=ic+1
       tforce(group(j,m),1)=tforce(group(j,m),1)+temp(1)-&
                                         pforce_global(0,0,0,j,ic,1)
       tforce(group(j,m),2)=tforce(group(j,m),2)+temp(2)-&
                                         pforce_global(0,0,0,j,ic,2)
       tforce(group(j,m),3)=tforce(group(j,m),3)+temp(3)-&
                                         pforce_global(0,0,0,j,ic,3)
       else
        exit
       endif
       enddo
       ic=0
       do m=1,20
       if(group(l,m).ne.0)then
       ic=ic+1
       read(8,'(23x,3f15.9)')(temp(k),k=1,3)
       bgc2(group(l,m))=.false.
       tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)-&
                                         pforce_global(0,0,0,l,ic,1)
       tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)-&
                                         pforce_global(0,0,0,l,ic,2)
       tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)-&
                                         pforce_global(0,0,0,l,ic,3)
       else
        exit
       endif
       enddo
       normal=.true.
       exit
      endif
      enddo
      close(8)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g09 '//cmd(i)(1:k-3)//'gjf')
       goto 2019
      endif

1991     open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
         normal=.false.
         do while(.true.)
           read(6,'(a70)',iostat=ist)fline
           if(ist/=0)exit
           if(fline(42:55).eq.'Electric Field')then
          read(6,*) 
          read(6,*)
           do m=1,natom
            if((.not.bgc(m)).or.(.not.bgc2(m)))then
             read(6,*)
            endif
           enddo
           do m=1,natom
           if(bgc(m).and.bgc2(m))then
           read(6,'(28x,3f18.10)')(temp(k),k=1,3)
           tfield(0,0,0,m,1)=tfield(0,0,0,m,1)+temp(1)
           tfield(0,0,0,m,2)=tfield(0,0,0,m,2)+temp(2)
           tfield(0,0,0,m,3)=tfield(0,0,0,m,3)+temp(3)
           endif
           enddo
!          do ic=-nxm,nxm
!           do jc=-nym,nym
!            do kc=-nzm,nzm
!             if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!              do m=1,natom
!                  read(6,'(28x,3f18.10)')(temp(k),k=1,3)
!                  tfield(ic,jc,kc,m,1)=tfield(ic,jc,kc,m,1)+temp(1)-&
!                         pfield_global(j,ic,jc,kc,m,1)-&
!                         pfield_global(l,ic,jc,kc,m,1)
!                  tfield(ic,jc,kc,m,2)=tfield(ic,jc,kc,m,2)+temp(2)-&
!                         pfield_global(j,ic,jc,kc,m,2)-&
!                         pfield_global(l,ic,jc,kc,m,2)
!                  tfield(ic,jc,kc,m,3)=tfield(ic,jc,kc,m,3)+temp(3)-&
!                         pfield_global(j,ic,jc,kc,m,3)-&
!                         pfield_global(l,ic,jc,kc,m,3)
!              enddo
!             endif
!            enddo
!           enddo
!          enddo
           do m=1,natom
           if(bgc(m))then
          tfield(0,0,0,m,1)=tfield(0,0,0,m,1)-pfield_global(j,0,0,0,m,1)
          tfield(0,0,0,m,2)=tfield(0,0,0,m,2)-pfield_global(j,0,0,0,m,2)
          tfield(0,0,0,m,3)=tfield(0,0,0,m,3)-pfield_global(j,0,0,0,m,3)
           endif
           enddo
           do m=1,natom
           if(bgc2(m))then
          tfield(0,0,0,m,1)=tfield(0,0,0,m,1)-pfield_global(l,0,0,0,m,1)
          tfield(0,0,0,m,2)=tfield(0,0,0,m,2)-pfield_global(l,0,0,0,m,2)
          tfield(0,0,0,m,3)=tfield(0,0,0,m,3)-pfield_global(l,0,0,0,m,3)
           endif
           enddo
           normal=.true.
           exit
           endif
         enddo
         close(6)
         if(.not.normal)then
!         call sleep(999999)
          k=len(trim(adjustl(cmd(i))))
          call system('g09 field_'//cmd(i)(1:k-3)//'gjf')
          goto 1991
         endif

      else

2020  open(8,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(3I1,4x,I3,8x,I3)')ic,jc,kc,j,l
      ic=ic-nxm*2
      jc=jc-nym*2
      kc=kc-nzm*2
      normal=.false.
      bgc=.true.
      bgc2=.true.
          do while(.true.)
          read(8,'(a80)',iostat=ist)fline
          if(ist/=0)exit 
          if(fline(28:32).eq.'EUMP2')then
!          indexx=index(fline,'=')+2                                            
!          read(fline(indexx:indexx+18),*)energy2(i,j)    
           read(fline(36:59),*)tenergy2(ic,jc,kc,j,l)
           tenergy2(ic,jc,kc,j,l)=tenergy2(ic,jc,kc,j,l)-selfeng
          endif
          if(fline(2:27).eq.'Self energy of the charges')then  
           read(fline(30:50),*)selfeng                          
          endif
          if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
           read(8,*)  
           read(8,*)  
           ic1=0
           do m=1,20
           if(group(l,m).ne.0)then
           read(8,'(23x,3f15.9)')(temp(k),k=1,3)
           bgc(group(l,m))=.false.
           ic1=ic1+1
           tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)-&
                                            pforce_global(0,0,0,l,ic1,1)
           tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)-&
                                            pforce_global(0,0,0,l,ic1,2)
           tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)-&
                                            pforce_global(0,0,0,l,ic1,3)
           endif
           enddo
!          ic1=0
!          do m=1,natom
!          if(resnum(m)==j)then
!          ic1=ic1+1
!          read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!           da=da+(temp(1)-pforce_global(ic,jc,kc,j,ic1,1))*ic
!           db=db+(temp(2)-pforce_global(ic,jc,kc,j,ic1,2))*jc
!           dc=dc+(temp(3)-pforce_global(ic,jc,kc,j,ic1,3))*kc
!          endif
!          enddo
           normal=.true.
           exit
          endif
          enddo
         close(8)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g09 '//cmd(i)(1:k-3)//'gjf')
       goto 2020
      endif
!        open(6,file='field_'//trim(adjustl(cmd(i))),status='old')
!        do while(.true.)
!          read(6,'(a70)',iostat=ist)fline
!          if(ist/=0)exit
!          if(fline(42:55).eq.'Electric Field')then
!         read(6,*) 
!         read(6,*)
!          do m=1,natom
!           if((resnum(m)==j).or.(resnum(m)==l))then
!            read(6,*)
!           endif
!          enddo
!          do m=1,natom
!          if(resnum(m)/=l)then
!          read(6,'(28x,3f18.10)')(temp(k),k=1,3)
!          tfield(m,1)=tfield(m,1)+temp(1)
!          tfield(m,2)=tfield(m,2)+temp(2)
!          tfield(m,3)=tfield(m,3)+temp(3)
!          endif
!          enddo
!          do ic1=-nxm,nxm
!           do jc1=-nym,nym
!            do kc1=-nzm,nzm
!             if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
!              do m=1,natom
!               if((ic.ne.ic1).or.(jc.ne.jc1).or.(kc.ne.kc1).or.&
!                  (resnum(m).ne.j))then
!                  read(6,'(28x,3f18.10)')(temp(k),k=1,3)
!                  tfield(m,1)=tfield(m,1)+temp(1)
!                  tfield(m,2)=tfield(m,2)+temp(2)
!                  tfield(m,3)=tfield(m,3)+temp(3)
!               endif
!              enddo
!             endif
!            enddo
!           enddo
!          enddo
!          do m=1,natom
!          if(resnum(m)/=l)then
!          tfield(m,1)=tfield(m,1)-pfield_global(0,0,0,0,0,0,l,m,1)
!          tfield(m,2)=tfield(m,2)-pfield_global(0,0,0,0,0,0,l,m,2)
!          tfield(m,3)=tfield(m,3)-pfield_global(0,0,0,0,0,0,l,m,3)
!          endif
!          enddo
!          do ic1=-nxm,nxm
!           do jc1=-nym,nym
!            do kc1=-nzm,nzm
!             if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
!              do m=1,natom
!                  tfield(m,1)=tfield(m,1)-&
!                              pfield_global(ic1,jc1,kc1,0,0,0,l,m,1)
!                  tfield(m,2)=tfield(m,2)-&
!                              pfield_global(ic1,jc1,kc1,0,0,0,l,m,2)
!                  tfield(m,3)=tfield(m,3)-&
!                              pfield_global(ic1,jc1,kc1,0,0,0,l,m,3)
!              enddo
!             endif
!            enddo
!           enddo
!          enddo
!          do m=1,natom
!          tfield(m,1)=tfield(m,1)-pfield_global(0,0,0,ic,jc,kc,j,m,1)
!          tfield(m,2)=tfield(m,2)-pfield_global(0,0,0,ic,jc,kc,j,m,2)
!          tfield(m,3)=tfield(m,3)-pfield_global(0,0,0,ic,jc,kc,j,m,3)
!          enddo
!          do ic1=-nxm,nxm
!           do jc1=-nym,nym
!            do kc1=-nzm,nzm
!             if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
!              do m=1,natom
!               if((ic.ne.ic1).or.(jc.ne.jc1).or.(kc.ne.kc1).or.&
!                  (resnum(m).ne.j))then
!                  tfield(m,1)=tfield(m,1)-&
!                              pfield_global(ic1,jc1,kc1,ic,jc,kc,j,m,1)
!                  tfield(m,2)=tfield(m,2)-&
!                              pfield_global(ic1,jc1,kc1,ic,jc,kc,j,m,2)
!                  tfield(m,3)=tfield(m,3)-&
!                              pfield_global(ic1,jc1,kc1,ic,jc,kc,j,m,3)
!               endif
!              enddo
!             endif
!            enddo
!           enddo
!          enddo
!          endif
!        enddo
!        close(6)

      endif
      enddo
      
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call MPI_Allreduce(tforce,force,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tfield,field,3*3*3*natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tefragment,efragment,27*nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(tenergy2,energy2,nwat*nwat*27,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
      call MPI_Allreduce(txd,xd,nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tyd,yd,nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tzd,zd,nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tdip,dip,nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call MPI_Allreduce(da,delta_a,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
!     call MPI_Allreduce(db,delta_b,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
!     call MPI_Allreduce(dc,delta_c,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       

          charge=charge*18.2223
          fele=0.0d0
          fvdw=0.0d0
          tforce=0.0d0
          ele=0.0d0
          evdw=0.0d0
          e2body=0.0d0
          do i=mytid+1,ngroup-1,numprocs
           do j=i+1,ngroup
            if(.not.connect(0,0,0,i,j))then
             do k=1,20
             if(group(i,k).ne.0)then
              do l=1,20
              if(group(j,l).ne.0)then
            xa=coord(group(i,k),1)
            ya=coord(group(i,k),2)
            za=coord(group(i,k),3)
            xc=coord(group(j,l),1)
            yc=coord(group(j,l),2)
            zc=coord(group(j,l),3)
            dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
            dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
            ele=ele+charge(group(i,k))*charge(group(j,l))/dist
            fele(group(i,k),1)=fele(group(i,k),1)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            fele(group(i,k),2)=fele(group(i,k),2)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            fele(group(i,k),3)=fele(group(i,k),3)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
            fele(group(j,l),1)=fele(group(j,l),1)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            fele(group(j,l),2)=fele(group(j,l),2)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            fele(group(j,l),3)=fele(group(j,l),3)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
            iacj=ntype*(iac(group(i,k))-1)
            ic=ico(iacj+iac(group(j,l)))
            dinv6=1.0/(dist2*dist2*dist2)
            dinv12=dinv6**2
            dinv8=dinv6/dist2
            dinv14=dinv12/dist2
            evdw=evdw+cn1(ic)*dinv12-cn2(ic)*dinv6
            fvdw(group(i,k),1)=fvdw(group(i,k),1)+(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(xa-xc)
            fvdw(group(i,k),2)=fvdw(group(i,k),2)+(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(ya-yc)
            fvdw(group(i,k),3)=fvdw(group(i,k),3)+(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(za-zc)
            fvdw(group(j,l),1)=fvdw(group(j,l),1)-(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(xa-xc)
            fvdw(group(j,l),2)=fvdw(group(j,l),2)-(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(ya-yc)
            fvdw(group(j,l),3)=fvdw(group(j,l),3)-(12.0*cn1(ic)*dinv14-&
                              6.0*cn2(ic)*dinv8)*(za-zc)
              endif
              enddo
             endif
             enddo
            endif
           enddo
          enddo  

      da=0.0
      db=0.0
      dc=0.0
      do j=mytid+1,ngroup,numprocs
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!!         do ii=1,natom
!!          da=da+field(ic,jc,kc,ii,1)*charge(ii)/18.2223*ic
!!          db=db+field(ic,jc,kc,ii,2)*charge(ii)/18.2223*jc
!!          dc=dc+field(ic,jc,kc,ii,3)*charge(ii)/18.2223*kc
!!         enddo 
           do i=1,ngroup
             if(connect(ic,jc,kc,i,j))then
              e2body=e2body+(energy2(ic,jc,kc,i,j)-&
              efragment(ic,jc,kc,i)-efragment(0,0,0,j))*0.5d0
             else
              do ii=1,20
              if(group(i,ii).ne.0)then
               do jj=1,20
               if(group(j,jj).ne.0)then
                xa=coord(group(j,jj),1)
                ya=coord(group(j,jj),2)
                za=coord(group(j,jj),3)
                xc=cellcoord(ic,jc,kc,group(i,ii),1)
                yc=cellcoord(ic,jc,kc,group(i,ii),2)
                zc=cellcoord(ic,jc,kc,group(i,ii),3)
                dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
                dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
              ele=ele+0.5d0*charge(group(i,ii))*charge(group(j,jj))/dist
            iacj=ntype*(iac(group(i,ii))-1)
            ic1=ico(iacj+iac(group(j,jj)))
            dinv6=1.0/(dist2*dist2*dist2)
            dinv12=dinv6**2
            dinv8=dinv6/dist2
            dinv14=dinv12/dist2
            evdw=evdw+0.5d0*(cn1(ic1)*dinv12-cn2(ic1)*dinv6)
         fvdw(group(j,jj),1)=fvdw(group(j,jj),1)+(12.0*cn1(ic1)*dinv14-&
                              6.0*cn2(ic1)*dinv8)*(xa-xc)
         fvdw(group(j,jj),2)=fvdw(group(j,jj),2)+(12.0*cn1(ic1)*dinv14-&
                              6.0*cn2(ic1)*dinv8)*(ya-yc)
         fvdw(group(j,jj),3)=fvdw(group(j,jj),3)+(12.0*cn1(ic1)*dinv14-&
                              6.0*cn2(ic1)*dinv8)*(za-zc)
               endif
               enddo
              endif
              enddo
             endif
           enddo
          endif
         enddo
        enddo
       enddo
      enddo

      call MPI_Allreduce(fele,tforce,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
      call MPI_Allreduce(fvdw,fvdw2,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
      call MPI_Allreduce(ele,ele2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
      call MPI_Allreduce(evdw,evdw2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
      call MPI_Allreduce(e2body,energy2body,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    

          tforce=tforce*0.5291771d0/627.51
          fvdw2=fvdw2*0.5291771d0/627.51

      do i=1,natom
      force(i,1)=force(i,1)+field(0,0,0,i,1)*charge(i)/18.2223-&
                                           tforce(i,1)
      force(i,2)=force(i,2)+field(0,0,0,i,2)*charge(i)/18.2223-&
                                           tforce(i,2)
      force(i,3)=force(i,3)+field(0,0,0,i,3)*charge(i)/18.2223-&
                                           tforce(i,3)
      enddo
!==================end QM energy calculation=====================
!==============two body QM energy calculation====================
      if(master)then

       open(99,access='append',file='out.data') 
       do i=1,ngroup
!       write(99,2000)' 0 0 0 cell',i,'QM energy=',efragment(0,0,0,i)
        etotal=etotal+efragment(0,0,0,i)
       enddo
       
2000   format(a11,I4,1x,a12,f16.8)
2001   format(3I2,1x,a4,I4,1x,a12,f16.8)


      do i=1,ngroup-1
       do j=i+1,ngroup
        if(connect(0,0,0,i,j))then
        energy2body=energy2body+energy2(0,0,0,i,j)-efragment(0,0,0,i)-&
                    efragment(0,0,0,j)  
        endif
       enddo
      enddo

      energy2body=energy2body-ele2/627.51

!     write(99,*)'The self-energy of the cell (au): ',etotal+energy2body

      endif
!=========================Long range interaction===================================     
       ele2=0.0
       ele=0.0
       fele=0.0
       tforce=0.0

      do jj=mytid+1,natom,numprocs
       do ic=-nxl,nxl
        do jc=-nyl,nyl
         do kc=-nzl,nzl
          if((abs(ic).ge.2).or.(abs(jc).ge.2).or.(abs(kc).ge.2))then
              do ii=1,natom
                xa=coord(jj,1)
                ya=coord(jj,2)
                za=coord(jj,3)
                xc=coord(ii,1)+real(ic)*clata
                yc=coord(ii,2)+real(jc)*clatb
                zc=coord(ii,3)+real(kc)*clatc 
                dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
                ele2=ele2+0.5d0*charge(ii)*charge(jj)/dist
           fele(jj,1)=fele(jj,1)+charge(ii)*charge(jj)/(dist**3)*(xa-xc)
           fele(jj,2)=fele(jj,2)+charge(ii)*charge(jj)/(dist**3)*(ya-yc)
           fele(jj,3)=fele(jj,3)+charge(ii)*charge(jj)/(dist**3)*(za-zc)
!               da=da+(charge(ii)*charge(jj)/(dist**3)*(xc-xa))*ic*0.5291771d0/627.51
!               db=db+(charge(ii)*charge(jj)/(dist**3)*(yc-ya))*jc*0.5291771d0/627.51
!               dc=dc+(charge(ii)*charge(jj)/(dist**3)*(zc-za))*kc*0.5291771d0/627.51
              enddo
          endif
         enddo
        enddo
       enddo
      enddo

       call MPI_Allreduce(fele,tforce,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       call MPI_Allreduce(ele2,ele,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      tforce=tforce*0.5291771d0/627.51

      do i=1,natom
       force(i,1)=force(i,1)+tforce(i,1)
       force(i,2)=force(i,2)+tforce(i,2)
       force(i,3)=force(i,3)+tforce(i,3)
      enddo

!      delta_a=(delta_a+da)*0.5d0
!      delta_b=(delta_b+db)*0.5d0
!      delta_c=(delta_c+dc)*0.5d0

!      delta_a=delta_a-p*clatb*clatc*1.21374312/10000.0
!      delta_b=delta_b-p*clata*clatc*1.21374312/10000.0
!      delta_c=delta_c-p*clata*clatb*1.21374312/10000.0

2011  format(a11,1x,2I4,a20,2f16.10)
2012  format(3I2,1x,a4,1x,2I4,a20,2f16.10)
      
      etotal=etotal+energy2body+ele/627.51
!     enthalpy=etotal+p*clata*clatb*clatc*2.29364191/10000.0

      if(master)then

     open(001,file='force_all.dat',access='append')
      do i=1,natom
       write(001,'(I5,3x,3E16.8)')i,(force(i,k),k=1,3)
      enddo
      write(001,'(a)')'======='
     close(001)

     write(99,*)'Total energy of the cell (au):',etotal

       close(99)

      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

       end subroutine
	
