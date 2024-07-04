      subroutine energy_cal
      use variables
      use ewald
      implicit none
      include 'parallel.h'
      include 'mpif.h'

      real(kind=8)::enthalpy,xdipo,ydipo,zdipo,dipo,&
                    xdipo2,ydipo2,zdipo2,dipo2
                    

      energy2=0.d0
      tenergy2=0.d0
      efragment=0.d0
      tefragment=0.d0
      etotal=0.d0
      force=0.d0
      tforce=0.d0
      pforce=0.d0
      pforce_global=0.d0
      field=0.d0
      tfield=0.d0
      pfield=0.d0
      pfield_global=0.d0
      f2b=0.d0
      ff2b=0.d0
      da=0.d0
      db=0.d0
      dc=0.d0
      xd=0.d0
      txd=0.d0
      yd=0.d0
      tyd=0.d0
      zd=0.d0
      tzd=0.d0
      dip=0.d0
      tdip=0.d0
      xdipo=0.d0
      ydipo=0.d0
      zdipo=0.d0
      dipo=0.d0
      e3b=0.d0
      e3b2=0.d0

      n=natom
      xbox=clata
      ybox=clatb
      zbox=clatc
      x(1:n)=coord(1:n,1)
      y(1:n)=coord(1:n,2)
      z(1:n)=coord(1:n,3)

      charge=charge*18.2223
      njob=0

       do i=1,nwat
        njob=njob+1
        cmd(njob)='g16 000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo

        do ic=-nxm,nxm
         do jc=-nym,nym
          do kc=-nzm,nzm
           if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
            do i=1,nwat
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

       do i=1,nwat-1
        do j=i+1,nwat
         if(connect(0,0,0,i,j))then
          njob=njob+1
          cmd(njob)='g16 000cell'//char(48+i/100)//&
                             char(48+(i-i/100*100)/10)//&
          char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
          char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.gjf'
         endif
        enddo
       enddo

       do j=1,nwat
        do ic=-nxm,nxm
         do jc=-nym,nym
          do kc=-nzm,nzm
           if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
            do i=1,nwat
             if(connect(ic,jc,kc,i,j))then
              njob=njob+1
              cmd(njob)='g16 '//char(48+ic+nxm*2)//&
         char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
         char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))//'-000cell'//&
         char(48+j/100)//char(48+(j-j/100*100)/10)//&
                         char(48+(j-j/10*10))//'.gjf'
             endif
            enddo
           endif
          enddo
         enddo
        enddo
       enddo

       do j=1,nr
        njob=njob+1
        cmd(njob)='g16 '//char(48+rbx(j)+2)//char(48+rby(j)+2)//&
                  char(48+rbz(j)+2)//'cell'//char(48+rb0(j)/100)//&
                             char(48+(rb0(j)-rb0(j)/100*100)/10)//&
                             char(48+(rb0(j)-rb0(j)/10*10))//&
                 '-'//char(48+rcx(j)+2)//char(48+rcy(j)+2)//&
                  char(48+rcz(j)+2)//'cell'//char(48+rc0(j)/100)//&
                             char(48+(rc0(j)-rc0(j)/100*100)/10)//&
                             char(48+(rc0(j)-rc0(j)/10*10))//'.gjf'
       enddo

       do j=1,n3b
        if((tbx(j)==0).and.(tby(j)==0).and.(tbz(j)==0).and.&
           (tcx(j)==0).and.(tcy(j)==0).and.(tcz(j)==0))then
         njob=njob+1
         cmd(njob)='g16 '//'000cell'//char(48+ta0(j)/100)//&
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&
                               char(48+(ta0(j)-ta0(j)/10*10))//&
                   '-000cell'//char(48+tb0(j)/100)//&
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&
                               char(48+(tb0(j)-tb0(j)/10*10))//&
                   '-000cell'//char(48+tc0(j)/100)//&
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&
                               char(48+(tc0(j)-tc0(j)/10*10))//'.gjf'
        endif
        if((tbx(j)==0).and.(tby(j)==0).and.(tbz(j)==0))then
         if((tcx(j)/=0).or.(tcy(j)/=0).or.(tcz(j)/=0))then
          njob=njob+1
          cmd(njob)='g16 '//'000cell'//char(48+ta0(j)/100)//&
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&
                               char(48+(ta0(j)-ta0(j)/10*10))//&
                   '-000cell'//char(48+tb0(j)/100)//&
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&
                               char(48+(tb0(j)-tb0(j)/10*10))//&
                   '-'//char(48+tcx(j)+2)//char(48+tcy(j)+2)//&
                    char(48+tcz(j)+2)//'cell'//char(48+tc0(j)/100)//&
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&
                               char(48+(tc0(j)-tc0(j)/10*10))//'.gjf'
         endif
        endif
        if((tbx(j)/=0).or.(tby(j)/=0).or.(tbz(j)/=0))then
         if((tcx(j)/=0).or.(tcy(j)/=0).or.(tcz(j)/=0))then
          njob=njob+1
          cmd(njob)='g16 '//'000cell'//char(48+ta0(j)/100)//&
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&
                               char(48+(ta0(j)-ta0(j)/10*10))//&
                   '-'//char(48+tbx(j)+2)//char(48+tby(j)+2)//&
                    char(48+tbz(j)+2)//'cell'//char(48+tb0(j)/100)//&
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&
                               char(48+(tb0(j)-tb0(j)/10*10))//&
                   '-'//char(48+tcx(j)+2)//char(48+tcy(j)+2)//&
                    char(48+tcz(j)+2)//'cell'//char(48+tc0(j)/100)//&
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&
                               char(48+(tc0(j)-tc0(j)/10*10))//'.gjf'
         endif
        endif
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
     
       do i=mytid+1,njob,numprocs
         call system(cmd(i))
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

       njob=0

       do i=1,nwat
        njob=njob+1
        cmd(njob)='g16 field_000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.gjf'
       enddo
       do i=1,nwat-1
        do j=i+1,nwat
         if(connect(0,0,0,i,j))then
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

!      njob=0

!      do i=1,nwat
!       njob=njob+1
!       cmd(njob)='formchk 000cell'//char(48+i/100)//&
!                        char(48+(i-i/100*100)/10)//&
!                        char(48+(i-i/10*10))//'.chk'
!      enddo

!       do ic=-nxm,nxm
!        do jc=-nym,nym
!         do kc=-nzm,nzm
!          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!           do i=1,nwat
!            if(onebody(ic,jc,kc,i))then
!             njob=njob+1
!             cmd(njob)='formchk '//char(48+ic+nxm*2)//&
!       char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
!                        char(48+(i-i/100*100)/10)//&
!                        char(48+(i-i/10*10))//'.chk'
!            endif
!           enddo
!          endif
!         enddo
!        enddo
!       enddo

!       do i=1,nwat-1
!       do j=i+1,nwat
!        if(connect(0,0,0,i,j))then
!         njob=njob+1
!         cmd(njob)='formchk 000cell'//char(48+i/100)//&
!                            char(48+(i-i/100*100)/10)//&
!         char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
!         char(48+(j-j/100*100)/10)//char(48+(j-j/10*10))//'.chk'
!        endif
!       enddo
!      enddo

!      do j=1,nwat
!       do ic=-nxm,nxm
!        do jc=-nym,nym
!         do kc=-nzm,nzm
!          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!           do i=j,nwat
!            if(connect(ic,jc,kc,i,j))then
!             njob=njob+1
!             cmd(njob)='formchk '//char(48+ic+nxm*2)//&
!        char(48+jc+nym*2)//char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
!        char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))//'-000cell'//&
!        char(48+j/100)//char(48+(j-j/100*100)/10)//&
!                        char(48+(j-j/10*10))//'.chk'
!            endif
!           enddo
!          endif
!         enddo
!        enddo
!       enddo
!      enddo

!      call MPI_Barrier(MPI_COMM_WORLD,ierr)

!      do i=mytid+1,njob,numprocs
!       call system(cmd(i))
!      enddo

!      call MPI_Barrier(MPI_COMM_WORLD,ierr)
!====================== QM energy calculation=====================
      njob=0
      do i=1,nwat
       njob=njob+1
       cmd(njob)='000cell'//char(48+i/100)//&
                   char(48+(i-i/100*100)/10)//&
                   char(48+(i-i/10*10))//'.log'
      enddo
      do ic=-nxm,nxm
       do jc=-nym,nym
        do kc=-nzm,nzm
         if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
          do i=1,nwat
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
2017  open(108,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3)')l
      normal=.false.
      do while(.true.)
      read(108,'(a80)',iostat=ist)fline
      if(ist/=0)exit
      if(fline(28:35).eq.'E(CORR)=')then
      read(fline(37:53),*)tefragment(0,0,0,l)
      tefragment(0,0,0,l)=tefragment(0,0,0,l)-selfeng
      endif
      if(fline(2:27).eq.'Self energy of the charges')then
      read(fline(30:50),*)selfeng
      endif
      if(fline(2:14).eq.'Dipole moment')then
       read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')txd(0,0,0,l),&
           tyd(0,0,0,l),tzd(0,0,0,l),tdip(0,0,0,l)
       xdipo=xdipo+txd(0,0,0,l)
       ydipo=ydipo+tyd(0,0,0,l)
       zdipo=zdipo+tzd(0,0,0,l)
       dipo=dipo+tdip(0,0,0,l)
      endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
         read(108,*)
         read(108,*)
         do m=1,3
         read(108,'(23x,3f15.9)')(temp(k),k=1,3)
         bgc(group(l,m))=.false.
         tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)
         tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)
         tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)
         pforce(0,0,0,l,m,1)=temp(1)
         pforce(0,0,0,l,m,2)=temp(2)
         pforce(0,0,0,l,m,3)=temp(3)
         enddo
         normal=.true.
         exit
      endif
      enddo
      close(108)
!     open(9,file=cmd(i)(1:10)//'.fchk')
!      do while(.true.)
!      read(9,'(a80)')fline
!      if(fline(1:18).eq.'Cartesian Gradient')then
!       read(fline(54:61),'(I8)')m
!       read(9,'(5E16.8)')(temp(k),k=1,m)
!       temp=-temp
!       ic=0
!       do m=1,7
!       if(group(l,m).ne.0)then
!        bgc(group(l,m))=.false.
!        ic=ic+1
!        tforce(group(l,m),1)=tforce(group(l,m),1)+temp(ic*3-2)
!        tforce(group(l,m),2)=tforce(group(l,m),2)+temp(ic*3-1)
!        tforce(group(l,m),3)=tforce(group(l,m),3)+temp(ic*3)
!        pforce(0,0,0,l,ic,1)=temp(ic*3-2)
!        pforce(0,0,0,l,ic,2)=temp(ic*3-1)
!        pforce(0,0,0,l,ic,3)=temp(ic*3)        
!       endif
!       enddo
!       normal=.true.
!      endif
!      if(fline(1:18).eq.'Dipole Moment')then                                           
!       read(9,'(3E16.8)')txd(l),tyd(l),tzd(l)          
!       exit                                      
!      endif
!      enddo
!     close(9)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g16 '//cmd(i)(1:k-3)//'gjf')
       goto 2017
      endif
      call system("rm "//cmd(i)(1:10)//".*")

1990   open(106,file='field_'//trim(adjustl(cmd(i))),status='old')
       normal=.false.
       do while(.true.)
        read(106,'(a70)',iostat=ist)fline                                                 
        if(ist/=0)exit
        if(fline(42:55).eq.'Electric Field')then
        read(106,*)
        read(106,*)
        do m=1,3
          read(106,*)
        enddo
        do m=1,natom
         if(bgc(m))then
        read(106,'(28x,3f18.10)')(temp(k),k=1,3)
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
       close(106)
       if(.not.normal)then
!       call sleep(999999)
        k=len(trim(adjustl(cmd(i))))
        call system('g16 field_'//cmd(i)(1:k-3)//'gjf')
        goto 1990
       endif
       call system("rm field_"//cmd(i)(1:10)//".*")

      endif
      enddo

      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).ne.'000')then
          
2018  open(108,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(3I1,4x,I3)')ic,jc,kc,l
      ic=ic-nxm*2
      jc=jc-nym*2
      kc=kc-nzm*2
      normal=.false.
       do while(.true.)
       read(108,'(a80)',iostat=ist)fline
       if(ist/=0)exit
       if(fline(28:35).eq.'E(CORR)=')then
        read(fline(37:53),*)tefragment(ic,jc,kc,l)
        tefragment(ic,jc,kc,l)=tefragment(ic,jc,kc,l)-selfeng
       endif
       if(fline(2:27).eq.'Self energy of the charges')then
        read(fline(30:50),*)selfeng
       endif
       if(fline(2:14).eq.'Dipole moment')then
      read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')txd(ic,jc,kc,l),&
           tyd(ic,jc,kc,l),tzd(ic,jc,kc,l),tdip(ic,jc,kc,l)
       endif
       if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
        read(108,*)
        read(108,*)
        do m=1,3
         read(108,'(23x,3f15.9)')(temp(k),k=1,3)
         pforce(ic,jc,kc,l,m,1)=temp(1)
         pforce(ic,jc,kc,l,m,2)=temp(2)
         pforce(ic,jc,kc,l,m,3)=temp(3)
        enddo
        normal=.true.
        exit
       endif
       enddo
      close(108)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g16 '//cmd(i)(1:k-3)//'gjf')
       goto 2018
      endif
      call system("rm "//cmd(i)(1:10)//".*")
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
      call MPI_Allreduce(txd,xd,nwat*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(tyd,yd,nwat*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(tzd,zd,nwat*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(tdip,dip,nwat*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(pforce,pforce_global,nwat*3*7*27,& 
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 
      call MPI_Allreduce(pfield,pfield_global,nwat*27*natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      njob=0
      do i=1,nwat-1
       do j=i+1,nwat
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
      do j=1,nwat
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do i=1,nwat
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
          endif
         enddo
        enddo
       enddo
      enddo

      do i=mytid+1,njob,numprocs
      if(cmd(i)(1:3).eq.'000')then

2019  open(108,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(7x,I3,8x,I3)')j,l
      normal=.false.
      bgc=.true.
      bgc2=.true.
      do while(.true.)
      read(108,'(a80)',iostat=ist)fline
      if(ist/=0)exit 
      if(fline(28:35).eq.'E(CORR)=')then
       read(fline(37:53),*)tenergy2(0,0,0,j,l)
       tenergy2(0,0,0,j,l)=tenergy2(0,0,0,j,l)-selfeng
      endif
      if(fline(2:27).eq.'Self energy of the charges')then  
       read(fline(30:50),*)selfeng                          
      endif
      if(fline(2:14).eq.'Dipole moment')then
      read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')tx2d(0,0,0,j,l),&
           ty2d(0,0,0,j,l),tz2d(0,0,0,j,l),t2dip(0,0,0,j,l)
        xdipo=xdipo+tx2d(0,0,0,j,l)-xd(0,0,0,j)-xd(0,0,0,l)
        ydipo=ydipo+ty2d(0,0,0,j,l)-yd(0,0,0,j)-yd(0,0,0,l)
        zdipo=zdipo+tz2d(0,0,0,j,l)-zd(0,0,0,j)-zd(0,0,0,l)
        dipo=dipo+t2dip(0,0,0,j,l)-dip(0,0,0,j)-dip(0,0,0,l)
      endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(108,*)  
       read(108,*)  
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       f2b(0,0,0,j,l,m,1:3)=temp(1:3)
       bgc(group(j,m))=.false.
       tforce(group(j,m),1)=tforce(group(j,m),1)+temp(1)-&
                                         pforce_global(0,0,0,j,m,1)
       tforce(group(j,m),2)=tforce(group(j,m),2)+temp(2)-&
                                         pforce_global(0,0,0,j,m,2)
       tforce(group(j,m),3)=tforce(group(j,m),3)+temp(3)-&
                                         pforce_global(0,0,0,j,m,3)
       enddo
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       f2b(0,0,0,j,l,3+m,1:3)=temp(1:3)
       bgc2(group(l,m))=.false.
       tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)-&
                                         pforce_global(0,0,0,l,m,1)
       tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)-&
                                         pforce_global(0,0,0,l,m,2)
       tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)-&
                                         pforce_global(0,0,0,l,m,3)
       enddo
       normal=.true.
       exit
      endif
      enddo
      close(108)
!     open(9,file=cmd(i)(1:21)//'.fchk')
!      do while(.true.)
!      read(9,'(a80)')fline
!      if(fline(1:18).eq.'Cartesian Gradient')then                                             
!       read(fline(54:61),'(I8)')m
!       read(9,'(5E16.8)')(temp(k),k=1,m)                                                      
!       temp=-temp                                                                             
!       ic=0
!       ii=0
!       jj=0
!       do m=1,7
!       if(group(j,m).ne.0)then
!        bgc(group(j,m))=.false.
!        ic=ic+1
!        ii=ii+1
!        tforce(group(j,m),1)=tforce(group(j,m),1)+temp(ic*3-2)-&
!                             pforce_global(0,0,0,j,ii,1)
!        tforce(group(j,m),2)=tforce(group(j,m),2)+temp(ic*3-1)-&
!                             pforce_global(0,0,0,j,ii,2)
!        tforce(group(j,m),3)=tforce(group(j,m),3)+temp(ic*3)-&
!                             pforce_global(0,0,0,j,ii,3)
!       else
!        exit
!       endif
!       enddo
!       do m=1,7
!       if(group(l,m).ne.0)then
!        bgc2(group(l,m))=.false.
!        ic=ic+1
!        jj=jj+1
!        tforce(group(l,m),1)=tforce(group(l,m),1)+temp(ic*3-2)-&
!                             pforce_global(0,0,0,l,jj,1)
!        tforce(group(l,m),2)=tforce(group(l,m),2)+temp(ic*3-1)-&
!                             pforce_global(0,0,0,l,jj,2)
!        tforce(group(l,m),3)=tforce(group(l,m),3)+temp(ic*3)-&
!                             pforce_global(0,0,0,l,jj,3)
!       else
!        exit
!       endif
!       enddo
!       normal=.true.
!      endif
!      if(fline(1:18).eq.'Dipole Moment')then                                           
!       read(9,'(3E16.8)')tx2d(j,l),ty2d(j,l),tz2d(j,l)
!       exit
!      endif
!      enddo
!     close(9)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g16 '//cmd(i)(1:k-3)//'gjf')
       goto 2019
      endif
      call system("rm "//cmd(i)(1:21)//".*")

1991     open(106,file='field_'//trim(adjustl(cmd(i))),status='old')
         normal=.false.
         do while(.true.)
           read(106,'(a70)',iostat=ist)fline
           if(ist/=0)exit
           if(fline(42:55).eq.'Electric Field')then
          read(106,*) 
          read(106,*)
           do m=1,3
             read(106,*)
           enddo
           do m=1,3
             read(106,*)
           enddo
           do m=1,natom
           if(bgc(m).and.bgc2(m))then
           read(106,'(28x,3f18.10)')(temp(k),k=1,3)
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
         close(106)
         if(.not.normal)then
!         call sleep(999999)
          k=len(trim(adjustl(cmd(i))))
          call system('g16 field_'//cmd(i)(1:k-3)//'gjf')
          goto 1991
         endif
         call system("rm field_"//cmd(i)(1:21)//".*")

      else

2020  open(108,file=trim(adjustl(cmd(i))),status='OLD')
      read(cmd(i),'(3I1,4x,I3,8x,I3)')ic,jc,kc,j,l
      ic=ic-nxm*2
      jc=jc-nym*2
      kc=kc-nzm*2
      normal=.false.
      bgc=.true.
      bgc2=.true.
          do while(.true.)
          read(108,'(a80)',iostat=ist)fline
          if(ist/=0)exit 
          if(fline(28:35).eq.'E(CORR)=')then
           read(fline(37:53),*)tenergy2(ic,jc,kc,j,l)
           tenergy2(ic,jc,kc,j,l)=tenergy2(ic,jc,kc,j,l)-selfeng
          endif
          if(fline(2:27).eq.'Self energy of the charges')then  
           read(fline(30:50),*)selfeng                          
          endif
      if(fline(2:14).eq.'Dipole moment')then
       read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')&
           tx2d(ic,jc,kc,j,l),&
           ty2d(ic,jc,kc,j,l),tz2d(ic,jc,kc,j,l),t2dip(ic,jc,kc,j,l)
        xdipo=xdipo+tx2d(ic,jc,kc,j,l)-xd(ic,jc,kc,j)-xd(0,0,0,l)
        ydipo=ydipo+ty2d(ic,jc,kc,j,l)-yd(ic,jc,kc,j)-yd(0,0,0,l)
        zdipo=zdipo+tz2d(ic,jc,kc,j,l)-zd(ic,jc,kc,j)-zd(0,0,0,l)
        dipo=dipo+t2dip(ic,jc,kc,j,l)-dip(ic,jc,kc,j)-dip(0,0,0,l)
      endif
          if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
           read(108,*)  
           read(108,*)  
           do m=1,3
           read(108,'(23x,3f15.9)')(temp(k),k=1,3)
           f2b(ic,jc,kc,j,l,m,1:3)=temp(1:3)
           bgc(group(l,m))=.false.
           tforce(group(l,m),1)=tforce(group(l,m),1)+temp(1)-&
                                            pforce_global(0,0,0,l,m,1)
           tforce(group(l,m),2)=tforce(group(l,m),2)+temp(2)-&
                                            pforce_global(0,0,0,l,m,2)
           tforce(group(l,m),3)=tforce(group(l,m),3)+temp(3)-&
                                            pforce_global(0,0,0,l,m,3)
           enddo
           do m=1,3
           read(108,'(23x,3f15.9)')(temp(k),k=1,3)
           f2b(ic,jc,kc,j,l,3+m,1:3)=temp(1:3)
           bgc2(group(j,m))=.false.
!          tforce(group(j,m),1)=tforce(group(j,m),1)+temp(1)-&
!                                        pforce_global(ic,jc,kc,j,m,1)
!          tforce(group(j,m),2)=tforce(group(j,m),2)+temp(2)-&
!                                        pforce_global(ic,jc,kc,j,m,2)
!          tforce(group(j,m),3)=tforce(group(j,m),3)+temp(3)-&
!                                        pforce_global(ic,jc,kc,j,m,3)
           enddo
 
!          do mi=group(j,1),group(j,3)
!           do mj=group(l,1),group(l,3)
!            xa=cellcoord(ic,jc,kc,mi,1)
!            ya=cellcoord(ic,jc,kc,mi,2)
!            za=cellcoord(ic,jc,kc,mi,3)
!            xc=coord(mj,1)
!            yc=coord(mj,2)
!            zc=coord(mj,3)
!            dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
!            dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
!            e3b=e3b+charge(mi)*charge(mj)/dist/627.51
!            tforce(mi,1)=tforce(mi,1)+&
!                charge(mi)*charge(mj)/(dist**3)*(xa-xc)*0.5291771d0/627.51
!            tforce(mi,2)=tforce(mi,2)+&
!                charge(mi)*charge(mj)/(dist**3)*(ya-yc)*0.5291771d0/627.51
!            tforce(mi,3)=tforce(mj,3)+&
!                charge(mi)*charge(mj)/(dist**3)*(za-zc)*0.5291771d0/627.51
!           enddo
!          enddo

!          do m=1,3
!          read(8,'(23x,3f15.9)')(temp(k),k=1,3)
!           da=da+(temp(1)-pforce_global(ic,jc,kc,j,m,1))*ic
!           db=db+(temp(2)-pforce_global(ic,jc,kc,j,m,2))*jc
!           dc=dc+(temp(3)-pforce_global(ic,jc,kc,j,m,3))*kc
!          enddo
           normal=.true.
           exit
          endif
          enddo
         close(108)
!     open(9,file=cmd(i)(1:21)//'.fchk')
!      do while(.true.)
!      read(9,'(a80)')fline
!      if(fline(1:18).eq.'Cartesian Gradient')then
!       read(fline(54:61),'(I8)')m
!       read(9,'(5E16.8)')(temp(k),k=1,m)
!       temp=-temp
!       ic1=0
!       do m=1,7
!       if(group(l,m).ne.0)then
!        bgc(group(l,m))=.false.
!        ic1=ic1+1
!        tforce(group(l,m),1)=tforce(group(l,m),1)+temp(ic1*3-2)-&
!                             pforce_global(0,0,0,l,ic1,1)
!        tforce(group(l,m),2)=tforce(group(l,m),2)+temp(ic1*3-1)-&
!                             pforce_global(0,0,0,l,ic1,2)
!        tforce(group(l,m),3)=tforce(group(l,m),3)+temp(ic1*3)-&
!                             pforce_global(0,0,0,l,ic1,3)
!       endif
!       enddo
!       normal=.true.
!       exit
!      endif
!      enddo
!     close(9)
      if(.not.normal)then
!      call sleep(999999)
       k=len(trim(adjustl(cmd(i))))
       call system('g16 '//cmd(i)(1:k-3)//'gjf')
       goto 2020
      endif
      call system("rm "//cmd(i)(1:21)//".*")
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

!     call MPI_Allreduce(tforce,force,natom*3,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tfield,field,3*3*3*natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tefragment,efragment,27*nwat,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(tenergy2,energy2,nwat*nwat*27,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(xdipo,xdipo2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(ydipo,ydipo2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(zdipo,zdipo2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(dipo,dipo2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(e3b,e3b2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(f2b,f2_global,27*128*128*6*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call MPI_Allreduce(da,delta_a,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
!     call MPI_Allreduce(db,delta_b,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)       
!     call MPI_Allreduce(dc,delta_c,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr) 

      do j=1,n3b
        if((tbx(j)==0).and.(tby(j)==0).and.(tbz(j)==0).and.&                         
           (tcx(j)==0).and.(tcy(j)==0).and.(tcz(j)==0))then
         cmd(j)='000cell'//char(48+ta0(j)/100)//&
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&
                               char(48+(ta0(j)-ta0(j)/10*10))//&
                   '-000cell'//char(48+tb0(j)/100)//&
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&                
                               char(48+(tb0(j)-tb0(j)/10*10))//&                     
                   '-000cell'//char(48+tc0(j)/100)//&                                
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&                
                               char(48+(tc0(j)-tc0(j)/10*10))//'.log'                
        endif
        if((tbx(j)==0).and.(tby(j)==0).and.(tbz(j)==0))then                          
         if((tcx(j)/=0).or.(tcy(j)/=0).or.(tcz(j)/=0))then
          cmd(j)='000cell'//char(48+ta0(j)/100)//&                        
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&                
                               char(48+(ta0(j)-ta0(j)/10*10))//&                     
                   '-000cell'//char(48+tb0(j)/100)//&                                
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&                
                               char(48+(tb0(j)-tb0(j)/10*10))//&                     
                   '-'//char(48+tcx(j)+2)//char(48+tcy(j)+2)//&                      
                    char(48+tcz(j)+2)//'cell'//char(48+tc0(j)/100)//&                
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&                
                               char(48+(tc0(j)-tc0(j)/10*10))//'.log'                
         endif
        endif
        if((tbx(j)/=0).or.(tby(j)/=0).or.(tbz(j)/=0))then
         if((tcx(j)/=0).or.(tcy(j)/=0).or.(tcz(j)/=0))then
          cmd(j)='000cell'//char(48+ta0(j)/100)//&
                               char(48+(ta0(j)-ta0(j)/100*100)/10)//&
                               char(48+(ta0(j)-ta0(j)/10*10))//&
                   '-'//char(48+tbx(j)+2)//char(48+tby(j)+2)//&
                    char(48+tbz(j)+2)//'cell'//char(48+tb0(j)/100)//&
                               char(48+(tb0(j)-tb0(j)/100*100)/10)//&
                               char(48+(tb0(j)-tb0(j)/10*10))//&
                   '-'//char(48+tcx(j)+2)//char(48+tcy(j)+2)//&
                    char(48+tcz(j)+2)//'cell'//char(48+tc0(j)/100)//&
                               char(48+(tc0(j)-tc0(j)/100*100)/10)//&
                               char(48+(tc0(j)-tc0(j)/10*10))//'.log'
         endif
        endif
      enddo

      do j=mytid+1,nr,numprocs
       open(108,file=char(48+rbx(j)+2)//char(48+rby(j)+2)//&
                    char(48+rbz(j)+2)//'cell'//char(48+rb0(j)/100)//&
                               char(48+(rb0(j)-rb0(j)/100*100)/10)//&
                               char(48+(rb0(j)-rb0(j)/10*10))//&
                   '-'//char(48+rcx(j)+2)//char(48+rcy(j)+2)//&
                    char(48+rcz(j)+2)//'cell'//char(48+rc0(j)/100)//&
                               char(48+(rc0(j)-rc0(j)/100*100)/10)//&
                               char(48+(rc0(j)-rc0(j)/10*10))//'.log')
       do k=1,27
        if((posx(k)==rbx(j)).and.(posy(k)==rby(j)).and.&
           (posz(k)==rbz(j)))then
         m1=k
        endif
        if((posx(k)==rcx(j)).and.(posy(k)==rcy(j)).and.&
           (posz(k)==rcz(j)))then
         m2=k
        endif
       enddo
       do while(.true.)
        read(108,'(a80)',iostat=ist)fline
        if(ist/=0)exit
        if(fline(28:35).eq.'E(CORR)=')then
         read(fline(37:53),*)temp(30)
         ee2b(m1,m2,rb0(j),rc0(j))=temp(30)-selfeng
        endif
        if(fline(2:27).eq.'Self energy of the charges')then
         read(fline(30:50),*)selfeng
        endif
        if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(108,*)
       read(108,*)
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       ff2b(m1,m2,rb0(j),rc0(j),m,1:3)=temp(1:3)
       enddo
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       ff2b(m1,m2,rb0(j),rc0(j),3+m,1:3)=temp(1:3)
       enddo
       exit
        endif
       enddo
       close(108)
       call system("rm "//char(48+rbx(j)+2)//char(48+rby(j)+2)//&
                          char(48+rbz(j)+2)//'cell'//char(48+rb0(j)/100)//&
                          char(48+(rb0(j)-rb0(j)/100*100)/10)//&
                          char(48+(rb0(j)-rb0(j)/10*10))//&
                     '-'//char(48+rcx(j)+2)//char(48+rcy(j)+2)//&
                          char(48+rcz(j)+2)//'cell'//char(48+rc0(j)/100)//&
                          char(48+(rc0(j)-rc0(j)/100*100)/10)//&
                          char(48+(rc0(j)-rc0(j)/10*10))//'.*')
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(ff2b,ff2_global,27*27*128*128*6*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(ee2b,ee2b_global,27*27*128*128,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      e3b_qm=0.d0
      do i=mytid+1,n3b,numprocs
       open(108,file=trim(adjustl(cmd(i))),status='OLD')

       if((tbx(i)==0).and.(tby(i)==0).and.(tbz(i)==0).and.&                         
           (tcx(i)==0).and.(tcy(i)==0).and.(tcz(i)==0))then
        do while(.true.)
        read(108,'(a80)',iostat=ist)fline
        if(ist/=0)exit 
        if(fline(28:35).eq.'E(CORR)=')then
         read(fline(37:53),*)temp(30)
         temp(30)=temp(30)-selfeng
        endif
        if(fline(2:27).eq.'Self energy of the charges')then  
         read(fline(30:50),*)selfeng                          
        endif
!       if(fline(2:14).eq.'Dipole moment')then
!     read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')tx2d(0,0,0,j,l),&
!          ty2d(0,0,0,j,l),tz2d(0,0,0,j,l),t2dip(0,0,0,j,l)
!        xdipo=xdipo+tx2d(0,0,0,j,l)-xd(0,0,0,j)-xd(0,0,0,l)
!        ydipo=ydipo+ty2d(0,0,0,j,l)-yd(0,0,0,j)-yd(0,0,0,l)
!        zdipo=zdipo+tz2d(0,0,0,j,l)-zd(0,0,0,j)-zd(0,0,0,l)
!        dipo=dipo+t2dip(0,0,0,j,l)-dip(0,0,0,j)-dip(0,0,0,l)
!       endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(108,*)  
       read(108,*)  
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(ta0(i),m),1)=tforce(group(ta0(i),m),1)+temp(1)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,1)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),m,1)+&
                                 pforce_global(0,0,0,ta0(i),m,1)
       tforce(group(ta0(i),m),2)=tforce(group(ta0(i),m),2)+temp(2)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,2)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),m,2)+&
                                 pforce_global(0,0,0,ta0(i),m,2)
       tforce(group(ta0(i),m),3)=tforce(group(ta0(i),m),3)+temp(3)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,3)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),m,3)+&
                                 pforce_global(0,0,0,ta0(i),m,3)
       enddo
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(tb0(i),m),1)=tforce(group(tb0(i),m),1)+temp(1)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,1)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),m,1)+&
                                 pforce_global(0,0,0,tb0(i),m,1)
       tforce(group(tb0(i),m),2)=tforce(group(tb0(i),m),2)+temp(2)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,2)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),m,2)+&
                                 pforce_global(0,0,0,tb0(i),m,2)
       tforce(group(tb0(i),m),3)=tforce(group(tb0(i),m),3)+temp(3)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,3)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),m,3)+&
                                 pforce_global(0,0,0,tb0(i),m,3)
       enddo
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(tc0(i),m),1)=tforce(group(tc0(i),m),1)+temp(1)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),3+m,1)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),3+m,1)+&
                                 pforce_global(0,0,0,tc0(i),m,1)
       tforce(group(tc0(i),m),2)=tforce(group(tc0(i),m),2)+temp(2)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),3+m,2)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),3+m,2)+&
                                 pforce_global(0,0,0,tc0(i),m,2)
       tforce(group(tc0(i),m),3)=tforce(group(tc0(i),m),3)+temp(3)-&
                                 f2_global(0,0,0,ta0(i),tc0(i),3+m,3)-&
                                 f2_global(0,0,0,tb0(i),tc0(i),3+m,3)+&
                                 pforce_global(0,0,0,tc0(i),m,3)
       enddo
       exit
      endif
        enddo
         e3b_qm=e3b_qm+temp(30)-energy2(0,0,0,ta0(i),tb0(i))-&
                               energy2(0,0,0,ta0(i),tc0(i))-&
                               energy2(0,0,0,tb0(i),tc0(i))+&
              efragment(0,0,0,ta0(i))+efragment(0,0,0,tb0(i))+&
              efragment(0,0,0,tc0(i))
        call system("rm "//cmd(i)(1:32)//".*") 
       endif

       if((tbx(i)==0).and.(tby(i)==0).and.(tbz(i)==0))then                          
         if((tcx(i)/=0).or.(tcy(i)/=0).or.(tcz(i)/=0))then
        do while(.true.)
        read(108,'(a80)',iostat=ist)fline
        if(ist/=0)exit 
        if(fline(28:35).eq.'E(CORR)=')then
         read(fline(37:53),*)temp(30)
         temp(30)=temp(30)-selfeng
        endif
        if(fline(2:27).eq.'Self energy of the charges')then  
         read(fline(30:50),*)selfeng                          
        endif
!       if(fline(2:14).eq.'Dipole moment')then
!     read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')tx2d(0,0,0,j,l),&
!          ty2d(0,0,0,j,l),tz2d(0,0,0,j,l),t2dip(0,0,0,j,l)
!        xdipo=xdipo+tx2d(0,0,0,j,l)-xd(0,0,0,j)-xd(0,0,0,l)
!        ydipo=ydipo+ty2d(0,0,0,j,l)-yd(0,0,0,j)-yd(0,0,0,l)
!        zdipo=zdipo+tz2d(0,0,0,j,l)-zd(0,0,0,j)-zd(0,0,0,l)
!        dipo=dipo+t2dip(0,0,0,j,l)-dip(0,0,0,j)-dip(0,0,0,l)
!       endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(108,*)  
       read(108,*)  
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(ta0(i),m),1)=tforce(group(ta0(i),m),1)+temp(1)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,1)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,1)+&
                                 pforce_global(0,0,0,ta0(i),m,1)
       tforce(group(ta0(i),m),2)=tforce(group(ta0(i),m),2)+temp(2)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,2)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,2)+&
                                 pforce_global(0,0,0,ta0(i),m,2)
       tforce(group(ta0(i),m),3)=tforce(group(ta0(i),m),3)+temp(3)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),m,3)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,3)+&
                                 pforce_global(0,0,0,ta0(i),m,3)
       enddo
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(tb0(i),m),1)=tforce(group(tb0(i),m),1)+temp(1)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,1)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),m,1)+&
                                 pforce_global(0,0,0,tb0(i),m,1)
       tforce(group(tb0(i),m),2)=tforce(group(tb0(i),m),2)+temp(2)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,2)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),m,2)+&
                                 pforce_global(0,0,0,tb0(i),m,2)
       tforce(group(tb0(i),m),3)=tforce(group(tb0(i),m),3)+temp(3)-&
                                 f2_global(0,0,0,ta0(i),tb0(i),3+m,3)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),m,3)+&
                                 pforce_global(0,0,0,tb0(i),m,3)
       enddo
!      do m=1,3
!      read(108,'(23x,3f15.9)')(temp(k),k=1,3)
!      tforce(group(tc0(i),m),1)=tforce(group(tc0(i),m),1)+temp(1)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,1)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),3+m,1)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,1)
!      tforce(group(tc0(i),m),2)=tforce(group(tc0(i),m),2)+temp(2)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,2)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),3+m,2)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,2)
!      tforce(group(tc0(i),m),3)=tforce(group(tc0(i),m),3)+temp(3)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,3)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i),3+m,3)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,3)
!      enddo
       exit
      endif
        enddo
         e3b_qm=e3b_qm+(temp(30)-energy2(0,0,0,ta0(i),tb0(i))-&
                energy2(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i))-&
                energy2(tcx(i),tcy(i),tcz(i),tc0(i),tb0(i))+&
              efragment(0,0,0,ta0(i))+efragment(0,0,0,tb0(i))+&
              efragment(tcx(i),tcy(i),tcz(i),tc0(i)))*2.d0/3.d0
         endif
        call system("rm "//cmd(i)(1:32)//".*")
       endif

       if((tbx(i)/=0).or.(tby(i)/=0).or.(tbz(i)/=0))then                          
         if((tcx(i)/=0).or.(tcy(i)/=0).or.(tcz(i)/=0))then
        do while(.true.)
        read(108,'(a80)',iostat=ist)fline
        if(ist/=0)exit 
        if(fline(28:35).eq.'E(CORR)=')then
         read(fline(37:53),*)temp(30)
         temp(30)=temp(30)-selfeng
        endif
        if(fline(2:27).eq.'Self energy of the charges')then  
         read(fline(30:50),*)selfeng                          
        endif
!       if(fline(2:14).eq.'Dipole moment')then
!     read(108,'(7x,f19.4,7x,f19.4,7x,f19.4,7x,f19.4)')tx2d(0,0,0,j,l),&
!          ty2d(0,0,0,j,l),tz2d(0,0,0,j,l),t2dip(0,0,0,j,l)
!        xdipo=xdipo+tx2d(0,0,0,j,l)-xd(0,0,0,j)-xd(0,0,0,l)
!        ydipo=ydipo+ty2d(0,0,0,j,l)-yd(0,0,0,j)-yd(0,0,0,l)
!        zdipo=zdipo+tz2d(0,0,0,j,l)-zd(0,0,0,j)-zd(0,0,0,l)
!        dipo=dipo+t2dip(0,0,0,j,l)-dip(0,0,0,j)-dip(0,0,0,l)
!       endif
      if(fline(38:59).eq.'Forces (Hartrees/Bohr)')then
       read(108,*)  
       read(108,*)  
       do m=1,3
       read(108,'(23x,3f15.9)')(temp(k),k=1,3)
       tforce(group(ta0(i),m),1)=tforce(group(ta0(i),m),1)+temp(1)-&
              f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),m,1)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,1)+&
                                 pforce_global(0,0,0,ta0(i),m,1)
       tforce(group(ta0(i),m),2)=tforce(group(ta0(i),m),2)+temp(2)-&
              f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),m,2)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,2)+&
                                 pforce_global(0,0,0,ta0(i),m,2)
       tforce(group(ta0(i),m),3)=tforce(group(ta0(i),m),3)+temp(3)-&
              f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),m,3)-&
              f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),m,3)+&
                                 pforce_global(0,0,0,ta0(i),m,3)
       enddo
!      do m=1,3
!      read(108,'(23x,3f15.9)')(temp(k),k=1,3)
!      tforce(group(tb0(i),m),1)=tforce(group(tb0(i),m),1)+temp(1)-&
!             f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),3+m,1)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),m,1)+&
!             pforce_global(tbx(i),tby(i),tbz(i),tb0(i),m,1)
!      tforce(group(tb0(i),m),2)=tforce(group(tb0(i),m),2)+temp(2)-&
!             f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),3+m,2)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),m,2)+&
!             pforce_global(tbx(i),tby(i),tbz(i),tb0(i),m,2)
!      tforce(group(tb0(i),m),3)=tforce(group(tb0(i),m),3)+temp(3)-&
!             f2_global(tbx(i),tby(i),tbz(i),tb0(i),ta0(i),3+m,3)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),m,3)+&
!             pforce_global(tbx(i),tby(i),tbz(i),tb0(i),m,3)
!      enddo
!      do m=1,3
!      read(108,'(23x,3f15.9)')(temp(k),k=1,3)
!      tforce(group(tc0(i),m),1)=tforce(group(tc0(i),m),1)+temp(1)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,1)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),3+m,1)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,1)
!      tforce(group(tc0(i),m),2)=tforce(group(tc0(i),m),2)+temp(2)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,2)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),3+m,2)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,2)
!      tforce(group(tc0(i),m),3)=tforce(group(tc0(i),m),3)+temp(3)-&
!             f2_global(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i),3+m,3)-&
!             ff2_global(m1,m2,tb0(i),tc0(i),3+m,3)+&
!                    pforce_global(tcx(i),tcy(i),tcz(i),tc0(i),m,3)
!      enddo
       exit
      endif
        enddo
          do k=1,27
           if((posx(k)==tbx(i)).and.(posy(k)==tby(i)).and.&
              (posz(k)==tbz(i)))then
            m1=k
           endif
           if((posx(k)==tcx(i)).and.(posy(k)==tcy(i)).and.&
              (posz(k)==tcz(i)))then
            m2=k
           endif
          enddo
        e3b_qm=e3b_qm+(temp(30)-&
                energy2(tbx(i),tby(i),tbz(i),tb0(i),ta0(i))-&
                energy2(tcx(i),tcy(i),tcz(i),tc0(i),ta0(i))-&
                ee2b_global(m1,m2,tb0(i),tc0(i))+&
        efragment(0,0,0,ta0(i))+efragment(tbx(i),tby(i),tbz(i),tb0(i))+&
              efragment(tcx(i),tcy(i),tcz(i),tc0(i)))/3.d0
         endif
        call system("rm "//cmd(i)(1:32)//".*")
       endif

       close(108)
      enddo
      
!==============================================================================

      call MPI_Allreduce(tforce,force,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(e3b_qm,e3b2_qm,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

          fele=0.0d0
          l_fele=0.d0
          l_fele2=0.d0
          fvdw=0.0d0
          tforce=0.0d0
          ele=0.0d0
          l_ele=0.d0
          l_ele2=0.d0
          evdw=0.0d0
          e2body=0.0d0
          do i=mytid+1,nwat-1,numprocs
           do j=i+1,nwat
            if(.not.connect(0,0,0,i,j))then
             do k=1,3
              do l=1,3
            xa=coord(group(i,k),1)
            ya=coord(group(i,k),2)
            za=coord(group(i,k),3)
            xc=coord(group(j,l),1)
            yc=coord(group(j,l),2)
            zc=coord(group(j,l),3)
            dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
            dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
            ele=ele+charge(group(i,k))*charge(group(j,l))/dist
            l_ele=l_ele+charge(group(i,k))*charge(group(j,l))/dist
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
            l_fele(group(i,k),1)=l_fele(group(i,k),1)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            l_fele(group(i,k),2)=l_fele(group(i,k),2)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            l_fele(group(i,k),3)=l_fele(group(i,k),3)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
            l_fele(group(j,l),1)=l_fele(group(j,l),1)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            l_fele(group(j,l),2)=l_fele(group(j,l),2)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            l_fele(group(j,l),3)=l_fele(group(j,l),3)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
!           iacj=ntype*(iac(group(i,k))-1)
!           ic=ico(iacj+iac(group(j,l)))
!           dinv6=1.0/(dist2*dist2*dist2)
!           dinv12=dinv6**2
!           dinv8=dinv6/dist2
!           dinv14=dinv12/dist2
!           evdw=evdw+cn1(ic)*dinv12-cn2(ic)*dinv6
!           fvdw(group(i,k),1)=fvdw(group(i,k),1)+(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(xa-xc)
!           fvdw(group(i,k),2)=fvdw(group(i,k),2)+(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(ya-yc)
!           fvdw(group(i,k),3)=fvdw(group(i,k),3)+(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(za-zc)
!           fvdw(group(j,l),1)=fvdw(group(j,l),1)-(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(xa-xc)
!           fvdw(group(j,l),2)=fvdw(group(j,l),2)-(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(ya-yc)
!           fvdw(group(j,l),3)=fvdw(group(j,l),3)-(12.0*cn1(ic)*dinv14-&
!                             6.0*cn2(ic)*dinv8)*(za-zc)
              enddo
             enddo
            endif
            if(connect(0,0,0,i,j))then
             do k=1,3
              do l=1,3
            xa=coord(group(i,k),1)
            ya=coord(group(i,k),2)
            za=coord(group(i,k),3)
            xc=coord(group(j,l),1)
            yc=coord(group(j,l),2)
            zc=coord(group(j,l),3)
            dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
            dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
            l_ele=l_ele+charge(group(i,k))*charge(group(j,l))/dist
            l_fele(group(i,k),1)=l_fele(group(i,k),1)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            l_fele(group(i,k),2)=l_fele(group(i,k),2)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            l_fele(group(i,k),3)=l_fele(group(i,k),3)+&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
            l_fele(group(j,l),1)=l_fele(group(j,l),1)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(xa-xc)
            l_fele(group(j,l),2)=l_fele(group(j,l),2)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(ya-yc)
            l_fele(group(j,l),3)=l_fele(group(j,l),3)-&
                 charge(group(i,k))*charge(group(j,l))/(dist**3)*(za-zc)
              enddo
             enddo
            endif
           enddo
          enddo  

      da=0.0
      db=0.0
      dc=0.0
      do j=mytid+1,nwat,numprocs
       do ic=-nxl,nxl
        do jc=-nyl,nyl
         do kc=-nzl,nzl
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
!          do ii=1,natom
!           da=da+field(ic,jc,kc,ii,1)*charge(ii)/18.2223*ic
!           db=db+field(ic,jc,kc,ii,2)*charge(ii)/18.2223*jc
!           dc=dc+field(ic,jc,kc,ii,3)*charge(ii)/18.2223*kc
!          enddo 
           do i=1,nwat
             if(connect(ic,jc,kc,i,j))then
              e2body=e2body+(energy2(ic,jc,kc,i,j)-&
              efragment(ic,jc,kc,i)-efragment(0,0,0,j))*0.5d0
             else
              do ii=1,3
               do jj=1,3
                xa=coord(group(j,jj),1)
                ya=coord(group(j,jj),2)
                za=coord(group(j,jj),3)
                xc=cellcoord(ic,jc,kc,group(i,ii),1)
                yc=cellcoord(ic,jc,kc,group(i,ii),2)
                zc=cellcoord(ic,jc,kc,group(i,ii),3)
                dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
                dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
              ele=ele+0.5d0*charge(group(i,ii))*charge(group(j,jj))/dist
               enddo
              enddo
             endif
             do ii=1,3
              do jj=1,3
               xa=coord(group(j,jj),1)
               ya=coord(group(j,jj),2)
               za=coord(group(j,jj),3)
               xc=cellcoord(ic,jc,kc,group(i,ii),1)
               yc=cellcoord(ic,jc,kc,group(i,ii),2)
               zc=cellcoord(ic,jc,kc,group(i,ii),3)
               dist2=(xa-xc)**2+(ya-yc)**2+(za-zc)**2
               dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
          l_ele=l_ele+0.5d0*charge(group(i,ii))*charge(group(j,jj))/dist
          l_fele(group(j,jj),1)=l_fele(group(j,jj),1)+&
             charge(group(i,ii))*charge(group(j,jj))/(dist**3)*(xa-xc)
          l_fele(group(j,jj),2)=l_fele(group(j,jj),2)+&
             charge(group(i,ii))*charge(group(j,jj))/(dist**3)*(ya-yc)
          l_fele(group(j,jj),3)=l_fele(group(j,jj),3)+&
             charge(group(i,ii))*charge(group(j,jj))/(dist**3)*(za-zc)
              enddo
             enddo
           enddo
          endif
         enddo
        enddo
       enddo
      enddo

      call MPI_Allreduce(l_fele,l_fele2,natom*3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
      call MPI_Allreduce(l_ele,l_ele2,1,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)    
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

          l_fele2=l_fele2*0.5291771d0/627.51
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

       open(96,access='append',file='Dipole.dat')
        write(96,'(4E16.8)')xdipo2,ydipo2,zdipo2,dipo2
       close(96)

!      open(99,access='append',file='out.data') 
       do i=1,nwat
        etotal=etotal+efragment(0,0,0,i)
       enddo
       
2000   format(a11,I4,1x,a12,f16.8)
2001   format(3I2,1x,a4,I4,1x,a12,f16.8)

      do i=1,nwat-1
       do j=i+1,nwat
        if(connect(0,0,0,i,j))then
        energy2body=energy2body+energy2(0,0,0,i,j)-efragment(0,0,0,i)-&
                    efragment(0,0,0,j)  
        endif
       enddo
      enddo

      energy2body=energy2body-ele2/627.51

      open(3001, access='append', file='out_withoutEwald.data')
       write(3001,*)etotal+energy2body+e3b2_qm
      close(3001)

      open(3002,file='force_withoutEwald.dat',access='append')
       do i=1,natom
        write(3002,'(I5,3x,3E16.8)')i,(force(i,k),k=1,3)
       enddo
       write(3002,'(a)')'========'
      close(3002)

      endif
!=========================Long range interaction===================================     
       ele2=0.0
       ele=0.0
       fele=0.0
       tforce=0.0
!     do jj=mytid+1,natom,numprocs
!      do ic=-nxl,nxl
!       do jc=-nyl,nyl
!        do kc=-nzl,nzl
!         if((abs(ic).ne.0).or.(abs(jc).ne.0).or.(abs(kc).ne.0))then
!             do ii=jj,natom
!               xa=coord(jj,1)
!               ya=coord(jj,2)
!               za=coord(jj,3)
!               xc=coord(ii,1)+real(ic)*clata
!               yc=coord(ii,2)+real(jc)*clatb
!               zc=coord(ii,3)+real(kc)*clatc 
!               dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
!               if(ii/=jj)then
!               ele2=ele2+charge(ii)*charge(jj)/dist
!               else
!               ele2=ele2+0.5d0*charge(ii)*charge(jj)/dist
!               endif
!          fele(jj,1)=fele(jj,1)+charge(ii)*charge(jj)/(dist**3)*(xa-xc)
!          fele(jj,2)=fele(jj,2)+charge(ii)*charge(jj)/(dist**3)*(ya-yc)
!          fele(jj,3)=fele(jj,3)+charge(ii)*charge(jj)/(dist**3)*(za-zc)
!          fele(ii,1)=fele(ii,1)-charge(ii)*charge(jj)/(dist**3)*(xa-xc)
!          fele(ii,2)=fele(ii,2)-charge(ii)*charge(jj)/(dist**3)*(ya-yc)
!          fele(ii,3)=fele(ii,3)-charge(ii)*charge(jj)/(dist**3)*(za-zc)
!             enddo
!         endif
!        enddo
!       enddo
!      enddo
!     enddo

!     do jj=mytid+1,natom,numprocs
!      do ii=jj+1,natom
!       if(resnum(ii).ne.resnum(jj))then
!               xa=coord(jj,1)
!               ya=coord(jj,2)
!               za=coord(jj,3)
!               xc=coord(ii,1)
!               yc=coord(ii,2)
!               zc=coord(ii,3)
!               dist=sqrt((xa-xc)**2+(ya-yc)**2+(za-zc)**2)
!               ele2=ele2+charge(ii)*charge(jj)/dist
!          fele(jj,1)=fele(jj,1)+charge(ii)*charge(jj)/(dist**3)*(xa-xc)
!          fele(jj,2)=fele(jj,2)+charge(ii)*charge(jj)/(dist**3)*(ya-yc)
!          fele(jj,3)=fele(jj,3)+charge(ii)*charge(jj)/(dist**3)*(za-zc)
!          fele(ii,1)=fele(ii,1)-charge(ii)*charge(jj)/(dist**3)*(xa-xc)
!          fele(ii,2)=fele(ii,2)-charge(ii)*charge(jj)/(dist**3)*(ya-yc)
!          fele(ii,3)=fele(ii,3)-charge(ii)*charge(jj)/(dist**3)*(za-zc)
!       endif
!      enddo
!     enddo

!      call MPI_Allreduce(fele,tforce,natom*3,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!      call MPI_Allreduce(ele2,ele,1,&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!         tforce=tforce*0.5291771d0/627.51

!     call echarge1d
!     ec=ec/627.51
!     dec=dec*0.5291771d0/627.51

      do i=1,natom
       force(i,1)=force(i,1)-dec(1,i)-l_fele2(i,1)
       force(i,2)=force(i,2)-dec(2,i)-l_fele2(i,2)
       force(i,3)=force(i,3)-dec(3,i)-l_fele2(i,3)
      enddo

2011  format(a11,1x,2I4,a20,2f16.10)
2012  format(3I2,1x,a4,1x,2I4,a20,2f16.10)
      
      etotal=etotal+energy2body-l_ele2/627.51+ec+e3b2_qm

      if(master)then

!    open(2021,file='force_all.dat',access='append')
!     do i=1,natom
!      write(2021,'(I5,3x,3E16.8)')i,(force(i,k),k=1,3)
!     enddo
!     write(2021,'(a)')'========'
!    close(2021)

!    write(99,*)'Total electrical energy of the unit cell (au):',etotal

!      close(99)

      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

       end subroutine
	
