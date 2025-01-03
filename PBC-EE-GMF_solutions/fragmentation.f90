           subroutine fragment_withchg
           use variables
           implicit none
           include 'parallel.h'
           include 'mpif.h'
       
          logical::withion(300),ics,iiod,ilable(100),withion2(300) 
          character(len=10)::ionmask
          integer(kind=4)::igroup(20,10),niongroup,marka,markb,markc

      open(105,file='H2O.inpcrd')
      read(105,*)
      read(105,*)
      read(105,'(6f12.7)')((coord(i,j),j=1,3),i=1,natom)
      read(105,'(3f12.7)')clata,clatb,clatc
      close(105)
 
       do i=-nxm,nxm
         do j=-nym,nym
           do k=-nzm,nzm
              do m=1,natom
                cellcoord(i,j,k,m,1)=coord(m,1)+real(i)*clata
                cellcoord(i,j,k,m,2)=coord(m,2)+real(j)*clatb
                cellcoord(i,j,k,m,3)=coord(m,3)+real(k)*clatc
              enddo
           enddo
         enddo
       enddo

       ilable=.false.
       do i=1,nion/2
        ionpair(i,1)=i
        dist2=999.d0
        do j=nion/2+1,nion
         if(.not.ilable(j))then
          dist=sqrt((coord(i,1)-coord(j,1))**2+&
                    (coord(i,2)-coord(j,2))**2+&
                    (coord(i,3)-coord(j,3))**2)
          if(dist.lt.dist2)then
           dist2=dist
           mark=j
          endif
         endif
        enddo
        ionpair(i,2)=mark
        ilable(mark)=.true.
       enddo

       ilable=.false.
       k=0
       igroup=0
       do i=1, nion/2
        if(.not.ilable(i))then
         ilable(i)=.true.
         k=k+1
         ii=2
         igroup(k,1)=ionpair(i,1)
         igroup(k,2)=ionpair(i,2)
!        do j=i+1, nion/2
!         if(.not.ilable(j))then
!         dist=sqrt((coord(ionpair(i,1),1)-coord(ionpair(j,1),1))**2+&
!                   (coord(ionpair(i,1),2)-coord(ionpair(j,1),2))**2+&
!                   (coord(ionpair(i,1),3)-coord(ionpair(j,1),3))**2)
!         temp(1)=dist
!         dist=sqrt((coord(ionpair(i,1),1)-coord(ionpair(j,2),1))**2+&
!                   (coord(ionpair(i,1),2)-coord(ionpair(j,2),2))**2+&
!                   (coord(ionpair(i,1),3)-coord(ionpair(j,2),3))**2)
!         temp(2)=dist
!         dist=sqrt((coord(ionpair(i,2),1)-coord(ionpair(j,1),1))**2+&
!                   (coord(ionpair(i,2),2)-coord(ionpair(j,1),2))**2+&
!                   (coord(ionpair(i,2),3)-coord(ionpair(j,1),3))**2)
!         temp(3)=dist
!         dist=sqrt((coord(ionpair(i,2),1)-coord(ionpair(j,2),1))**2+&
!                   (coord(ionpair(i,2),2)-coord(ionpair(j,2),2))**2+&
!                   (coord(ionpair(i,2),3)-coord(ionpair(j,2),3))**2)
!         temp(4)=dist
!          if(minval(temp(1:4)).lt.4.0d0)then
!          ii=ii+1
!          igroup(k,ii)=ionpair(j,1)
!          ii=ii+1
!          igroup(k,ii)=ionpair(j,2)
!          ilable(j)=.true.
!          endif
!         endif
!        enddo
        endif
       enddo
       niongroup=k

       withion=.false.
       group=0
       do i=1,niongroup
        withion2=.false.
        do ii=1,10
         if(igroup(i,ii)==0)exit
         dist2=999.0
         do j=1,nwat
          if(.not.withion(j))then
          do k=qmo(j,0),qmo(j,0)+2
          dist=sqrt((coord(igroup(i,ii),1)-coord(k,1))**2+&
                    (coord(igroup(i,ii),2)-coord(k,2))**2+&
                    (coord(igroup(i,ii),3)-coord(k,3))**2)
           if((dist.lt.dist2))then
            dist2=dist
            mark=j
           endif
          enddo
          endif
         enddo
         if(.not.withion2(mark))then
         do j=1,20
          if(group(i,j).eq.0)then
           group(i,j)=qmo(mark,0)
           group(i,j+1)=qmo(mark,0)+1
           group(i,j+2)=qmo(mark,0)+2
           exit
          endif
         enddo
         withion2(mark)=.true.
         endif
        enddo
        do j=1,nwat
         if(withion2(j))then
          withion(j)=.true.
         endif
        enddo
        do ii=1,10
         if(igroup(i,ii)/=0)then
          do j=1,20
           if(group(i,j)==0)then
            group(i,j)=igroup(i,ii)
            exit
           endif
          enddo
         else
          exit
         endif
        enddo
       enddo

       ngroup=niongroup
       do i=1,nwat
        if(.not.withion(i))then
         ngroup=ngroup+1
         group(ngroup,1)=qmo(i,0)
         group(ngroup,2)=qmo(i,0)+1
         group(ngroup,3)=qmo(i,0)+2
        endif
       enddo

      connect=.false.
      onebody=.false.
      do i=1,ngroup-1
       do j=i+1,ngroup
        do ii=1,20
        if(group(i,ii).ne.0)then
         do jj=1,20
         if(group(j,jj).ne.0)then
       dist=sqrt((coord(group(i,ii),1)-coord(group(j,jj),1))**2+&
                 (coord(group(i,ii),2)-coord(group(j,jj),2))**2+&
                 (coord(group(i,ii),3)-coord(group(j,jj),3))**2)
          if(dist.le.r2cut2)then
           connect(0,0,0,i,j)=.true.
          endif
         endif
         enddo
        endif
        enddo
       enddo
      enddo

      do ic=-nxm,nxm
       do jc=-nym,nym
        do kc=-nzm,nzm
         do i=1,ngroup 
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
          do ii=1,20
          if(group(i,ii).ne.0)then
           do j=1,ngroup
            do jj=1,20
            if(group(j,jj).ne.0)then
            dist=sqrt((cellcoord(ic,jc,kc,group(i,ii),1)-&
                                             coord(group(j,jj),1))**2+&
                      (cellcoord(ic,jc,kc,group(i,ii),2)-&
                                             coord(group(j,jj),2))**2+&
                      (cellcoord(ic,jc,kc,group(i,ii),3)-&
                                             coord(group(j,jj),3))**2)
             if(dist.le.r2cut2)then
             connect(ic,jc,kc,i,j)=.true.
             onebody(ic,jc,kc,i)=.true.
             endif
            endif
            enddo
           enddo
          endif 
          enddo
          endif
         enddo
        enddo
       enddo
      enddo

      scal=.false.
      do i=nion+1,natom,3
       do j=nion+1,natom,3
        if(j/=i)then
         dist=sqrt((coord(i,1)-coord(j+1,1))**2+&
                   (coord(i,2)-coord(j+1,2))**2+&
                   (coord(i,3)-coord(j+1,3))**2)
         if(dist.le.1.75d0)then
          scal(j)=.true.
         endif
         dist=sqrt((coord(i,1)-coord(j+2,1))**2+&
                   (coord(i,2)-coord(j+2,2))**2+&
                   (coord(i,3)-coord(j+2,3))**2)
         if(dist.le.1.75d0)then
          scal(j)=.true.
         endif
        endif
       enddo
      enddo

      do i=1,nion
       charge(i)=charge(i)*0.7d0
      enddo
      do i=nion+1,natom,3
       if(scal(i))then
        charge(i) = charge(i) * 0.8d0
        charge(i+1) = charge(i+1) * 0.8d0
        charge(i+2) = charge(i+2) * 0.8d0
       else
        charge(i) = charge(i) * 0.9d0
        charge(i+1) = charge(i+1) * 0.9d0
        charge(i+2) = charge(i+2) * 0.9d0
       endif
      enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================fragmentation==================================
      njob=0
      do i=1,ngroup
       njob=njob+1
       cmd(njob)='000cell'//char(48+i/100)//&
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
            cmd(njob)=char(48+ic+nxm*2)//char(48+jc+nym*2)//&
                     char(48+kc+nzm*2)//'cell'//char(48+i/100)//&
                   char(48+(i-i/100*100)/10)//&
                   char(48+(i-i/10*10))//'.gjf'
           endif
          enddo
         endif
        enddo
       enddo
      enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)

      do ni=mytid+1,njob,numprocs
      if(cmd(ni)(1:3).eq.'000')then

       read(cmd(ni),'(7x,I3)')i
       ncharge=0
       bgc=.true.
       open(6,file=trim(adjustl(cmd(ni))),status='replace')
       open(7,file='field_'//trim(adjustl(cmd(ni))),status='replace')
       write(6,'(a)')'%chk=000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.chk'
       write(7,'(a)')'%chk=000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                         char(48+(i-i/10*10))//'.chk'
       write(6,'(a)')'%nproc=4'
       write(7,'(a)')'%nproc=4' 
       write(6,'(a)')'%mem=5GB' 
       write(7,'(a)')'%mem=5GB' 
       write(6,'(a)')'#p MP2(FULL)/aug-cc-pvdz nosymm charge force'//&
                     ' density'
       write(7,'(a)')'#p density=(check,mp2) nosymm geom=allcheck '//&
                 'guess=(read,only) prop=(read,field)'
       write(6,'(a)')
       write(6,'(a)')'Have a nice day'
       write(6,'(a)')
       write(7,'(a)')
       do k=1,20
        if(group(i,k)==0)then
         exit
        else
         if(resname(group(i,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(i,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do j=1,20
        if(group(i,j).ne.0)then
        write(6,1004)atomname(group(i,j))(2:3),(coord(group(i,j),k),k=1,3)
        bgc(group(i,j))=.false.
        else
         exit
        endif
       enddo
1002  format(a1,4x,3f14.8)    
1004  format(a2,3x,3f14.8)  
       write(6,*)
       do j=1,natom
        if(bgc(j))then
         write(6,1003)coord(j,1),coord(j,2),coord(j,3),charge(j)
         write(7,'(5x,3f14.8)')coord(j,1),coord(j,2),coord(j,3)
        endif
       enddo
1003  format(3x,3f14.8,f16.8)
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do j=1,natom
            write(6,1003)cellcoord(ic,jc,kc,j,1),&
            cellcoord(ic,jc,kc,j,2),cellcoord(ic,jc,kc,j,3),charge(j)
            write(7,'(5x,3f14.8)')cellcoord(ic,jc,kc,j,1),&
            cellcoord(ic,jc,kc,j,2),cellcoord(ic,jc,kc,j,3)
           enddo
          endif
         enddo
        enddo
       enddo
       write(6,*)
       write(7,*)
       close(6)
       close(7)

      else

         read(cmd(ni),'(3I1,4x,I3)')ii,jj,kk,i
         ic=ii-nxm*2
         jc=jj-nym*2
         kc=kk-nzm*2
         ncharge=0
         bgc=.true.
         open(6,file=trim(adjustl(cmd(ni))),status='replace')
!!       open(7,file='field_'//trim(adjustl(cmd(ni))),status='replace')
         write(6,'(a)')'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
                       'cell'//char(48+i/100)//&
                     char(48+(i-i/100*100)/10)//&
                     char(48+(i-i/10*10))//'.chk'
!!       write(7,*)'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
!!                 'cell'//char(48+i/100)//&
!!                   char(48+(i-i/100*100)/10)//&
!!                   char(48+(i-i/10*10))//'.chk'
         write(6,'(a)')'%nproc=4'
!!       write(7,*)'%nproc=4'
         write(6,'(a)')'%mem=5GB'
!!       write(7,*)'%mem=5GB'
         write(6,'(a)')'#p MP2(FULL)/aug-cc-pVDZ nosymm charge '
!!       write(7,*)'#p density=(check,mp2) nosymm geom=allcheck '//&
!!                 'guess=(read,only) prop=(read,field)'
             write(6,'(a)')
!!           write(7,*)
             write(6,'(a)')'Have a nice day'
             write(6,'(a)')
       do k=1,20
        if(group(i,k)==0)then
         exit
        else
         if(resname(group(i,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(i,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=1,20
        if(group(i,l).ne.0)then
        write(6,1004)atomname(group(i,l))(2:3),&
                     (cellcoord(ic,jc,kc,group(i,l),k),k=1,3)
        bgc(group(i,l))=.false.
        else
         exit
        endif
       enddo
       write(6,*)
       do ic1=-nxm,nxm
        do jc1=-nym,nym
         do kc1=-nzm,nzm
          do l=1,natom
           if((ic1.ne.ic).or.(jc1.ne.jc).or.(kc1.ne.kc).or.&
             bgc(l))then
            write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
!!          write(7,'(5x,3f14.8)')cellcoord(ic1,jc1,kc1,l,1),&
!!          cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3)
           endif
          enddo
         enddo
        enddo
       enddo
       write(6,*)
!!     write(7,*)
       close(6)
!!     close(7)
       
      endif
      enddo
!=======================================================================
!=========================two body interaction==========================
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
                char(48+(j-j/10*10))//'.gjf'
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
                char(48+(j-j/10*10))//'.gjf'
            endif
           enddo
          enddo
         endif
        enddo
       enddo
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      do ni=mytid+1,njob,numprocs
      if(cmd(ni)(1:3).eq.'000')then
       
       read(cmd(ni),'(7x,I3,8x,I3)')i,j
         ncharge=0
         bgc=.true.
         open(6,file=trim(adjustl(cmd(ni))),status='replace')
         open(7,file='field_'//trim(adjustl(cmd(ni))),status='replace')
         write(6,'(a)')'%chk=000cell'//char(48+i/100)//& 
                   char(48+(i-i/100*100)/10)//& 
           char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                            char(48+(j-j/100*100)/10)//&
                            char(48+(j-j/10*10))//'.chk'
         write(7,'(a)')'%chk=000cell'//char(48+i/100)//& 
                   char(48+(i-i/100*100)/10)//& 
           char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                            char(48+(j-j/100*100)/10)//&
                            char(48+(j-j/10*10))//'.chk'
         write(6,'(a)')'%nproc=4'
         write(7,'(a)')'%nproc=4'
         write(6,'(a)')'%mem=5GB'
         write(7,'(a)')'%mem=5GB'
         write(6,'(a)')'#p MP2(FULL)/aug-cc-pVDZ nosymm charge force'
         write(7,'(a)')'#p density=(check,mp2) nosymm geom=allcheck '//&
                   'guess=(read,only) prop=(read,field)'
         write(6,'(a)')
         write(6,'(a)')'Have a nice day'
         write(6,'(a)')
         write(7,'(a)')
       do k=1,20
        if(group(i,k)==0)then
         exit
        else
         if(resname(group(i,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(i,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       do k=1,20
        if(group(j,k)==0)then
         exit
        else
         if(resname(group(j,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(j,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=1,20
        if(group(i,l).ne.0)then
        write(6,1004)atomname(group(i,l))(2:3),(coord(group(i,l),k),k=1,3)
        bgc(group(i,l))=.false.
        else
         exit
        endif
       enddo
       do l=1,20
        if(group(j,l).ne.0)then
        write(6,1004)atomname(group(j,l))(2:3),(coord(group(j,l),k),k=1,3)
        bgc(group(j,l))=.false.
        else
         exit
        endif
       enddo
       write(6,*)
       do l=1,natom
        if(bgc(l))then
         write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
         write(7,'(5x,3f14.8)')coord(l,1),coord(l,2),coord(l,3) 
        endif
       enddo
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do l=1,natom
            write(6,1003)cellcoord(ic,jc,kc,l,1),&
            cellcoord(ic,jc,kc,l,2),cellcoord(ic,jc,kc,l,3),charge(l)
            write(7,'(5x,3f14.8)')cellcoord(ic,jc,kc,l,1),&
            cellcoord(ic,jc,kc,l,2),cellcoord(ic,jc,kc,l,3)
           enddo
          endif
         enddo
        enddo
       enddo
       write(6,*)
       write(7,*)
       close(6)
       close(7)
     
      else

        read(cmd(ni),'(3I1,4x,I3,8x,I3)')ii,jj,kk,i,j
        ic=ii-nxm*2
        jc=jj-nym*2
        kc=kk-nzm*2
        ncharge=0
        bgc=.true.
        bgc2=.true.
        open(6,file=trim(adjustl(cmd(ni))),status='replace')
!       open(7,file='field_'//trim(adjustl(cmd(ni))),status='replace')
        write(6,'(a)')'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
               'cell'//char(48+i/100)//char(48+(i-i/100*100)/10)//&
            char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                       char(48+(j-j/100*100)/10)//&
                       char(48+(j-j/10*10))//'.chk'
!       write(7,*)'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
!              'cell'//char(48+i/100)//char(48+(i-i/100*100)/10)//&
!           char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
!                      char(48+(j-j/100*100)/10)//&
!                      char(48+(j-j/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
!       write(7,*)'%nproc=4'
        write(6,'(a)')'%mem=5GB'
!       write(7,*)'%mem=5GB'
       write(6,'(a)')'#p MP2(FULL)/aug-cc-pVDZ nosymm charge force'
!      write(7,*)'#p density=(check,mp2) nosymm geom=allcheck '//&
!                'guess=(read,only) prop=(read,field)'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
!       write(7,*)
       do k=1,20
        if(group(i,k)==0)then
         exit
        else
         if(resname(group(i,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(i,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       do k=1,20
        if(group(j,k)==0)then
         exit
        else
         if(resname(group(j,k)).eq.'Na+')then
          ncharge=ncharge+1
         endif
         if(resname(group(j,k)).eq.'Cl-')then
          ncharge=ncharge-1
         endif
        endif
       enddo
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=1,20
        if(group(j,l).ne.0)then
        write(6,1004)atomname(group(j,l))(2:3),(coord(group(j,l),k),k=1,3)
        bgc(group(j,l))=.false.
        else
         exit
        endif
       enddo
       do l=1,20
        if(group(i,l).ne.0)then
        write(6,1004)atomname(group(i,l))(2:3),&
                     (cellcoord(ic,jc,kc,group(i,l),k),k=1,3)
        bgc2(group(i,l))=.false.
        else
         exit
        endif
       enddo
       write(6,*)
       do l=1,natom
        if(bgc(l))then
         write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
!        write(7,'(5x,3f14.8)')coord(l,1),coord(l,2),coord(l,3)
        endif
       enddo
       do ic1=-nxm,nxm
        do jc1=-nym,nym
         do kc1=-nzm,nzm
          if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
           do l=1,natom
            if((ic.ne.ic1).or.(jc.ne.jc1).or.(kc.ne.kc1).or.bgc2(l))then
             write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
!            write(7,'(5x,3f14.8)')cellcoord(ic1,jc1,kc1,l,1),&
!           cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3)
            endif
           enddo
          endif
         enddo
        enddo
       enddo
       write(6,*)
!      write(7,*)
       close(6)
!      close(7)
        

      endif
      enddo

!=============================================================================

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

         end subroutine 
