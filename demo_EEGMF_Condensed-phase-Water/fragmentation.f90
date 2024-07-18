           subroutine fragment_withchg
           use variables
           use ewald
           implicit none
           include 'parallel.h'
           include 'mpif.h'

           integer(kind=4)::ic2,jc2,kc2
          
      open(105,file='H2O.inpcrd')
      read(105,*)
      read(105,*)
      read(105,'(6f12.7)')((coord(i,j),j=1,3),i=1,natom)
      read(105,'(3f12.7)')clata,clatb,clatc
      close(105)
 
       m=0
       do i=-1,1
         do j=-1,1
           do k=-1,1
            m=m+1
            posx(m)=i
            posy(m)=j
            posz(m)=k
           enddo
         enddo
       enddo

       do i=-nxl,nxl
         do j=-nyl,nyl
           do k=-nzl,nzl
              do m=1,natom
                cellcoord(i,j,k,m,1)=coord(m,1)+real(i)*clata
                cellcoord(i,j,k,m,2)=coord(m,2)+real(j)*clatb
                cellcoord(i,j,k,m,3)=coord(m,3)+real(k)*clatc
              enddo
           enddo
         enddo
       enddo

      connect=.false.
      connect3=.false.
      onebody=.false.
      do i=1,nwat-1
       do j=i+1,nwat
        do ii=1,3
         do jj=1,3
       dist=sqrt((coord(group(i,ii),1)-coord(group(j,jj),1))**2+&
                 (coord(group(i,ii),2)-coord(group(j,jj),2))**2+&
                 (coord(group(i,ii),3)-coord(group(j,jj),3))**2)
          if(dist.le.r2cut2)then
           connect(0,0,0,i,j)=.true.
          endif
          if(dist.le.rcut3)then
           connect3(0,0,0,i,j)=.true.
          endif
         enddo
        enddo
       enddo
      enddo

      do j=1,nwat
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          do i=j,nwat 
           if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
            do ii=1,3
             do jj=1,3
            dist=sqrt((cellcoord(ic,jc,kc,group(i,ii),1)-&
                                             coord(group(j,jj),1))**2+&
                      (cellcoord(ic,jc,kc,group(i,ii),2)-&
                                             coord(group(j,jj),2))**2+&
                      (cellcoord(ic,jc,kc,group(i,ii),3)-&
                                             coord(group(j,jj),3))**2)
             if(dist.le.r2cut2)then
              if(i/=j)then
               connect(ic,jc,kc,i,j)=.true.
               onebody(ic,jc,kc,i)=.true.
                dnx=0-ic
                dny=0-jc
                dnz=0-kc
                ic1=0+dnx
                jc1=0+dny
                kc1=0+dnz
                connect(ic1,jc1,kc1,j,i)=.true.
                onebody(ic1,jc1,kc1,j)=.true.
               if(dist.le.rcut3)then
                connect3(ic,jc,kc,i,j)=.true.
                dnx=0-ic
                dny=0-jc
                dnz=0-kc
                ic1=0+dnx
                jc1=0+dny
                kc1=0+dnz
                connect3(ic1,jc1,kc1,j,i)=.true.
               endif
              else
               if(.not.connect(-ic,-jc,-kc,i,j))then
                connect(ic,jc,kc,i,j)=.true.
                onebody(ic,jc,kc,i)=.true.
                dnx=0-ic
                dny=0-jc
                dnz=0-kc
                ic1=0+dnx
                jc1=0+dny
                kc1=0+dnz
                connect(ic1,jc1,kc1,j,i)=.true.
                onebody(ic1,jc1,kc1,j)=.true.
                if(dist.le.rcut3)then
                 connect3(ic,jc,kc,i,j)=.true.
                 dnx=0-ic
                 dny=0-jc
                 dnz=0-kc
                 ic1=0+dnx
                 jc1=0+dny
                 kc1=0+dnz
                 connect3(ic1,jc1,kc1,j,i)=.true.
                endif
               endif
              endif
             endif
             enddo
            enddo
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

      n3b=0
      do i=1,nwat-2
       do j=i+1,nwat-1
        do k=j+1,nwat
         if((connect3(0,0,0,i,j).and.connect3(0,0,0,i,k).and.&
             connect(0,0,0,j,k)).or.(connect3(0,0,0,i,j).and.&
             connect(0,0,0,i,k).and.connect3(0,0,0,j,k)).or.&
            (connect(0,0,0,i,j).and.connect3(0,0,0,i,k).and.&
             connect3(0,0,0,j,k)))then
          n3b=n3b+1
          ta0(n3b)=i
          tb0(n3b)=j
          tbx(n3b)=0
          tby(n3b)=0
          tbz(n3b)=0
          tc0(n3b)=k
          tcx(n3b)=0
          tcy(n3b)=0
          tcz(n3b)=0
         endif
        enddo
       enddo
      enddo

      do i=1,nwat-1
       do j=i+1,nwat
        if(connect(0,0,0,i,j))then
         do ic=-nxm,nxm
          do jc=-nym,nym
           do kc=-nzm,nzm
            if((ic/=0).or.(jc/=0).or.(kc/=0))then
             do k=1,nwat
              if((connect(0,0,0,i,j).and.connect3(ic,jc,kc,k,i).and.&
                  connect3(ic,jc,kc,k,j)).or.(connect3(0,0,0,i,j).and.&
                  connect(ic,jc,kc,k,i).and.connect3(ic,jc,kc,k,j)).or.&
                 (connect3(0,0,0,i,j).and.connect3(ic,jc,kc,k,i).and.&
                  connect(ic,jc,kc,k,j)))then
               n3b=n3b+1
               ta0(n3b)=i
               tb0(n3b)=j
               tbx(n3b)=0
               tby(n3b)=0
               tbz(n3b)=0
               tc0(n3b)=k
               tcx(n3b)=ic
               tcy(n3b)=jc
               tcz(n3b)=kc
              endif
             enddo
            endif
           enddo
          enddo
         enddo
        endif
       enddo
      enddo

      do i=1,nwat
       do ic=-nxm,nxm
        do jc=-nym,nym
         do kc=-nzm,nzm
          if((ic/=0).or.(jc/=0).or.(kc/=0))then
          do j=1,nwat
           do k=j+1,nwat
           if((connect(ic,jc,kc,j,i).and.connect3(ic,jc,kc,k,i).and.&
              connect3(0,0,0,j,k)).or.(connect3(ic,jc,kc,j,i).and.&
              connect(ic,jc,kc,k,i).and.connect3(0,0,0,j,k)).or.&
             (connect3(ic,jc,kc,j,i).and.connect3(ic,jc,kc,k,i).and.&
              connect(0,0,0,j,k)))then
                 n3b=n3b+1
                 ta0(n3b)=i
                 tbx(n3b)=ic
                 tby(n3b)=jc
                 tbz(n3b)=kc
                 tb0(n3b)=j
                 tcx(n3b)=ic
                 tcy(n3b)=jc
                 tcz(n3b)=kc
                 tc0(n3b)=k
           endif
           enddo
          enddo
          endif
         enddo
        enddo
       enddo
      enddo

      do i=1,nwat
       do ii=1,26
        ic=posx(ii)
        jc=posy(ii)
        kc=posz(ii)
          if((ic/=0).or.(jc/=0).or.(kc/=0))then
          do j=1,nwat
           if(connect3(ic,jc,kc,j,i))then
            do jj=ii+1,27
             ic1=posx(jj)
             jc1=posy(jj)
             kc1=posz(jj)
               if((ic1/=0).or.(jc1/=0).or.(kc1/=0))then
               do k=1,nwat
                if(connect3(ic1,jc1,kc1,k,i))then
                 dnx=0-ic
                 dny=0-jc
                 dnz=0-kc
                 ic2=ic1+dnx
                 jc2=jc1+dny
                 kc2=kc1+dnz
                 if(connect(ic2,jc2,kc2,k,j))then
                 n3b=n3b+1
                 ta0(n3b)=i
                 tbx(n3b)=ic
                 tby(n3b)=jc
                 tbz(n3b)=kc
                 tb0(n3b)=j
                 tcx(n3b)=ic1
                 tcy(n3b)=jc1
                 tcz(n3b)=kc1
                 tc0(n3b)=k
                 endif
                endif
               enddo
               endif
            enddo
           endif
          enddo
          endif
       enddo
      enddo

      nr=0
      x2b=.false.
       do j=1,n3b
        if((tbx(j)/=0).or.(tby(j)/=0).or.(tbz(j)/=0))then
         if((tcx(j)/=0).or.(tcy(j)/=0).or.(tcz(j)/=0))then
          do k=1,27
         if((posx(k)==tbx(j)).and.(posy(k)==tby(j)).and.&
            (posz(k)==tbz(j)))then
          m1=k
         endif
         if((posx(k)==tcx(j)).and.(posy(k)==tcy(j)).and.&
            (posz(k)==tcz(j)))then
          m2=k
         endif
          enddo
          if(.not.x2b(m1,m2,tb0(j),tc0(j)))then
          x2b(m1,m2,tb0(j),tc0(j))=.true.
          nr=nr+1
          rbx(nr)=tbx(j)
          rby(nr)=tby(j)
          rbz(nr)=tbz(j)
          rb0(nr)=tb0(j)
          rcx(nr)=tcx(j)
          rcy(nr)=tcy(j)
          rcz(nr)=tcz(j)
          rc0(nr)=tc0(j)
          endif
         endif
        endif
       enddo

       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================fragmentation==================================
      njob=0
      do i=1,nwat
       njob=njob+1
       cmd(njob)='000cell'//char(48+i/100)//&
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
       write(6,'(a)')'%mem=10GB' 
       write(7,'(a)')'%mem=10GB' 
       write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force'//&
                     ' density'
       write(7,'(a)')'#p chkbasis density=(check,cc) nosymm'//&
                 ' geom=allcheck '//&
                 'guess=(read,only) prop=(read,field)'
       write(6,'(a)')
       write(6,'(a)')'Have a nice day'
       write(6,'(a)')
       write(7,'(a)')
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do j=group(i,1),group(i,3)
         write(6,1002)atomname(j)(2:2),(coord(j,k),k=1,3)
         bgc(j)=.false.
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
       do ic=-nxl,nxl
        do jc=-nyl,nyl
         do kc=-nzl,nzl
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do j=1,natom
            write(6,1003)cellcoord(ic,jc,kc,j,1),&
            cellcoord(ic,jc,kc,j,2),cellcoord(ic,jc,kc,j,3),charge(j)
!           write(7,'(5x,3f14.8)')cellcoord(ic,jc,kc,j,1),&
!           cellcoord(ic,jc,kc,j,2),cellcoord(ic,jc,kc,j,3)
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
         write(6,'(a)')'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
                       'cell'//char(48+i/100)//&
                     char(48+(i-i/100*100)/10)//&
                     char(48+(i-i/10*10))//'.chk'
         write(6,'(a)')'%nproc=4'
         write(6,'(a)')'%mem=10GB'
         write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
         write(6,'(a)')
         write(6,'(a)')'Have a nice day'
         write(6,'(a)')
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=group(i,1),group(i,3)
         write(6,1002)atomname(l)(2:2),&
                     (cellcoord(ic,jc,kc,l,k),k=1,3)
         bgc(l)=.false.
       enddo
       write(6,*)
       do ic1=-nxl,nxl
        do jc1=-nyl,nyl
         do kc1=-nzl,nzl
          do l=1,natom
           if((ic1.ne.ic).or.(jc1.ne.jc).or.(kc1.ne.kc).or.&
             bgc(l))then
            write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
           endif
          enddo
         enddo
        enddo
       enddo
       write(6,*)
       close(6)
       
      endif
      enddo
!=======================================================================
!=========================two body interaction==========================
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
                char(48+(j-j/10*10))//'.gjf'
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
                char(48+(j-j/10*10))//'.gjf'
            endif
           enddo
          endif
         enddo
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
         write(6,'(a)')'%mem=10GB'
         write(7,'(a)')'%mem=10GB'
         write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
         write(7,'(a)')'#p chkbasis density=(check,cc) nosymm'//&
                       ' geom=allcheck '//&
                   'guess=(read,only) prop=(read,field)'
         write(6,'(a)')
         write(6,'(a)')'Have a nice day'
         write(6,'(a)')
         write(7,'(a)')
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=group(i,1),group(i,3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
       enddo
       do l=group(j,1),group(j,3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
       enddo
       write(6,*)
       do l=1,natom
        if(bgc(l))then
         write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
         write(7,'(5x,3f14.8)')coord(l,1),coord(l,2),coord(l,3) 
        endif
       enddo
       do ic=-nxl,nxl
        do jc=-nyl,nyl
         do kc=-nzl,nzl
          if((ic.ne.0).or.(jc.ne.0).or.(kc.ne.0))then
           do l=1,natom
            write(6,1003)cellcoord(ic,jc,kc,l,1),&
            cellcoord(ic,jc,kc,l,2),cellcoord(ic,jc,kc,l,3),charge(l)
!           write(7,'(5x,3f14.8)')cellcoord(ic,jc,kc,l,1),&
!           cellcoord(ic,jc,kc,l,2),cellcoord(ic,jc,kc,l,3)
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
        write(6,'(a)')'%chk='//char(48+ii)//char(48+jj)//char(48+kk)//&
               'cell'//char(48+i/100)//char(48+(i-i/100*100)/10)//&
            char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                       char(48+(j-j/100*100)/10)//&
                       char(48+(j-j/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
        write(6,'(a)')'%mem=10GB'
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
       write(6,'(2x,I2,a4)')ncharge,'  1 '
       do l=group(j,1),group(j,3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
       enddo
       do l=group(i,1),group(i,3)
       write(6,1002)atomname(l)(2:2),(cellcoord(ic,jc,kc,l,k),k=1,3)
       bgc2(l)=.false.
       enddo
       write(6,*)
       do l=1,natom
        if(bgc(l))then
         write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
        endif
       enddo
       do ic1=-nxl,nxl
        do jc1=-nyl,nyl
         do kc1=-nzl,nzl
          if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
           do l=1,natom
            if((ic.ne.ic1).or.(jc.ne.jc1).or.(kc.ne.kc1).or.bgc2(l))then
             write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
            endif
           enddo
          endif
         enddo
        enddo
       enddo
       write(6,*)
       close(6)

      endif
      enddo

!=============================================================================

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

!===============================3B interaction================================

      do ni=mytid+1,nr,numprocs
        ncharge=0
        bgc2=.true.
        bgc3=.true.
       open(6,file=char(48+rbx(ni)+2)//char(48+rby(ni)+2)//&
                    char(48+rbz(ni)+2)//'cell'//char(48+rb0(ni)/100)//&
                               char(48+(rb0(ni)-rb0(ni)/100*100)/10)//&
                               char(48+(rb0(ni)-rb0(ni)/10*10))//&
                   '-'//char(48+rcx(ni)+2)//char(48+rcy(ni)+2)//&
                    char(48+rcz(ni)+2)//'cell'//char(48+rc0(ni)/100)//&
                               char(48+(rc0(ni)-rc0(ni)/100*100)/10)//&
                               char(48+(rc0(ni)-rc0(ni)/10*10))//'.gjf',status='replace')
        write(6,'(a)')'%chk='//char(48+rbx(ni)+2)//char(48+rby(ni)+2)//&
                    char(48+rbz(ni)+2)//'cell'//char(48+rb0(ni)/100)//&
                               char(48+(rb0(ni)-rb0(ni)/100*100)/10)//&
                               char(48+(rb0(ni)-rb0(ni)/10*10))//&
                   '-'//char(48+rcx(ni)+2)//char(48+rcy(ni)+2)//&
                    char(48+rcz(ni)+2)//'cell'//char(48+rc0(ni)/100)//&
                               char(48+(rc0(ni)-rc0(ni)/100*100)/10)//&
                               char(48+(rc0(ni)-rc0(ni)/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
        write(6,'(a)')'%mem=10GB'
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
        write(6,'(2x,I2,a4)')ncharge,'  1 '
        do l=group(rb0(ni),1),group(rb0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(rbx(ni),rby(ni),rbz(ni),l,k),k=1,3)
         bgc2(l)=.false.
        enddo
        do l=group(rc0(ni),1),group(rc0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(rcx(ni),rcy(ni),rcz(ni),l,k),k=1,3)
         bgc3(l)=.false.
        enddo
        write(6,*)
        do ic1=-nxl,nxl
         do jc1=-nyl,nyl
          do kc1=-nzl,nzl
            do l=1,natom
             if((rbx(ni).ne.ic1).or.(rby(ni).ne.jc1).or.&
                (rbz(ni).ne.kc1).or.bgc2(l))then
              if((rcx(ni).ne.ic1).or.(rcy(ni).ne.jc1).or.&
                 (rcz(ni).ne.kc1).or.bgc3(l))then
               write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
              endif
             endif
            enddo
          enddo
         enddo
        enddo
        write(6,*)
        close(6)
      enddo

      do ni=mytid+1,n3b,numprocs
       
        ncharge=0
        bgc=.true.
        bgc2=.true.
        bgc3=.true.

        if((tbx(ni)==0).and.(tby(ni)==0).and.(tbz(ni)==0).and.&
           (tcx(ni)==0).and.(tcy(ni)==0).and.(tcz(ni)==0))then
        open(6,file='000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-000cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-000cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.gjf',status='replace')
        write(6,'(a)')'%chk='//'000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-000cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-000cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
        write(6,'(a)')'%mem=10GB'
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
        write(6,'(2x,I2,a4)')ncharge,'  1 '
        do l=group(ta0(ni),1),group(ta0(ni),3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
        enddo
        do l=group(tb0(ni),1),group(tb0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tbx(ni),tby(ni),tbz(ni),l,k),k=1,3)
         bgc(l)=.false.
        enddo
        do l=group(tc0(ni),1),group(tc0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tcx(ni),tcy(ni),tcz(ni),l,k),k=1,3)
         bgc(l)=.false.
        enddo
        write(6,*)
        do l=1,natom
         if(bgc(l))then
          write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
         endif
        enddo
        do ic1=-nxl,nxl
         do jc1=-nyl,nyl
          do kc1=-nzl,nzl
           if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
            do l=1,natom
             write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
            enddo
           endif
          enddo
         enddo
        enddo
        write(6,*)
        close(6)
        endif

        if((tbx(ni)==0).and.(tby(ni)==0).and.(tbz(ni)==0))then
         if((tcx(ni)/=0).or.(tcy(ni)/=0).or.(tcz(ni)/=0))then
        open(6,file='000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-000cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-'//char(48+tcx(ni)+2)//char(48+tcy(ni)+2)//&
                    char(48+tcz(ni)+2)//'cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.gjf',status='replace')
        write(6,'(a)')'%chk='//'000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-000cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-'//char(48+tcx(ni)+2)//char(48+tcy(ni)+2)//&
                    char(48+tcz(ni)+2)//'cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
        write(6,'(a)')'%mem=10GB'
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
        write(6,'(2x,I2,a4)')ncharge,'  1 '
        do l=group(ta0(ni),1),group(ta0(ni),3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
        enddo
        do l=group(tb0(ni),1),group(tb0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tbx(ni),tby(ni),tbz(ni),l,k),k=1,3)
         bgc(l)=.false.
        enddo
        do l=group(tc0(ni),1),group(tc0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tcx(ni),tcy(ni),tcz(ni),l,k),k=1,3)
         bgc2(l)=.false.
        enddo
        write(6,*)
        do l=1,natom
         if(bgc(l))then
          write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
         endif
        enddo
        do ic1=-nxl,nxl
         do jc1=-nyl,nyl
          do kc1=-nzl,nzl
           if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
            do l=1,natom
             if((tcx(ni).ne.ic1).or.(tcy(ni).ne.jc1).or.&
                (tcz(ni).ne.kc1).or.bgc2(l))then
              write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
             endif
            enddo
           endif
          enddo
         enddo
        enddo
        write(6,*)
        close(6)
         endif
        endif

        if((tbx(ni)/=0).or.(tby(ni)/=0).or.(tbz(ni)/=0))then
         if((tcx(ni)/=0).or.(tcy(ni)/=0).or.(tcz(ni)/=0))then
        open(6,file='000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-'//char(48+tbx(ni)+2)//char(48+tby(ni)+2)//&
                    char(48+tbz(ni)+2)//'cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-'//char(48+tcx(ni)+2)//char(48+tcy(ni)+2)//&
                    char(48+tcz(ni)+2)//'cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.gjf',status='replace')
        write(6,'(a)')'%chk='//'000cell'//char(48+ta0(ni)/100)//&
                               char(48+(ta0(ni)-ta0(ni)/100*100)/10)//&
                               char(48+(ta0(ni)-ta0(ni)/10*10))//&
                   '-'//char(48+tbx(ni)+2)//char(48+tby(ni)+2)//&
                    char(48+tbz(ni)+2)//'cell'//char(48+tb0(ni)/100)//&
                               char(48+(tb0(ni)-tb0(ni)/100*100)/10)//&
                               char(48+(tb0(ni)-tb0(ni)/10*10))//&
                   '-'//char(48+tcx(ni)+2)//char(48+tcy(ni)+2)//&
                    char(48+tcz(ni)+2)//'cell'//char(48+tc0(ni)/100)//&
                               char(48+(tc0(ni)-tc0(ni)/100*100)/10)//&
                               char(48+(tc0(ni)-tc0(ni)/10*10))//'.chk'
        write(6,'(a)')'%nproc=4'
        write(6,'(a)')'%mem=10GB'
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm charge force density'
        write(6,'(a)')
        write(6,'(a)')'Have a nice day'
        write(6,'(a)')
        write(6,'(2x,I2,a4)')ncharge,'  1 '
        do l=group(ta0(ni),1),group(ta0(ni),3)
         write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
         bgc(l)=.false.
        enddo
        do l=group(tb0(ni),1),group(tb0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tbx(ni),tby(ni),tbz(ni),l,k),k=1,3)
         bgc2(l)=.false.
        enddo
        do l=group(tc0(ni),1),group(tc0(ni),3)
         write(6,1002)atomname(l)(2:2),&
                      (cellcoord(tcx(ni),tcy(ni),tcz(ni),l,k),k=1,3)
         bgc3(l)=.false.
        enddo
        write(6,*)
        do l=1,natom
        if(bgc(l))then
         write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
        endif
        enddo
        do ic1=-nxl,nxl
         do jc1=-nyl,nyl
          do kc1=-nzl,nzl
           if((ic1.ne.0).or.(jc1.ne.0).or.(kc1.ne.0))then
            do l=1,natom
             if((tbx(ni).ne.ic1).or.(tby(ni).ne.jc1).or.&
                (tbz(ni).ne.kc1).or.bgc2(l))then
              if((tcx(ni).ne.ic1).or.(tcy(ni).ne.jc1).or.&
                 (tcz(ni).ne.kc1).or.bgc3(l))then
               write(6,1003)cellcoord(ic1,jc1,kc1,l,1),&
         cellcoord(ic1,jc1,kc1,l,2),cellcoord(ic1,jc1,kc1,l,3),charge(l)
              endif
             endif
            enddo
           endif
          enddo
         enddo
        enddo
        write(6,*)
        close(6)
         endif
        endif

      enddo

!=============================================================================

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

         end subroutine 
