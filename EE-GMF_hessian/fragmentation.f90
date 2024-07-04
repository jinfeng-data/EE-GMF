           subroutine fragment_withchg
           use variables
           implicit none
           include 'parallel.h'
           include 'mpif.h'
          
           real(kind=8)::pr(natom*3+3)

      existence=.False.
      inquire(file='min.rst',exist=existence)
      if(existence)then
      open(15,file='min.rst')
      read(15,*)
      read(15,*)
      read(15,'(6f15.8)')((coord(i,j),j=1,3),i=1,natom)
      close(15)
      else
      open(15,file='cluster.inpcrd')
      read(15,*)
      read(15,*)
      read(15,'(6f15.8)')((coord(i,j),j=1,3),i=1,natom)
      close(15)
      endif
 
       call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================fragmentation==================================
      do i=mytid+1,nres,numprocs
       open(6,file='000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                char(48+(i-i/100*100-(i-i/100*100)/10*10))//'.gjf',status='replace')
       open(7,file='field_000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
                char(48+(i-i/100*100-(i-i/100*100)/10*10))//'.gjf',status='replace')
       write(6,'(a)')'%chk=000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
               char(48+(i-i/100*100-(i-i/100*100)/10*10))//'.chk'
       write(7,'(a)')'%chk=000cell'//char(48+i/100)//&
                         char(48+(i-i/100*100)/10)//&
               char(48+(i-i/100*100-(i-i/100*100)/10*10))//'.chk'
       write(6,'(a)')'%nproc=4'
       write(7,'(a)')'%nproc=4' 
       write(6,'(a)')'%mem=5GB' 
       write(7,'(a)')'%mem=5GB' 
       write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm SCF=tight charge'//&
                 ' freq'
       write(7,'(a)')'#p density=(check,ccd) nosymm geom=allcheck '//&
                 'guess=(read,only) prop=(read,field)'
       write(6,*)
       write(6,'(a)')'Have a nice day'
       write(6,*)
       write(7,*)
       if(i /= 1)then
        write(6,'(a7)')'  0  1 '
       else
        write(6,'(a7)')'  1  1 '
       endif
       do j=1,natom
        if(resnum(j)==i)then
         write(6,1002)atomname(j)(2:2),(coord(j,k),k=1,3)
        endif
       enddo
1002  format(a1,4x,3f14.8)      
       write(6,*)
       do j=1,natom
        if(resnum(j)/=i)then
         write(6,1003)coord(j,1),coord(j,2),coord(j,3),charge(j)
         write(7,'(5x,3f14.8)')coord(j,1),coord(j,2),coord(j,3)
        endif
       enddo
1003  format(3x,3f14.8,f16.8)
       write(6,*)
       write(7,*)
       close(6)
       close(7)
      enddo
!=======================================================================
!=========================two body interaction==========================

      connect=.false.
      do i=1,nres-1
       do j=i+1,nres
        do ii=1,natom
        if((resnum(ii)==i).and.(atomname(ii)(2:2).eq.'O'))then
         do jj=1,natom
         if((resnum(jj)==j).and.(atomname(jj)(2:2).eq.'O'))then
       dist=sqrt((coord(ii,1)-coord(jj,1))**2+&
                 (coord(ii,2)-coord(jj,2))**2+&
                 (coord(ii,3)-coord(jj,3))**2)
          if(dist.le.10.0d0)then
           connect(i,j)=.true.
          endif
         endif
         enddo
        endif
        enddo
       enddo
      enddo

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      do i=mytid+1,nres,numprocs
       do j=i+1,nres
       if(connect(i,j))then
         open(6,file='000cell'//char(48+i/100)//&
                              char(48+(i-i/100*100)/10)//&
           char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                              char(48+(j-j/100*100)/10)//&
                        char(48+(j-j/10*10))//'.gjf',status='replace')
         open(7,file='field_000cell'//char(48+i/100)//&
                              char(48+(i-i/100*100)/10)//&
           char(48+(i-i/10*10))//'-000cell'//char(48+j/100)//&
                              char(48+(j-j/100*100)/10)//&
                        char(48+(j-j/10*10))//'.gjf',status='replace')
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
        write(6,'(a)')'#p CCD/aug-cc-pVDZ nosymm SCF=tight charge'//&
                  ' freq'
        write(7,'(a)')'#p density=(check,ccd) nosymm geom=allcheck '//&
                   'guess=(read,only) prop=(read,field)'
         write(6,*)
         write(6,'(a)')'Have a nice day'
         write(6,*)
         write(7,*)
         if(i /= 1)then
          write(6,'(a7)')'  0  1 '
         else
          write(6,'(a7)')'  1  1 '
         endif
         do l=1,natom
          if((resnum(l)==i).or.(resnum(l)==j))then
           write(6,1002)atomname(l)(2:2),(coord(l,k),k=1,3)
          endif
         enddo
         write(6,*)
         do l=1,natom
          if((resnum(l)/=i).and.(resnum(l)/=j))then
           write(6,1003)coord(l,1),coord(l,2),coord(l,3),charge(l)
           write(7,'(5x,3f14.8)')coord(l,1),coord(l,2),coord(l,3) 
          endif
         enddo
         write(6,*)
         write(7,*)
         close(6)
         close(7)
       endif
       enddo
      enddo
!=============================================================================

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

         end subroutine 
