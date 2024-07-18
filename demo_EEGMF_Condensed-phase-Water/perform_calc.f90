      program main 
      implicit none

      integer(kind=4)::i,j,k,natom
      real(kind=8)::box(3),coord(1000,3)

      natom=64 * 3

!     box(1)=9.0711937d0
!     box(2)=15.6308467d0
!     box(3)=14.7152310d0
      open(123,file='struc.crd')
       do k=1,1000
        write(*,*)k
        read(123,'(6f12.7)')((coord(i,j),j=1,3),i=1,natom)
        read(123,'(3f12.7)')box(1:3)
        read(123,*)
        if(k > 0)then
        open(126,file='H2O.inpcrd',status='replace')
         write(126,'(a)')'Liquid Water Dynamics'
         write(126,'(I6)')natom
         write(126,'(6f12.7)')((coord(i,j),j=1,3),i=1,natom)
         write(126,'(6f12.7)')box(1:3)
        close(126)
        call system("./run_frag.sh")
        endif
       enddo
      close(123)

      end program
