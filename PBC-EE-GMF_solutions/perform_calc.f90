      program main 
      implicit none

      integer(kind=4)::i,j,k,natom,ist,n
      real(kind=8)::box(3),coord(1000,3)

      natom=427
      n=0

      box(1)=16.3961921d0
      box(2)=16.2746506d0
      box(3)=16.2495958d0
      open(123,file='trj_2.crd')
       do while(.true.)
        read(123,'(6f12.7)',iostat=ist)((coord(i,j),j=1,3),i=1,natom)
        if(ist/=0)exit
        read(123,*)
        read(123,*)
        n=n+1
        if(n > 0)then
        open(126,file='H2O.inpcrd',status='replace')
         write(126,'(a)')'Liquid Water Dynamics'
         write(126,'(I6)')natom
         write(126,'(6f12.7)')((coord(i,j),j=1,3),i=1,natom)
         write(126,'(6f12.7)')box(1:3)
        close(126)
        call system("./liquid_water_dynamic.x")
        endif
       enddo
      close(123)

      end program
