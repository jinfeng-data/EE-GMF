           module variables
           implicit none
            
           integer(kind=4)::i,j,k,m,l,natom,resnum(9999),ist,&
                            indexx,ntype,ic,nxm,nym,nzm,ii,n,&
                            jj,kk,nxl,nyl,nzl,jc,kc,ic1,jc1,kc1,&
                            njob,nres,nn,iter,nxll,nyll,nzll
           real(kind=8),allocatable::charge(:),coord(:,:),force(:,:),&
             efragment(:),energy2(:,:),fele(:,:),&
             tfield(:,:),field(:,:),tforce(:,:),&
             tefragment(:),tenergy2(:,:),&
             pforce_global(:,:,:),&
             pfield_global(:,:,:),&
             pforce(:,:,:),&
             pfield(:,:,:)
           real(kind=8)::dist,dist2,dinv6,dinv12,etotal,energy2body,&
                         temp(999),vdw_6,vdw_12,ele,evdw,ele2,evdw2,&
                         dinv8,dinv14,ele1,ebond,emm,clata,clatb,clatc,&
                         xa,ya,za,xc,yc,zc,cut,rqm,selfeng,da,db,dc,&
                         delta_a,delta_b,delta_c,gtol,fret,p
           character(len=4)::atomname(9999)
           character(len=3)::resname(9999)
           character(len=80)::fline,cmdstr
           character(len=60)::cmd(5000)
           logical,allocatable::connect(:,:)
           logical::existence,normal
           end module

           program fragment_optimization
           use variables
           implicit none 
           include 'parallel.h'
           include 'mpif.h' 


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mytid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

           master=mytid.eq.0

      i=0
      open(15,file='cluster.pdb')
      do while(.true.)
       read(15,'(A80)',iostat=ist)fline
       if(ist/=0)exit
       if(fline(1:4).eq.'ATOM')then
        i=i+1
        read(fline,'(12x,a4)')atomname(i)
        read(fline,'(17x,a3)')resname(i)
        read(fline,'(22x,I4)')resnum(i)
       endif
      enddo
      close(15)
      
      natom=i
      nres=resnum(natom)
      allocate(charge(natom))
      allocate(coord(natom,3))
      allocate(connect(nres,nres))
      allocate(energy2(nres,nres))
      allocate(tenergy2(nres,nres))
      allocate(efragment(nres))
      allocate(tefragment(nres))
      allocate(force(natom,3))
      allocate(tforce(natom,3))
      allocate(field(natom,3))
      allocate(tfield(natom,3))
      allocate(pforce(nres,50,3))
      allocate(pfield(nres,natom,3))
      allocate(pforce_global(nres,50,3))
      allocate(pfield_global(nres,natom,3))
      allocate(fele(natom,3))
      
      do i=1,natom
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
      
         call fragment_withchg
         call energy_cal

         deallocate(charge)
         deallocate(coord)
         deallocate(connect)
         deallocate(energy2)
         deallocate(tenergy2)
         deallocate(efragment)
         deallocate(tefragment)
         deallocate(force)
         deallocate(tforce)
         deallocate(pforce)
         deallocate(pforce_global)
         deallocate(field)
         deallocate(tfield)
         deallocate(pfield)
         deallocate(pfield_global)
         deallocate(fele)
         
        call MPI_FINALIZE(ierr)

         end program
