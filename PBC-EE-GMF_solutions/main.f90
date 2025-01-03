           module variables
           implicit none
            
           integer(kind=4)::i,j,k,m,l,natom,resnum(9999),ist,&
                            indexx,ntype,ic,nxm,nym,nzm,ii,ni,&
                            jj,kk,nxl,nyl,nzl,jc,kc,ic1,jc1,kc1,&
                            njob,nres,nn,iter,nxll,nyll,nzll,nion,nwat,&
                          qmo(200,0:20),ncharge,group(200,20),mark,iacj,&
                          ionpair(200,2),ngroup
           real(kind=8),allocatable::charge(:),coord(:,:),force(:,:),&
             efragment(:,:,:,:),energy2(:,:,:,:,:),fele(:,:),fvdw(:,:),&
             cellcoord(:,:,:,:,:),tfield(:,:,:,:,:),field(:,:,:,:,:),&
             tforce(:,:),cn1(:),cn2(:),fvdw2(:,:),charge2(:),&
             txd(:),tyd(:),tzd(:),tdip(:),xd(:),yd(:),zd(:),dip(:),&
             tefragment(:,:,:,:),tenergy2(:,:,:,:,:),&
             pforce_global(:,:,:,:,:,:),&
             pfield_global(:,:,:,:,:,:),&
             pforce(:,:,:,:,:,:),&
             pfield(:,:,:,:,:,:)
           integer(kind=4),allocatable::iac(:),ico(:)
           real(kind=8)::dist,dist2,dinv6,dinv12,etotal,energy2body,&
                         temp(10),vdw_6,vdw_12,ele,evdw,ele2,evdw2,&
                         dinv8,dinv14,ele1,ebond,emm,clata,clatb,clatc,&
                         xa,ya,za,xc,yc,zc,cut,rqm,selfeng,da,db,dc,&
                         delta_a,delta_b,delta_c,gtol,fret,p,r2cut1,&
                         r2cut2,e2body,pi
           character(len=4)::atomname(9999)
           character(len=3)::resname(9999)
           character(len=80)::fline,cmdstr
           character(len=60)::cmd(5000)
           logical,allocatable::connect(:,:,:,:,:),onebody(:,:,:,:)
           logical::existence,normal,bgc(500),bgc2(500),scal(500)
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
           pi=3.141592653589793238d0

!====The periodic boundary condation control========
         nxm=1
         nym=1
         nzm=1
         nxl=2
         nyl=2
         nzl=2
         nxll=1
         nyll=1
         nzll=1
         p=1.0
!===================================================
      r2cut1=5.0d0
      r2cut2=5.0d0
      nion=0
      nwat=0
      qmo=0
      group=0
      open(106,file='H2O.pdb')
      do while(.true.)
       read(106,'(A80)',iostat=ist)fline
       if(ist/=0)exit
       if(fline(1:4).eq.'ATOM')then
        read(fline,'(6x,I5)')i
        read(fline,'(12x,a4)')atomname(i)
        if((atomname(i).eq.' Na+').or.(atomname(i).eq.' Cl-'))then
         nion=nion+1
        endif
        read(fline,'(17x,a3)')resname(i)
        if((resname(i).eq.'WAT').and.(atomname(i).eq.' O  '))then
         nwat=nwat+1
         qmo(nwat,0)=i
         group(nwat,1)=i
         group(nwat,2)=i+1
         group(nwat,3)=i+2
        endif
        read(fline,'(22x,I4)')resnum(i)
       endif
      enddo
      close(106)
      
      open(106,file='H2O.prmtop')                                                
      do while(.true.)
       read(106,'(A80)',iostat=ist)fline                                                  
       if(ist/=0)exit                                                                   
        if(fline(1:14).eq.'%FLAG POINTERS')then                                         
         read(106,*)                                                                      
         read(106,'(2I8)')natom,ntype                                                    
         nres=resnum(natom)
         allocate(charge(natom))
         allocate(charge2(natom))
         allocate(coord(natom,3))
         allocate(cellcoord(-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(onebody(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(connect(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(energy2(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(tenergy2(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(efragment(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(tefragment(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(force(natom,3))
         allocate(tforce(natom,3))
         allocate(field(-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(tfield(-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(pforce(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,10,3))
         allocate(pfield(nwat,-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(pforce_global(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,10,3))
        allocate(pfield_global(nwat,-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(fele(natom,3))
         allocate(fvdw(natom,3))
         allocate(fvdw2(natom,3))
         allocate(txd(nwat))
         allocate(tyd(nwat))
         allocate(tzd(nwat))
         allocate(tdip(nwat))
         allocate(xd(nwat))
         allocate(yd(nwat))
         allocate(zd(nwat))
         allocate(dip(nwat))
         allocate(iac(natom))
         allocate(ico(ntype*ntype))
         allocate(cn1(-1:ntype*(ntype+1)/2))
         allocate(cn2(-1:ntype*(ntype+1)/2))
         cn1=0.0d0
         cn2=0.0d0
        endif
        if(fline(1:12).eq.'%FLAG CHARGE')then                                           
         read(106,*)                                                                      
         read(106,1001)(charge(i),i=1,natom)                                              
        endif
        if(fline(2:21).eq.'FLAG ATOM_TYPE_INDEX') then
         read(106,*)
         read(106,1006)(iac(i),i=1,natom)
        endif
        if(fline(2:26).eq.'FLAG NONBONDED_PARM_INDEX') then
         read(106,*)
         read(106,1006)(ico(i),i=1,ntype*ntype)
        endif
        if(fline(2:25).eq.'FLAG LENNARD_JONES_ACOEF') then
         read(106,*)
         read(106,1001)(cn1(i),i=1,ntype*(ntype+1)/2)
        endif
        if(fline(2:25).eq.'FLAG LENNARD_JONES_BCOEF') then
         read(106,*)
         read(106,1001)(cn2(i),i=1,ntype*(ntype+1)/2)
        endif
      enddo
1001  format(5E16.8)                                                                    
1006  format(10I8)                                                                      
      close(106)
      
      charge=charge/18.2223

!     do i=1,natom
!      if((atomname(i).eq.' CS '))then
!       charge(i)=charge(i)*1.d0
!      endif
!     enddo

         call fragment_withchg
         call energy_cal

         deallocate(charge)
         deallocate(charge2)
         deallocate(coord)
         deallocate(cellcoord)
         deallocate(connect)
         deallocate(onebody)
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
         deallocate(fvdw)
         deallocate(fvdw2)
         deallocate(txd)
         deallocate(tyd)
         deallocate(tzd)
         deallocate(tdip)
         deallocate(xd)
         deallocate(yd)
         deallocate(zd)
         deallocate(dip)
         deallocate(iac)
         deallocate(ico)
         deallocate(cn1)
         deallocate(cn2)
         
        call MPI_FINALIZE(ierr)

         end program
