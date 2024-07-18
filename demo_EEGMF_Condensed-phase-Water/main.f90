           module variables
           implicit none
            
           integer(kind=4)::i,j,k,m,l,natom,resnum(9999),ist,&
                        nr,indexx,ntype,ic,nxm,nym,nzm,ii,ni,&
                         jj,kk,nxl,nyl,nzl,jc,kc,ic1,jc1,kc1,&
                  njob,nres,nn,iter,nxll,nyll,nzll,nion,nwat,&
             qmo(200,0:5),ncharge,group(200,7),mark,iacj,n3b,&
             ta0(99999),tb0(99999),tc0(99999),tcx(99999),&
             tcy(99999),tcz(99999),tbx(99999),tby(99999),&
             tbz(99999),mi,mj,dnx,dny,dnz,posx(27),posy(27),posz(27),&
             rbx(9999),rby(9999),rbz(9999),rb0(9999),&
             rcx(9999),rcy(9999),rcz(9999),rc0(9999),m1,m2
           real(kind=8),allocatable::charge(:),coord(:,:),force(:,:),&
             efragment(:,:,:,:),energy2(:,:,:,:,:),fele(:,:),fvdw(:,:),&
             cellcoord(:,:,:,:,:),tfield(:,:,:,:,:),field(:,:,:,:,:),&
             tforce(:,:),cn1(:),cn2(:),fvdw2(:,:),l_fele(:,:),&
             txd(:,:,:,:),tyd(:,:,:,:),tzd(:,:,:,:),tdip(:,:,:,:),&
             xd(:,:,:,:),yd(:,:,:,:),zd(:,:,:,:),dip(:,:,:,:),&
             tefragment(:,:,:,:),tenergy2(:,:,:,:,:),&
             pforce_global(:,:,:,:,:,:),&
             pfield_global(:,:,:,:,:,:),tx2d(:,:,:,:,:),&
             ty2d(:,:,:,:,:),tz2d(:,:,:,:,:),t2dip(:,:,:,:,:),&
             pforce(:,:,:,:,:,:),x2d(:,:),y2d(:,:),z2d(:,:),&
             pfield(:,:,:,:,:,:),l_fele2(:,:)
           integer(kind=4),allocatable::iac(:),ico(:)
           real(kind=8)::dist,dist2,dinv6,dinv12,etotal,energy2body,&
                         temp(30),vdw_6,vdw_12,ele,evdw,ele2,evdw2,&
                         dinv8,dinv14,ele1,ebond,emm,clata,clatb,clatc,&
                         xa,ya,za,xc,yc,zc,cut,rqm,selfeng,da,db,dc,&
                         delta_a,delta_b,delta_c,gtol,fret,p,r2cut1,&
                         r2cut2,e2body,rcut3,e3b,e3b2,e3b_qm,e3b2_qm,&
                         f2b(-1:1,-1:1,-1:1,128,128,6,3),&
                         f2_global(-1:1,-1:1,-1:1,128,128,6,3),&
                         ff2b(27,27,128,128,6,3),&
                         ff2_global(27,27,128,128,6,3),&
                        ee2b(27,27,128,128),ee2b_global(27,27,128,128),&
                        l_ele,l_ele2
           character(len=4)::atomname(9999)
           character(len=3)::resname(9999)
           character(len=80)::fline,cmdstr
           character(len=60)::cmd(99000)
           logical,allocatable::connect(:,:,:,:,:),onebody(:,:,:,:),&
                                connect3(:,:,:,:,:)
           logical::existence,normal,bgc(500),bgc2(500),bgc3(500),&
                    x2b(27,27,128,128)

           end module

           program fragment_optimization
           use variables
           use ewald
           implicit none 
           include 'parallel.h'
           include 'mpif.h' 


      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mytid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

           master=mytid.eq.0

!====The periodic boundary condation control========
         nxm = 1
         nym = 1
         nzm = 1
         nxl = 2
         nyl = 2
         nzl = 2
         nxll = 1
         nyll = 1
         nzll = 1
         p = 1.0
!===================================================
      r2cut1 = 7.0d0
      r2cut2 = 7.0d0
      rcut3 = 2.5d0
      nwat = 64
      natom = nwat * 3
      group = 0

         allocate(charge(natom))
         allocate(coord(natom,3))
         allocate(cellcoord(-nxl:nxl,-nyl:nyl,-nzl:nzl,natom,3))
         allocate(onebody(-nxl:nxl,-nyl:nyl,-nzl:nzl,nwat))
         allocate(connect(-nxl:nxl,-nyl:nyl,-nzl:nzl,nwat,nwat))
         allocate(connect3(-nxl:nxl,-nyl:nyl,-nzl:nzl,nwat,nwat))
         allocate(energy2(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(tenergy2(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(efragment(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(tefragment(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(force(natom,3))
         allocate(tforce(natom,3))
         allocate(field(-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(tfield(-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(pforce(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,7,3))
         allocate(pfield(nwat,-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(pforce_global(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,7,3))
        allocate(pfield_global(nwat,-nxm:nxm,-nym:nym,-nzm:nzm,natom,3))
         allocate(fele(natom,3))
         allocate(l_fele(natom,3))
         allocate(l_fele2(natom,3))
         allocate(fvdw(natom,3))
         allocate(fvdw2(natom,3))
         allocate(txd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(tyd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(tzd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(tx2d(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(ty2d(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(tz2d(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(t2dip(-nxm:nxm,-nym:nym,-nzm:nzm,nwat,nwat))
         allocate(tdip(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(xd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(yd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(zd(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))
         allocate(x2d(nwat,nwat))
         allocate(y2d(nwat,nwat))
         allocate(z2d(nwat,nwat))
         allocate(dip(-nxm:nxm,-nym:nym,-nzm:nzm,nwat))

      do i = 1, nwat
       group(i, 1) = i * 3 - 2
       group(i, 2) = i * 3 - 1
       group(i, 3) = i * 3 - 0
       atomname(i * 3 - 2) = ' O  '
       atomname(i * 3 - 1) = ' H1 '
       atomname(i * 3 - 0) = ' H2 '
       charge(i * 3 - 2) = -0.82d0
       charge(i * 3 - 1) =  0.41d0
       charge(i * 3 - 0) =  0.41d0
       resnum(i * 3 - 2) = i
       resnum(i * 3 - 1) = i
       resnum(i * 3 - 0) = i
      enddo
      
         call fragment_withchg
         call energy_cal

         deallocate(charge)
         deallocate(coord)
         deallocate(cellcoord)
         deallocate(connect)
         deallocate(connect3)
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
         deallocate(l_fele)
         deallocate(l_fele2)
         deallocate(fvdw)
         deallocate(fvdw2)
         deallocate(txd)
         deallocate(tyd)
         deallocate(tzd)
         deallocate(tx2d)
         deallocate(ty2d)
         deallocate(tz2d)
         deallocate(t2dip)
         deallocate(tdip)
         deallocate(xd)
         deallocate(yd)
         deallocate(zd)
         deallocate(x2d)
         deallocate(y2d)
         deallocate(z2d)
         deallocate(dip)
         
        call MPI_FINALIZE(ierr)

         end program
