      integer(kind=4)::mytid,numprocs,ierr
      logical master
      common /parallel/mytid,numprocs,ierr,master
