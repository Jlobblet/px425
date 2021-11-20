!==========================================================!
!  Cluster search algorithms for use with PX425 (2017)     !
!  assignment 3.                                           !
!                                                          !
!  Original code created by D.Quigley - November 2014      !
!==========================================================!
module cluster_search

  use iso_c_binding
  implicit none       ! Impose strong typing

  private             ! Everything private by default...

  ! kind parameters 
  integer,parameter :: dp = c_double
  integer,parameter :: it = c_long

  !... unless exposed here
  public :: find_clusters_recursive
  public :: find_clusters_eqclass

contains

  subroutine find_clusters_recursive(Nvert,Maxcon,Ncon,Lcon,lclus,nclus) bind(c)
    !-------------------------------------------------------------------!
    ! Subroutine to build a list of connected vertex clusters, and      !
    ! report the size of the largest cluster. Uses a recursive search   !
    ! over connections to each vertex.                                  !
    !-------------------------------------------------------------------!
    ! D. Quigley - University of Warwick                                !
    !===================================================================!
    implicit none

    integer(kind=it),value,intent(in) :: Nvert       ! Number of vertices
    integer(kind=it),value,intent(in) :: Maxcon      ! Max connections per vertex

    ! Number of edges per vertex - C style indexing
    integer(kind=it),dimension(0:Nvert-1) :: Ncon

    ! List of vertices joined to each vertex by edges
    integer(kind=it),dimension(0:Maxcon-1,0:Nvert-1),intent(in) :: Lcon

    ! Outputs - largest cluster and number of clusters
    integer(kind=it),intent(out) :: lclus,nclus

    ! Logical array indicating which vertices have already been counted
    logical,allocatable,dimension(:) :: lvisited

    ! Array holding size of each cluster
    integer(kind=it),allocatable,dimension(:) :: cluster_size

    ! Loop counters and error flags
    integer(kind=it) :: iv,ierr,ivcluster

    !-------------------------!
    ! Allocate memory         !
    !-------------------------!
    allocate(lvisited(0:Nvert-1),stat=ierr)
    if (ierr/=0) stop 'Error allocating lvisited in find_clusters_recursive'
    allocate(cluster_size(0:Nvert-1),stat=ierr)
    if (ierr/=0) stop 'Error allocating cluster_size in find_clusters_recursive'

    !-------------------------!
    ! Perform analysis        !
    !-------------------------!

    ! Initialise cluster count
    ivcluster = 0
    lvisited  = .false.
    cluster_size = 0

    
    ! loop over vertices
    do iv = 0,Nvert-1

       ! if we haven't already visited this vertex
       if (.not.lvisited(iv)) then

          ! .. then we have now ..
          lvisited(iv) = .true.

          ! New cluster at current ivcluster
          cluster_size(ivcluster)   = 0
          !cluster_list(:,ivcluster) = 0

          ! now add vertex iv to that cluster
          cluster_size(ivcluster) = cluster_size(ivcluster) + 1

          ! next search onward over all edges involving this vertex
          ! on return from this call, we should have followed
          ! all possible links originating on vertex iv, i.e. found
          ! all members of the cluster containing iv.
          call vertex_search(iv,ivcluster)

          ! Next cluster will be...
          ivcluster = ivcluster + 1

       end if

    end do

    ! Analysis complete. Report number of clusters and size of cluster
    lclus = maxval(cluster_size)
    nclus = ivcluster

    !-------------------------!
    ! Release memory          !
    !-------------------------!
    deallocate(lvisited,cluster_size,stat=ierr) ! cluster_list
    if (ierr/=0) stop 'Error releasing memory in find_clusters_recursive'


  contains

    recursive subroutine vertex_search(i,icluster)
      !-------------------------------------------------------------------!
      ! Recursive subroutine to identify all vertices in the same cluster !
      ! as vertix i, which resides in cluster number icluster.            !
      !-------------------------------------------------------------------!
      ! D. Quigley - University of Warwick                                !
      !===================================================================!
      implicit none
      integer(kind=it),intent(in) :: i,icluster
      integer(kind=it)            :: j,jlist

      if (Ncon(i)==0) return ! exit if no connections

      do jlist = 0,Ncon(i)-1

         j = Lcon(jlist,i)

         if (.not.lvisited(j) ) then

            lvisited(j) = .true.      ! j has been reached

            ! Increment the size of icluster, and
            ! add j to the list of vertices in icluster
            cluster_size(icluster) = cluster_size(icluster) + 1

            ! Search onward from here, calling this routine recursively
            call vertex_search(j,icluster)

         end if

      end do

      return

    end subroutine  vertex_search

  end subroutine find_clusters_recursive

  subroutine find_clusters_eqclass(Nvert,Maxcon,Ncon,Lcon,lclus,nclus) bind(c)
    !-------------------------------------------------------------------!
    ! Subroutine to build a list of connected vertex clusters, and      !
    ! report the size of the largest cluster. Identifies connected      !
    ! vertices as equivalent and assigns each to an equivalence class   !
    ! See Numerical Recipes by Press et al for details.                 !
    !-------------------------------------------------------------------!
    ! D. Quigley - University of Warwick                                !
    !===================================================================!
    implicit none

    integer(kind=it),value,intent(in) :: Nvert       ! Number of vertices
    integer(kind=it),value,intent(in) :: Maxcon      ! Max connections per vertex

    ! Number of edges per vertex
    integer(kind=it),dimension(0:Nvert-1) :: Ncon

    ! List of vertices jointed to each vertex by edges
    integer(kind=it),dimension(0:Maxcon-1,0:Nvert-1),intent(in) :: Lcon

    ! Outputs - largest cluster and number of clusters
    integer(kind=it),intent(out) :: lclus,nclus

    integer(kind=it),allocatable,dimension(:) :: cluster_size

    integer(kind=it),allocatable,dimension(:) :: lcl

    integer(kind=it) :: nvcluster

    integer(kind=it) :: ic,iv,jv,ierr


    allocate(lcl(0:Nvert-1),stat=ierr)
    if (ierr/=0) stop 'Error allocating lcl in find_clusters_eqclass'
   
    allocate(cluster_size(0:Nvert-1),stat=ierr)
    if (ierr/=0) stop 'Error allocating cluster_size in find_clusters_recursive'

    lcl = -1

    do iv = 0,nvert-1
       lcl(iv)   = iv                           
       do jv = 0,iv-2                               
          lcl(jv) = lcl(lcl(jv))                               
          do ic = 0,Ncon(iv)-1
             if (Lcon(ic,iv)==jv) then  ! if iv and jv equivalent
                lcl(lcl(lcl(jv))) = iv
             end if
          end do
       end do
    end do

    cluster_size = 0
    do iv = 0,nvert-1
       if (lcl(iv) /= -1) then
          lcl(iv) = lcl(lcl(iv))
          cluster_size(lcl(iv)) = cluster_size(lcl(iv)) + 1
       end if
    end do

    nvcluster    = 0
    do iv = 0,nvert-1
       if (cluster_size(iv)>0) nvcluster = nvcluster + 1
    end do

    nclus = nvcluster 
    lclus = maxval(cluster_size)

    deallocate(lcl,cluster_size,stat=ierr)  
    if (ierr/=0) stop 'Error releasing memory in find_clusters_eqclass'

  end subroutine find_clusters_eqclass

end module cluster_search
