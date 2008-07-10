!!
!! A _BoxArray_ is an array of boxes.
!!
module boxarray_module

  use bl_types
  use box_module
  use list_box_module
  use bl_mem_stat_module

  implicit none

  type boxarray
     integer :: dim = 0
     integer :: nboxes = 0
     type(box), pointer :: bxs(:) => Null()
  end type boxarray

  interface get_dim
     module procedure boxarray_dim
  end interface

  interface empty
     module procedure boxarray_empty
  end interface

  interface built_q
     module procedure boxarray_built_q
  end interface

  interface copy
     module procedure boxarray_build_copy
     module procedure boxarray_build_copy_l
  end interface

  interface build
     module procedure boxarray_build_v
     module procedure boxarray_build_l
     module procedure boxarray_build_bx
  end interface

  interface destroy
     module procedure boxarray_destroy
  end interface

  interface nboxes
     module procedure boxarray_nboxes
     module procedure boxlist_nboxes
  end interface

  interface volume
     module procedure boxarray_volume
  end interface

  interface dvolume
     module procedure boxarray_dvolume
  end interface

  interface set_box
     module procedure boxarray_set_box
  end interface

  interface get_box
     module procedure boxarray_get_box
  end interface

  interface boxarray_maxsize
     module procedure boxarray_maxsize_i
     module procedure boxarray_maxsize_v
  end interface

  interface boxarray_coarsen
     module procedure boxarray_coarsen_v
     module procedure boxarray_coarsen_i
  end interface

  interface boxarray_refine
     module procedure boxarray_refine_v
     module procedure boxarray_refine_i
  end interface

  interface boxarray_shift
     module procedure boxarray_shift_v
     module procedure boxarray_shift_i
  end interface

  interface boxarray_intersection
     module procedure boxarray_intersection_bx
  end interface

  interface boxarray_grow
     module procedure boxarray_grow_n
     module procedure boxarray_grow_n_f
     module procedure boxarray_grow_n_d_f
     module procedure boxarray_grow_v
     module procedure boxarray_grow_v_f
  end interface

  interface bbox
     module procedure boxarray_bbox
  end interface

  interface print
     module procedure boxarray_print
  end interface

  interface contains
     module procedure boxarray_box_contains
     module procedure boxarray_boxarray_contains
  end interface

  interface equal
     module procedure boxarray_equal
  end interface
  interface operator( .eq. )
     module procedure boxarray_equal
  end interface

  interface not_equal
     module procedure boxarray_not_equal
  end interface
  interface operator( .ne. )
     module procedure boxarray_not_equal
  end interface

  private :: boxlist_simplify
  private :: boxarray_maxsize_l
  private :: boxlist_box_diff
  private :: boxlist_build_a
  private :: boxlist_nboxes
  private :: boxlist_verify_dim

  type(mem_stats), private, save :: boxarray_ms

contains
  
  function boxarray_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = associated(ba1%bxs, ba2%bxs)
  end function boxarray_equal

  function boxarray_not_equal(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    r = .not. associated(ba1%bxs, ba2%bxs)
  end function boxarray_not_equal

  function boxarray_same_q(ba1, ba2) result(r)
    type(boxarray), intent(in) :: ba1, ba2
    logical :: r
    integer :: i
    if ( ba1 == ba2 ) then
       r = .true.
    else if ( ba1%dim /= ba2%dim .or. ba1%nboxes /= ba2%nboxes ) then
       r = .false.
    else
       do i = 1, ba1%nboxes
         if ( ba1%bxs(i) /= ba2%bxs(i) ) then
            r = .false.
            return
         end if
       end do
       r = .true.
    end if
  end function boxarray_same_q

  subroutine boxarray_set_mem_stats(ms)
    type(mem_stats), intent(in) :: ms
    boxarray_ms = ms
  end subroutine boxarray_set_mem_stats

  function boxarray_mem_stats() result(r)
    type(mem_stats) :: r
    r = boxarray_ms
  end function boxarray_mem_stats

  function boxarray_empty(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%nboxes == 0
  end function boxarray_empty

  function boxarray_built_q(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%dim /= 0
  end function boxarray_built_q

  function boxarray_dim(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%dim
  end function boxarray_dim

  subroutine boxarray_set_box(ba, i, bx)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: i
    type(box), intent(in) :: bx
    ba%bxs(i) = bx
  end subroutine boxarray_set_box

  function boxarray_get_box(ba, i) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: i
    type(box) :: r
    r = ba%bxs(i)
  end function boxarray_get_box

  subroutine boxarray_build_copy(ba, ba1)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(boxarray), intent(in) :: ba1
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_COPY: already built")
    if ( .not. built_q(ba1) ) return
    ba%nboxes = size(ba1%bxs)
    allocate(ba%bxs(size(ba1%bxs)))
    ba%bxs = ba1%bxs
    ba%dim = ba1%dim
    call boxarray_verify_dim(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_copy

  subroutine boxarray_build_copy_l(ba, bl)
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    if ( built_q(ba) ) call destroy(ba)
    call boxarray_build_l(ba, bl)
  end subroutine boxarray_build_copy_l

  subroutine boxarray_build_v(ba, bxs, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in), dimension(:) :: bxs
    logical, intent(in), optional :: sort
    logical :: lsort
    
    lsort = .false. ; if (present(sort)) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_V: already built")
    ba%nboxes = size(bxs)
    allocate(ba%bxs(size(bxs)))
    ba%bxs = bxs
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if (lsort) call boxarray_sort(ba) !! make sure all grids are sorted
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_v

  subroutine boxarray_build_bx(ba, bx)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_BX: already built")
    ba%nboxes = 1
    allocate(ba%bxs(1))
    ba%bxs(1) = bx
    ba%dim = bx%dim
    call boxarray_verify_dim(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_bx

  subroutine boxarray_build_l(ba, bl, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    logical, intent(in), optional :: sort
    type(list_box_node), pointer :: bln
    logical :: lsort
    integer :: i
    !
    ! Default is to sort.
    !
    lsort = .true. ; if ( present(sort) ) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_L: already built")
    ba%nboxes = size(bl)
    allocate(ba%bxs(ba%nboxes))
    bln => begin(bl)
    i = 1
    do while (associated(bln))
       ba%bxs(i) = value(bln)
       i = i + 1
       bln=>next(bln)
    end do
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if ( lsort ) call boxarray_sort(ba)
    call mem_stats_alloc(boxarray_ms, ba%nboxes)
  end subroutine boxarray_build_l

  subroutine boxarray_destroy(ba)
    type(boxarray), intent(inout) :: ba
    if ( associated(ba%bxs) ) then
       call mem_stats_dealloc(boxarray_ms, ba%nboxes)
       deallocate(ba%bxs) 
       ba%bxs => Null()
    end if
    ba%dim = 0
    ba%nboxes = 0
  end subroutine boxarray_destroy

  subroutine boxlist_build_a(bl, ba)
    type(boxarray), intent(in) :: ba
    type(list_box), intent(out) :: bl
    integer :: i
    do i = 1, ba%nboxes
       call push_back(bl, ba%bxs(i))
    end do
  end subroutine boxlist_build_a

  subroutine boxarray_sort(ba)
    use sort_box_module
    type(boxarray), intent(inout) :: ba
    call box_sort(ba%bxs)
  end subroutine boxarray_sort

  subroutine boxarray_verify_dim(ba, stat)
    use bl_error_module
    type(boxarray), intent(in) :: ba
    integer, intent(out), optional :: stat
    integer :: i, dm
    if ( present(stat) ) stat = 0
    if ( ba%nboxes < 1 ) return
    dm = ba%dim
    if ( dm == 0 ) then
       dm = ba%bxs(1)%dim
    end if
    if ( dm == 0 ) then
       call bl_error("BOXARRAY_VERIFY_DIM: dim is zero!")
    end if
    do i = 1, ba%nboxes
       if ( ba%dim /= ba%bxs(i)%dim ) then
          if ( present(stat) ) then
             stat = 1
             return
          else
             call bl_error("BOXARRAY_VERIFY_DIM: " // &
                  "ba%dim not equal to some boxes dim: ", ba%dim)
          end if
       end if
    end do
  end subroutine boxarray_verify_dim

  subroutine boxlist_verify_dim(bl, stat)
    use bl_error_module
    type(list_box), intent(in) :: bl
    integer, intent(out), optional :: stat
    type(list_box_node), pointer :: bln
    type(box) :: bx
    integer :: dm
    if ( present(stat) ) stat = 0
    if ( size(bl) < 1 ) return
    bln => begin(bl)
    bx = value(bln)
    dm = bx%dim
    do while (associated(bln))
       bx = value(bln)
       if ( bx%dim /= dm ) then
          if ( present(stat) ) then
             stat = 1
             return
          else
             call bl_error("BOXLIST_VERIFY_DIM:" // &
                  "some box's dim not equal to the first box's dim: ", dm)
          end if
       end if
    end do
  end subroutine boxlist_verify_dim

  subroutine boxarray_grow_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv)
    end do
  end subroutine boxarray_grow_v
  subroutine boxarray_grow_v_f(ba, rv, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:), face
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), rv, face)
    end do
  end subroutine boxarray_grow_v_f
  subroutine boxarray_grow_n(ba, n)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n)
    end do
  end subroutine boxarray_grow_n
  subroutine boxarray_grow_n_f(ba, n, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, face)
    end do
  end subroutine boxarray_grow_n_f
  subroutine boxarray_grow_n_d_f(ba, n, dim, face)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n, face, dim
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n, dim, face)
    end do
  end subroutine boxarray_grow_n_d_f

  subroutine boxarray_nodalize(ba, nodal)
    type(boxarray), intent(inout) :: ba
    logical, intent(in), optional :: nodal(:)
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = box_nodalize(ba%bxs(i), nodal)
    end do
  end subroutine boxarray_nodalize

  function boxarray_projectable(ba, rr) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: rr(:)
    integer :: i
    r = .true.
    do i = 1, nboxes(ba)
       if ( .not. box_projectable(ba%bxs(i), rr) ) then
          r = .false.
          exit
       end if
    end do
  end function boxarray_projectable

  subroutine boxarray_coarsen_v_m(ba, cv, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    logical, intent(in) :: mask(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv, mask)
    end do
  end subroutine boxarray_coarsen_v_m
  subroutine boxarray_coarsen_v(ba, cv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: cv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), cv)
    end do
  end subroutine boxarray_coarsen_v
  subroutine boxarray_coarsen_i_m(ba, ci, mask)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    logical, intent(in) :: mask(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci, mask)
    end do
  end subroutine boxarray_coarsen_i_m
  subroutine boxarray_coarsen_i(ba, ci)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ci
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = coarsen(ba%bxs(i), ci)
    end do
  end subroutine boxarray_coarsen_i

  subroutine boxarray_shift_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), rv)
    end do
  end subroutine boxarray_shift_v
  subroutine boxarray_shift_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = shift(ba%bxs(i), ri)
    end do
  end subroutine boxarray_shift_i
    
  subroutine boxarray_refine_v(ba, rv)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: rv(:)
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), rv)
    end do
  end subroutine boxarray_refine_v
  subroutine boxarray_refine_i(ba, ri)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: ri
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = refine(ba%bxs(i), ri)
    end do
  end subroutine boxarray_refine_i
  !
  ! This is a very naive implementation.
  !
  subroutine boxarray_intersection_bx(ba, bx)
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    integer :: i
    do i = 1, ba%nboxes
      ba%bxs(i) = intersection(ba%bxs(i), bx)
    end do
    call boxarray_simplify(ba)
  end subroutine boxarray_intersection_bx

  subroutine boxarray_box_boundary_v_f(bao, bx, nv, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: nv(:), face
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_v_f(bao, baa, nv, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_v_f
  subroutine boxarray_box_boundary_v(bao, bx, nv)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: nv(:)
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_v(bao, baa, nv)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_v
  subroutine boxarray_box_boundary_n_f(bao, bx, n, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n, face
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n_f(bao, baa, n, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n_f
  subroutine boxarray_box_boundary_n_d_f(bao, bx, n, dim, face)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n, face, dim
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n_d_f(bao, baa, n, dim, face)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n_d_f
  subroutine boxarray_box_boundary_n(bao, bx, n)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n(bao, baa, n)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n

  subroutine boxarray_boundary_v_f(bao, ba, nv, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: nv(:), face
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, nv, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_v_f
  subroutine boxarray_boundary_v(bao, ba, nv)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: nv(:)
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, nv)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_v
  subroutine boxarray_boundary_n_d_f(bao, ba, n, dim, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n, face, dim
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n, dim, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n_d_f
  subroutine boxarray_boundary_n_f(bao, ba, n, face)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer, intent(in) :: n, face
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n, face)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n_f
  subroutine boxarray_boundary_n(bao, ba, n)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer,        intent(in)  :: n
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n
  function boxarray_nboxes(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%nboxes
  end function boxarray_nboxes

  function boxlist_nboxes(bl) result(r)
    type(list_box), intent(in) :: bl
    integer :: r
    r = size(bl)
  end function boxlist_nboxes

  function boxarray_volume(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer(kind=ll_t) :: r
    integer :: i
    r = 0_ll_t
    do i = 1, ba%nboxes
       r = r + box_volume(ba%bxs(i))
    end do
  end function boxarray_volume

  function boxarray_dvolume(ba) result(r)
    type(boxarray), intent(in) :: ba
    real(dp_t) :: r
    integer :: i
    r = 0
    do i = 1, ba%nboxes
       r = r + box_dvolume(ba%bxs(i))
    end do
  end function boxarray_dvolume

  function boxarray_bbox(ba) result(r)
    type(boxarray), intent(in) :: ba
    type(box) :: r
    integer :: i
    r = nobox(ba%dim)
    do i = 1, ba%nboxes
       r = bbox(r, ba%bxs(i))
    end do
  end function boxarray_bbox

  subroutine boxarray_box_diff(ba, b1, b2)
    type(boxarray), intent(out) :: ba
    type(box), intent(in) :: b1, b2
    type(list_box) :: bl
    bl = boxlist_box_diff(b1, b2)
    call boxarray_build_l(ba, bl)
    call destroy(bl)
  end subroutine boxarray_box_diff

  ! a retro name
  subroutine boxarray_complementIn(ba, bx, ba1)
    type(boxarray), intent(out) :: ba
    type(boxarray), intent(in) :: ba1
    type(box),intent(in) :: bx
    call boxarray_boxarray_diff(ba, bx, ba1)
  end subroutine boxarray_complementIn

  subroutine boxarray_boxarray_diff(ba, bx, ba1)
    type(boxarray), intent(out) :: ba
    type(boxarray), intent(in) :: ba1
    type(box), intent(in) :: bx
    type(list_box) :: bl1, bl
    type(box) :: bx1
    integer :: i
    call build(bl1)
    do i = 1, nboxes(ba1)
       bx1 = intersection(bx, get_box(ba1,i))
       if ( empty(bx1) ) cycle
       call push_back(bl1, bx1)
    end do
    bl = boxlist_boxlist_diff(bx, bl1)
    call boxarray_build_l(ba, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine boxarray_boxarray_diff

  subroutine boxarray_diff(bao, ba)
    type(boxarray), intent(inout) :: bao
    type(boxarray), intent(in) :: ba
    type(list_box) :: bl, bl1, bl2
    integer :: i
    call build(bl1, ba%bxs)
    do i = 1, bao%nboxes
       bl2 = boxlist_boxlist_diff(bao%bxs(i), bl1)
       call splice(bl, bl2)
    end do
    call boxarray_destroy(bao)
    call boxarray_build_l(bao, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine boxarray_diff

  subroutine boxarray_maxsize_i(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in) :: chunk
    integer :: vchunk(bxa%dim) 
    vchunk = chunk
    call boxarray_maxsize_v(bxa, vchunk)
  end subroutine boxarray_maxsize_i

  subroutine boxarray_maxsize_v(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in), dimension(:) :: chunk
    type(list_box) :: bl
    bl = boxarray_maxsize_l(bxa, chunk)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bl)
    call destroy(bl)
  end subroutine boxarray_maxsize_v

  function boxarray_maxsize_l(bxa, chunk) result(r)
    type(list_box) :: r
    type(boxarray), intent(in) ::  bxa
    integer, intent(in), dimension(:) :: chunk
    integer :: i,k
    type(list_box_node), pointer :: li
    integer :: len(bxa%dim)
    integer :: nl, bs, rt, nblk, sz, ex, ks, ps
    type(box) :: bxr, bxl

    do i = 1, bxa%nboxes
       call push_back(r, bxa%bxs(i))
    end do
    li => begin(r)
    do while ( associated(li) )
       len = extent(value(li))
       do i = 1, bxa%dim
          if ( len(i) > chunk(i) ) then
             rt = 1
             bs = chunk(i)
             nl = len(i)
             do while ( mod(bs,2) == 0 .AND. mod(nl,2) == 0)
                rt = rt * 2
                bs = bs/2
                nl = nl/2
             end do
             nblk = nl/bs
             if ( mod(nl,bs) /= 0 ) nblk = nblk + 1
             sz   = nl/nblk
             ex   = mod(nl,nblk)
             do k = 0, nblk-2
                if ( k < ex ) then
                   ks = (sz+1)*rt
                else
                   ks = sz*rt
                end if
                ps = upb(value(li), i) - ks + 1
                call box_chop(value(li), bxr, bxl, i, ps)
                call set(li, bxr)
                call push_back(r, bxl)
             end do
          end if
       end do
       li => next(li)
    end do

  end function boxarray_maxsize_l

  function boxlist_box_diff(bx1, b2) result(r)
    type(box), intent(in) :: bx1, b2
    type(list_box) :: r
    type(box) :: b1, bn
    integer, dimension(bx1%dim) :: b2lo, b2hi, b1lo, b1hi
    integer :: i, dm

    dm = bx1%dim
    b1 = bx1
    if ( .not. contains(b2,b1) ) then
       if ( .not. intersects(b1,b2) ) then
          call push_back(r, b1)
       else
          b2lo = lwb(b2); b2hi = upb(b2)
          do i = 1, dm
             b1lo = lwb(b1); b1hi = upb(b1)
             if ( b1lo(i) < b2lo(i) .AND. b2lo(i) <= b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b1lo(i))
                call set_upb(bn, i, b2lo(i)-1)
                call push_back(r, bn)
                call set_lwb(b1, i, b2lo(i))
             end if
             if ( b1lo(i) <= b2hi(i) .AND. b2hi(i) < b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b2hi(i)+1)
                call set_upb(bn, i, b1hi(i))
                call push_back(r, bn)
                call set_upb(b1, i, b2hi(i))
             end if
          end do
       end if
    end if

  end function boxlist_box_diff

  ! r = bx - bxl
  function boxlist_boxlist_diff(bx, bxl) result(r)

    use bl_prof_module

    type(box), intent(in) :: bx
    type(list_box), intent(in) :: bxl
    type(list_box) :: r, bl
    type(list_box_node), pointer :: blp, bp
    type(bl_prof_timer), save :: bpt

    call build(bpt, "boxlist_boxlist_diff")

    call push_back(r, bx)
    blp => begin(bxl)
    do while ( associated(blp) .AND. .NOT. empty(r) )
       bp => begin(r)
       do while ( associated(bp) )
          if ( intersects(value(bp), value(blp)) ) then
             bl = boxlist_box_diff(value(bp), value(blp))
             call splice(r, bl)
             bp => erase(r, bp)
          else
             bp => next(bp)
          end if
       end do
       blp => next(blp)
    end do

    call destroy(bpt)

  end function boxlist_boxlist_diff

  subroutine boxarray_simplify(bxa)
    type(boxarray), intent(inout) :: bxa
    type(list_box) :: bxl
    call build(bxl, bxa%bxs)
    call boxlist_simplify(bxl)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bxl)
    call destroy(bxl)
  end subroutine boxarray_simplify

  subroutine boxlist_simplify(bxl)

    use bl_prof_module

    type(list_box), intent(inout) :: bxl
    integer :: dm
    type(bl_prof_timer), save :: bpt

    if ( size(bxl) == 0 ) return

    call build(bpt, "boxlist_simplify")

    dm = box_dim(front(bxl))
    !
    ! TODO -- limit number of simp() calls to 1, 2 or 3?
    !
    do while ( simp() > 0 )
    end do

    call destroy(bpt)

  contains

    function simp() result(cnt)
      use box_module
      integer :: cnt
      integer :: joincnt
      integer, dimension(dm) :: lo, hi, alo, ahi, blo, bhi
      type(list_box_node), pointer :: ba, bb
      type(box) :: bx
      logical match, canjoin
      integer :: i

      if ( size(bxl) == 0 ) return
      cnt = 0
      ba => begin(bxl)

      do while ( associated(ba) )
         alo = lwb(value(ba)); ahi = upb(value(ba))
         match = .FALSE.
         bb => next(ba)
         do while ( associated(bb) )
            blo = lwb(value(bb)); bhi = upb(value(bb))
            canjoin = .TRUE.
            joincnt = 0
            do i = 1, dm
               if ( alo(i) == blo(i) .AND. ahi(i)==bhi(i) ) then
                  lo(i) = alo(i)
                  hi(i) = ahi(i)
               else if ( alo(i)<=blo(i) .AND. blo(i)<=ahi(i)+1 ) then
                  lo(i) = alo(i)
                  hi(i) = max(ahi(i),bhi(i))
                  joincnt = joincnt + 1
               else if ( blo(i)<=alo(i) .AND. alo(i)<=bhi(i)+1 ) then
                  lo(i) = blo(i)
                  hi(i) = max(ahi(i),bhi(i))
                  joincnt = joincnt + 1
               else
                  canjoin = .FALSE.
                  exit
               end if
            end do
            if ( canjoin .AND. (joincnt <= 1) ) then
               ! Modify b and remove a from the list.
               call build(bx, lo, hi)
               call set(bb, bx)
               ba => erase(bxl, ba)
               cnt = cnt + 1
               match = .TRUE.
               exit
            else
               ! No match found, try next element.
               bb => next(bb)
            end if
         end do
         ! If a match was found, a was already advanced in the list.
         if ( .not. match ) then
            ba => next(ba)
         end if
      end do

    end function simp

  end subroutine boxlist_simplify

  subroutine boxarray_print(ba, str, unit, legacy, skip)
    use bl_IO_module
    type(boxarray), intent(in) :: ba
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: legacy
    integer, intent(in), optional :: skip
    integer :: i
    integer :: un
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("BOXARRAY[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i5)') ba%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i5)') ba%nboxes
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    do i = 1, ba%nboxes
       call print(ba%bxs(i), unit=unit, advance = 'NO', &
            legacy = legacy, skip = unit_get_skip(skip)+ 1)
       if ( i == ba%nboxes ) then
          write(unit=un, fmt='("}]")')
       else
          write(unit=un, fmt='(",")')
       end if
    end do
  end subroutine boxarray_print

  subroutine boxlist_print(bl, str, unit, legacy, skip)
    use bl_IO_module
    type(list_box), intent(in) :: bl
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: legacy
    integer, intent(in), optional :: skip
    type(list_box_node), pointer :: bn
    integer :: un, cnt
    un = unit_stdout(unit)
    call unit_skip(un, skip)
    write(unit=un, fmt = '("BOXLIST[(*")', advance = 'no')
    if ( present(str) ) then
       write(unit=un, fmt='(" ",A)') str
    else
       write(unit=un, fmt='()')
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NBOXES  = ",i5)') size(bl)
    call unit_skip(un, skip)
    write(unit=un, fmt='(" *) {")')
    cnt = 0
    bn => begin(bl)
    do while ( associated(bn) )
       call print(value(bn), unit=unit, advance = 'NO',  &
            legacy = legacy, skip = unit_get_skip(skip)+ 1)
       bn => next(bn)
       cnt = cnt + 1
       if ( cnt == size(bl) ) then
          write(unit=un, fmt='(",")')
       else
          write(unit=un, fmt='(",")')
       end if
    end do
  end subroutine boxlist_print

  function boxarray_clean(boxes) result(r)
    logical :: r
    type(box), intent(in), dimension(:) :: boxes
    integer :: i, j

    do i = 1, size(boxes)-1
       do j = i+1, size(boxes)
          if ( intersects(boxes(i),boxes(j)) ) then
             r = .FALSE.
             return
          end if
       end do
    end do
    r = .TRUE.

  end function boxarray_clean

  subroutine boxarray_add_clean(ba, bx)

    use bl_prof_module
    use list_box_module

    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp
    type(bl_prof_timer), save :: bpt

    if ( empty(ba) ) then
       call boxarray_build_bx(ba, bx)
       return
    end if
    call build(bpt, "ba_add_clean")
    call build(bl, ba%bxs)
    call push_back(check, bx)
    lp => begin(bl)
    do while ( associated(lp) )
       cp => begin(check)
       do while ( associated(cp) )
          if ( intersects(value(cp), value(lp)) ) then
             tmpbl = boxlist_box_diff(value(cp), value(lp))
             call splice(tmp, tmpbl)
             cp => erase(check, cp)
          else
             cp => next(cp)
          end if
       end do
       call splice(check, tmp)
       lp => next(lp)
    end do
    call splice(bl, check)
    call boxlist_simplify(bl)
    call boxarray_build_copy_l(ba, bl)
    call destroy(bl)
    call destroy(bpt)

  end subroutine boxarray_add_clean

  subroutine boxarray_add_clean_boxes(ba, bxs, simplify)

    use bl_prof_module
    use list_box_module

    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bxs(:)
    logical, intent(in), optional :: simplify
    logical :: lsimplify
    type(list_box) :: check, tmp, tmpbl, bl
    type(list_box_node), pointer :: cp, lp
    integer :: i
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ba_add_clean_boxes")

    lsimplify = .true.; if ( present(simplify) ) lsimplify = simplify

    if ( size(bxs) .eq. 0 ) return
    
    if ( empty(ba) ) call boxarray_build_bx(ba, bxs(1))

    call build(bl, ba%bxs)

    do i = 1, size(bxs)
       call build(check, bxs(i:i))
       lp => begin(bl)
       do while ( associated(lp) )
          cp => begin(check)
          do while ( associated(cp) )
             if ( intersects(value(cp), value(lp)) ) then
                tmpbl = boxlist_box_diff(value(cp), value(lp))
                call splice(tmp, tmpbl)
                cp => erase(check, cp)
             else
                cp => next(cp)
             end if
          end do
          call splice(check, tmp)
          lp => next(lp)
       end do
       call splice(bl, check)
    end do

    if ( lsimplify ) call boxlist_simplify(bl)

    call boxarray_build_copy_l(ba, bl)
    call destroy(bl)
    call destroy(bpt)

  end subroutine boxarray_add_clean_boxes

  subroutine boxarray_to_domain(ba)
    type(boxarray), intent(inout) :: ba
    type(boxarray) :: ba1
    call boxarray_add_clean_boxes(ba1, ba%bxs)
    call boxarray_destroy(ba)
    ba = ba1
  end subroutine boxarray_to_domain

  subroutine boxarray_box_corners(ba, bx, ng)
    use bl_error_module

    type(boxarray), intent(out) :: ba
    type(box),      intent(in)  :: bx
    integer,        intent(in)  :: ng
    integer                     :: i, len(1:bx%dim)
    type(boxarray)              :: tba

    if (ng < 0) call bl_error("BOXARRAY_BOX_CORNERS: ng must be >= 0!")

    call boxarray_build_bx(ba, grow(bx, ng))

    len = 0
    do i = 1, bx%dim
       len(i) = ng
       call boxarray_build_bx(tba, grow(bx, len))
       call boxarray_diff(ba, tba)
       call destroy(tba)
       len(i) = 0
    end do
  end subroutine boxarray_box_corners

  function boxarray_box_contains(ba, bx) result(r)
    use bl_error_module
    logical                    :: r
    type(boxarray), intent(in) :: ba
    type(box),      intent(in) :: bx

    type(list_box) :: bl1, bl
    type(box)      :: bx1
    integer        :: i

    if ( nboxes(ba) .eq. 0 ) &
       call bl_error('Empty boxarray in boxarray_box_contains')

    call build(bl1)
    do i = 1, nboxes(ba)
       bx1 = intersection(bx, get_box(ba,i))
       if ( empty(bx1) ) cycle
       call push_back(bl1, bx1)
    end do
    bl = boxlist_boxlist_diff(bx, bl1)
    r = empty(bl)
    call destroy(bl)
    call destroy(bl1)

  end function boxarray_box_contains

  function boxarray_boxarray_contains(ba1, ba2) result(r)
    use bl_error_module
    logical :: r
    type(boxarray), intent(in) :: ba1, ba2

    integer :: i

    if ( nboxes(ba1) .eq. 0 ) &
       call bl_error('Empty boxarray ba1 in boxarray_boxarray_contains')

    if ( nboxes(ba2) .eq. 0 ) &
       call bl_error('Empty boxarray ba2 in boxarray_boxarray_contains')

    do i = 1, nboxes(ba2)
       r = boxarray_box_contains(ba1, get_box(ba2,i)) 
       if ( .not. r ) return
    end do

    r = .true.

  end function boxarray_boxarray_contains

end module boxarray_module
