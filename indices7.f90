! Corrected spin-1 program. Build with:
! gfortran -O2 -ffree-line-length-none indices_corrected.f90 -llapack -lblas -o indices_corrected
! This source has not been compiler-tested in the audit environment.
! indices_verified.py independently certifies the SAME exact input tensors.
! Fixed inputs: NumPy default_rng(4), then preserved explicitly below.
! Normality means the span of products of ONE fixed length is M_D.
module cg_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains

  real(dp) function fact(n)
    integer, intent(in) :: n
    integer :: i
    fact = 1.0_dp
    if (n <= 1) return
    do i = 2, n
      fact = fact * real(i, dp)
    end do
  end function fact

  real(dp) function cg(j1, m1, j2, m2, Jtot, Mtot)
    real(dp), intent(in) :: j1, m1, j2, m2, Jtot, Mtot
    real(dp) :: delta, prefac, sumterm, term
    integer :: kmin, kmax, k
    integer :: a1, a2, a3, a4, a5, a6, a7, a8, a9

    cg = 0.0_dp
    if (abs(Mtot - (m1+m2)) > 1.0d-9) return
    if (Jtot > j1+j2+1.0d-9 .or. Jtot < abs(j1-j2)-1.0d-9) return

    a1 = nint(j1+j2-Jtot); a2 = nint(j1-j2+Jtot); a3 = nint(-j1+j2+Jtot)
    a4 = nint(j1+j2+Jtot+1.0_dp)
    if (a1 < 0 .or. a2 < 0 .or. a3 < 0) return

    delta = sqrt(fact(a1)*fact(a2)*fact(a3) / fact(a4))

    a5 = nint(Jtot+Mtot); a6 = nint(Jtot-Mtot)
    a7 = nint(j1-m1);     a8 = nint(j1+m1);    a9 = nint(j2-m2)
    if (a5 < 0 .or. a6 < 0 .or. a7 < 0 .or. a8 < 0 .or. a9 < 0) return

    prefac = delta * sqrt((2.0_dp*Jtot+1.0_dp) * fact(a5)*fact(a6)* &
             fact(a7)*fact(a8)*fact(a9)*fact(nint(j2+m2)))

    kmin = max(0, nint(j2-Jtot-m1), nint(j1-Jtot+m2))
    kmax = min(nint(j1+j2-Jtot), nint(j1-m1), nint(j2+m2))

    sumterm = 0.0_dp
    do k = kmin, kmax
      if ((nint(j1+j2-Jtot)-k) < 0) cycle
      if ((nint(j1-m1)-k) < 0) cycle
      if ((nint(j2+m2)-k) < 0) cycle
      if ((nint(Jtot-j2+m1)+k) < 0) cycle
      if ((nint(Jtot-j1-m2)+k) < 0) cycle

      term = 1.0_dp / (fact(k) * fact(nint(j1+j2-Jtot)-k) * &
             fact(nint(j1-m1)-k) * fact(nint(j2+m2)-k) * &
             fact(nint(Jtot-j2+m1)+k) * fact(nint(Jtot-j1-m2)+k))

      if (mod(k,2) == 0) then
        sumterm = sumterm + term
      else
        sumterm = sumterm - term
      end if
    end do

    cg = prefac * sumterm
  end function cg

  subroutine cg_matrix(jout, jin, q, s_phys, C, dout, din)
    real(dp), intent(in) :: jout, jin, s_phys
    integer, intent(in) :: q, dout, din
    complex(dp), intent(out) :: C(dout,din)
    integer :: a, b
    real(dp) :: mo, mi

    C = (0.0_dp, 0.0_dp)
    do a = 1, dout
      mo = jout - real(a-1, dp)
      do b = 1, din
        mi = jin - real(b-1, dp)
        C(a,b) = cmplx(cg(jin, mi, s_phys, real(q,dp), jout, mo), 0.0_dp, dp)
      end do
    end do
  end subroutine cg_matrix

end module cg_mod
module verified_mps
  use cg_mod
  implicit none
contains
  subroutine kron(A, na1, na2, B, nb1, nb2, C)
    integer, intent(in) :: na1, na2, nb1, nb2
    complex(dp), intent(in) :: A(na1,na2), B(nb1,nb2)
    complex(dp), intent(out) :: C(na1*nb1, na2*nb2)
    integer :: i, j, p, q

    do i = 1, na1
      do j = 1, na2
        do p = 1, nb1
          do q = 1, nb2
            C((i-1)*nb1+p, (j-1)*nb2+q) = A(i,j)*B(p,q)
          end do
        end do
      end do
    end do
  end subroutine kron

  subroutine time_reversal_matrix(j, Dmat, dim_d)
    real(dp), intent(in) :: j
    integer, intent(in) :: dim_d
    complex(dp), intent(out) :: Dmat(dim_d,dim_d)
    integer :: a, b
    real(dp) :: ma, mb

    Dmat = (0.0_dp, 0.0_dp)
    do a = 1, dim_d
      ma = j - real(a-1, dp)
      do b = 1, dim_d
        mb = j - real(b-1, dp)
        if (abs(mb + ma) < 1.0d-9) then
          Dmat(a,b) = cmplx((-1.0_dp)**nint(j-ma), 0.0_dp, dp)
        end if
      end do
    end do
  end subroutine time_reversal_matrix

  subroutine build_Aq(ja, jb, r, q, s_phys, Maa, Mab, Mba, Mbb, Aq, D)
    real(dp), intent(in) :: ja, jb, s_phys
    integer, intent(in) :: r, q, D
    complex(dp), intent(in) :: Maa(r,r), Mab(r,r), Mba(r,r), Mbb(r,r)
    complex(dp), intent(out) :: Aq(D,D)
    integer :: da, db
    complex(dp), allocatable :: Caa(:,:), Cab(:,:), Cba(:,:), Cbb(:,:)
    complex(dp), allocatable :: Baa(:,:), Bab(:,:), Bba(:,:), Bbb(:,:)

    da = nint(2*ja+1); db = nint(2*jb+1)
    allocate(Caa(da,da), Cab(da,db), Cba(db,da), Cbb(db,db))
    allocate(Baa(da*r,da*r), Bab(da*r,db*r), Bba(db*r,da*r), Bbb(db*r,db*r))

    call cg_matrix(ja, ja, q, s_phys, Caa, da, da)
    call cg_matrix(ja, jb, q, s_phys, Cab, da, db)
    call cg_matrix(jb, ja, q, s_phys, Cba, db, da)
    call cg_matrix(jb, jb, q, s_phys, Cbb, db, db)

    call kron(Caa, da, da, Maa, r, r, Baa)
    call kron(Cab, da, db, Mab, r, r, Bab)
    call kron(Cba, db, da, Mba, r, r, Bba)
    call kron(Cbb, db, db, Mbb, r, r, Bbb)

    Aq(1:da*r, 1:da*r)     = Baa
    Aq(1:da*r, da*r+1:D)   = Bab
    Aq(da*r+1:D, 1:da*r)   = Bba
    Aq(da*r+1:D, da*r+1:D) = Bbb

    deallocate(Caa, Cab, Cba, Cbb, Baa, Bab, Bba, Bbb)
  end subroutine build_Aq

  subroutine build_all_Aq(ja, jb, r, s_phys, Maa, Mab, Mba, Mbb, Aall, D)
    real(dp), intent(in) :: ja, jb, s_phys
    integer, intent(in) :: r, D
    complex(dp), intent(in) :: Maa(r,r), Mab(r,r), Mba(r,r), Mbb(r,r)
    complex(dp), intent(out) :: Aall(D,D,*)
    integer :: Sph, nq, idx, q

    Sph = nint(s_phys); nq = 2*Sph + 1
    do idx = 1, nq
      q = Sph - (idx-1)
      call build_Aq(ja, jb, r, q, s_phys, Maa, Mab, Mba, Mbb, Aall(:,:,idx), D)
    end do
  end subroutine build_all_Aq

  subroutine kron_add_block(A, na, B, nb, C, dimC, row_off, col_off)
    integer, intent(in) :: na, nb, dimC, row_off, col_off
    complex(dp), intent(in) :: A(na,na), B(nb,nb)
    complex(dp), intent(inout) :: C(dimC,dimC)
    integer :: i, j, p, q

    do i = 1, na
      do j = 1, na
        do p = 1, nb
          do q = 1, nb
            C(row_off + (i-1)*nb + p - 1, col_off + (j-1)*nb + q - 1) = A(i,j)*B(p,q)
          end do
        end do
      end do
    end do
  end subroutine kron_add_block

  subroutine build_spin_generators(j, Jz, Jp, Jm, dim_j)
    real(dp), intent(in) :: j
    integer, intent(in) :: dim_j
    complex(dp), intent(out) :: Jz(dim_j,dim_j), Jp(dim_j,dim_j), Jm(dim_j,dim_j)
    integer :: a, b
    real(dp) :: ma, mb

    Jz = (0.0_dp, 0.0_dp); Jp = (0.0_dp, 0.0_dp); Jm = (0.0_dp, 0.0_dp)
    do a = 1, dim_j
      ma = j - real(a-1, dp)
      Jz(a,a) = cmplx(ma, 0.0_dp, dp)
      do b = 1, dim_j
        mb = j - real(b-1, dp)
        if (abs(ma - (mb + 1.0_dp)) < 1.0d-9) then
          Jp(a,b) = cmplx(sqrt(j*(j+1.0_dp) - mb*(mb+1.0_dp)), 0.0_dp, dp)
        end if
        if (abs(ma - (mb - 1.0_dp)) < 1.0d-9) then
          Jm(a,b) = cmplx(sqrt(j*(j+1.0_dp) - mb*(mb-1.0_dp)), 0.0_dp, dp)
        end if
      end do
    end do
  end subroutine build_spin_generators

  subroutine check_su2_symmetry(ja, jb, r, s_phys, Aall, D, err)
    real(dp), intent(in) :: ja, jb, s_phys
    integer, intent(in) :: r, D
    complex(dp), intent(in) :: Aall(D,D,*)
    real(dp), intent(out) :: err
    integer :: dim_a, dim_b, k, Sph, nq, idx, q
    complex(dp), allocatable :: Jz_a(:,:), Jz_b(:,:), Jp_a(:,:), Jp_b(:,:), Jm_a(:,:), Jm_b(:,:)
    complex(dp), allocatable :: Ir(:,:)
    complex(dp) :: Jz(D,D), Jp(D,D), Jm(D,D), comm(D,D), expected(D,D)
    real(dp) :: max_e, coefp, coefm

    dim_a = nint(2*ja+1); dim_b = nint(2*jb+1)
    allocate(Jz_a(dim_a,dim_a), Jz_b(dim_b,dim_b))
    allocate(Jp_a(dim_a,dim_a), Jp_b(dim_b,dim_b))
    allocate(Jm_a(dim_a,dim_a), Jm_b(dim_b,dim_b))
    allocate(Ir(r,r))

    Ir = (0.0_dp, 0.0_dp)
    do k = 1, r; Ir(k,k) = (1.0_dp, 0.0_dp); end do

    call build_spin_generators(ja, Jz_a, Jp_a, Jm_a, dim_a)
    call build_spin_generators(jb, Jz_b, Jp_b, Jm_b, dim_b)

    Jz = (0.0_dp, 0.0_dp); Jp = (0.0_dp, 0.0_dp); Jm = (0.0_dp, 0.0_dp)
    call kron_add_block(Jz_a, dim_a, Ir, r, Jz, D, 1, 1)
    call kron_add_block(Jz_b, dim_b, Ir, r, Jz, D, dim_a*r+1, dim_a*r+1)
    call kron_add_block(Jp_a, dim_a, Ir, r, Jp, D, 1, 1)
    call kron_add_block(Jp_b, dim_b, Ir, r, Jp, D, dim_a*r+1, dim_a*r+1)
    call kron_add_block(Jm_a, dim_a, Ir, r, Jm, D, 1, 1)
    call kron_add_block(Jm_b, dim_b, Ir, r, Jm, D, dim_a*r+1, dim_a*r+1)

    Sph = nint(s_phys); nq = 2*Sph + 1
    max_e = 0.0_dp

    do idx = 1, nq
      q = Sph - (idx-1)

      comm = matmul(Jz, Aall(:,:,idx)) - matmul(Aall(:,:,idx), Jz)
      max_e = max(max_e, maxval(abs(comm - real(q,dp)*Aall(:,:,idx))))

      coefp = s_phys*(s_phys+1.0_dp) - real(q,dp)*real(q+1,dp)
      if (coefp < 0.0_dp) coefp = 0.0_dp
      comm = matmul(Jp, Aall(:,:,idx)) - matmul(Aall(:,:,idx), Jp)
      if (q < Sph) then
        expected = sqrt(coefp) * Aall(:,:,idx-1)
      else
        expected = (0.0_dp, 0.0_dp)
      end if
      max_e = max(max_e, maxval(abs(comm - expected)))

      coefm = s_phys*(s_phys+1.0_dp) - real(q,dp)*real(q-1,dp)
      if (coefm < 0.0_dp) coefm = 0.0_dp
      comm = matmul(Jm, Aall(:,:,idx)) - matmul(Aall(:,:,idx), Jm)
      if (q > -Sph) then
        expected = sqrt(coefm) * Aall(:,:,idx+1)
      else
        expected = (0.0_dp, 0.0_dp)
      end if
      max_e = max(max_e, maxval(abs(comm - expected)))
    end do

    err = max_e
    deallocate(Jz_a, Jz_b, Jp_a, Jp_b, Jm_a, Jm_b, Ir)
  end subroutine check_su2_symmetry

  subroutine fixed_inputs(Z)
    complex(dp), intent(out) :: Z(4,4,3)
    real(dp) :: re(4,4), im(4,4)
    re = reshape([ 2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, -3.0_dp, -1.0_dp, -2.0_dp ], [4,4])
    im = reshape([ 1.0_dp, -1.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 0.0_dp, 3.0_dp, -2.0_dp, -3.0_dp, -3.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, -1.0_dp ], [4,4])
    Z(:,:,1) = cmplx(re, im, dp)
    re = reshape([ 2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -3.0_dp, 3.0_dp, 1.0_dp, 3.0_dp, -3.0_dp ], [4,4])
    im = reshape([ 0.0_dp, 2.0_dp, 3.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 3.0_dp, -2.0_dp, 1.0_dp ], [4,4])
    Z(:,:,2) = cmplx(re, im, dp)
    re = reshape([ 3.0_dp, -3.0_dp, -3.0_dp, -3.0_dp, 3.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, -1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, -2.0_dp ], [4,4])
    im = reshape([ 2.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 1.0_dp, -2.0_dp, -1.0_dp, -2.0_dp, -2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, -2.0_dp, 3.0_dp ], [4,4])
    Z(:,:,3) = cmplx(re, im, dp)
  end subroutine fixed_inputs

  subroutine make_case(h, t, reflect, A, W, R, phase, ja, jb)
    integer, intent(in) :: h, t, reflect
    complex(dp), allocatable, intent(out) :: A(:,:,:), W(:,:), R(:,:)
    integer, intent(out) :: phase
    real(dp), intent(out) :: ja, jb
    complex(dp) :: J(4,4), K(4,4), Z(4,4,3), X(4,4,3)
    complex(dp) :: aa(4,4), ab(4,4), ba(4,4), bb(4,4)
    complex(dp), allocatable :: Ua(:,:), Ub(:,:)
    integer :: n, da, db, D, eps
    if (h == -1) then
      ja=0.5_dp; jb=1.5_dp
    else
      ja=1.0_dp; jb=2.0_dp
    end if
    eps=h*reflect; phase=-eps
    J=(0.0_dp,0.0_dp); K=J
    if (h*t == 1) then
      do n=1,4
        J(n,n)=1.0_dp
      end do
    else
      J(1,2)=1.0_dp; J(2,1)=-1.0_dp
      J(3,4)=1.0_dp; J(4,3)=-1.0_dp
    end if
    if (eps == 1) then
      K(1,1)=1.0_dp; K(2,2)=1.0_dp
      K(3,3)=-1.0_dp; K(4,4)=-1.0_dp
    else
      K(1,2)=1.0_dp; K(2,1)=-1.0_dp
      K(3,4)=-1.0_dp; K(4,3)=1.0_dp
    end if
    call fixed_inputs(Z)
    do n=1,3
      X(:,:,n)=Z(:,:,n)+matmul(J,matmul(conjg(Z(:,:,n)),transpose(J)))
    end do
    aa=X(:,:,1)-phase*matmul(K,matmul(transpose(X(:,:,1)),transpose(K)))
    bb=X(:,:,2)-phase*matmul(K,matmul(transpose(X(:,:,2)),transpose(K)))
    ab=X(:,:,3)
    ba=phase*sqrt((2*ja+1)/(2*jb+1))*matmul(K,matmul(transpose(ab),transpose(K)))
    da=nint(2*ja+1); db=nint(2*jb+1); D=4*(da+db)
    allocate(A(D,D,3),W(D,D),R(D,D),Ua(da,da),Ub(db,db))
    call build_all_Aq(ja,jb,4,1.0_dp,aa,ab,ba,bb,A,D)
    call time_reversal_matrix(ja,Ua,da)
    call time_reversal_matrix(jb,Ub,db)
    W=(0.0_dp,0.0_dp); R=W
    call kron_add_block(Ua,da,J,4,W,D,1,1)
    call kron_add_block(Ub,db,J,4,W,D,4*da+1,4*da+1)
    call kron_add_block(Ua,da,K,4,R,D,1,1)
    call kron_add_block(Ub,db,K,4,R,D,4*da+1,4*da+1)
  end subroutine make_case

  subroutine verify_indices(A,W,R,D,t,reflect,phase,err)
    integer, intent(in) :: D,t,reflect,phase
    complex(dp), intent(in) :: A(D,D,3),W(D,D),R(D,D)
    real(dp), intent(out) :: err
    complex(dp) :: Id(D,D), res(D,D)
    integer :: i
    Id=(0.0_dp,0.0_dp)
    do i=1,D
      Id(i,i)=1.0_dp
    end do
    err=maxval(abs(matmul(W,conjg(W))-t*Id))
    err=max(err,maxval(abs(transpose(R)-reflect*R)))
    do i=1,3
      res=matmul(W,matmul(conjg(A(:,:,i)),conjg(transpose(W))))-(-1)**(i-1)*A(:,:,4-i)
      err=max(err,maxval(abs(res)))
      ! Ordinary reflection, without q -> -q or an extra physical rotation.
      res=transpose(A(:,:,i))-phase*matmul(conjg(transpose(R)),matmul(A(:,:,i),R))
      err=max(err,maxval(abs(res)))
    end do
  end subroutine verify_indices

  subroutine fixed_length_rank(A,D,nmax,rank_out,L_out)
    integer, intent(in) :: D,nmax
    complex(dp), intent(in) :: A(D,D,3)
    integer, intent(out) :: rank_out,L_out
    complex(dp), allocatable :: basis(:,:), cand(:,:), umat(:,:), work(:)
    complex(dp) :: vt(1,1), query(1), identity(D,D), prod(D,D)
    real(dp), allocatable :: sv(:), rw(:)
    integer :: DD,nb,ncols,mn,L,k,q,i,info,lwork
    DD=D*D; nb=1; rank_out=0; L_out=-1
    allocate(basis(DD,1));identity=(0.0_dp,0.0_dp)
    do i=1,D
      identity(i,i)=1.0_dp/sqrt(real(D,dp))
    end do
    basis(:,1)=reshape(identity,[DD])
    do L=1,nmax
      ncols=3*nb;mn=min(DD,ncols)
      allocate(cand(DD,ncols),umat(DD,mn),sv(mn),rw(5*mn))
      do k=1,nb
        do q=1,3
          prod=matmul(reshape(basis(:,k),[D,D]),A(:,:,q))
          cand(:,3*(k-1)+q)=reshape(prod,[DD])
        end do
      end do
      ! No copy of basis is retained as a candidate: all words have length L.
      call zgesvd('S','N',DD,ncols,cand,DD,sv,umat,DD,vt,1,query,-1,rw,info)
      if (info /= 0) error stop 'SVD workspace query failed'
      lwork=max(1,nint(real(query(1),dp)));allocate(work(lwork))
      call zgesvd('S','N',DD,ncols,cand,DD,sv,umat,DD,vt,1,work,lwork,rw,info)
      if (info /= 0) error stop 'SVD failed'
      rank_out=0
      if (sv(1)>0) rank_out=count(sv>1.0d-10*sv(1))
      deallocate(basis);nb=rank_out
      allocate(basis(DD,nb));basis=umat(:,1:nb)
      deallocate(cand,umat,sv,rw,work)
      if (rank_out == DD) then
        L_out=L
        exit
      end if
      if (rank_out == 0) exit
    end do
  end subroutine fixed_length_rank

  subroutine spectrum(A,D,gap,xi)
    integer, intent(in) :: D
    complex(dp), intent(in) :: A(D,D,3)
    real(dp), intent(out) :: gap,xi
    complex(dp), allocatable :: E(:,:), term(:,:), vals(:),work(:)
    complex(dp) :: vl(1,1), vr(1,1), query(1)
    real(dp), allocatable :: rw(:), mags(:)
    real(dp) :: lead,second,ratio
    integer :: DD,q,info,lwork,k
    DD=D*D
    allocate(E(DD,DD),term(DD,DD),vals(DD),rw(2*DD),mags(DD))
    E=(0.0_dp,0.0_dp)
    do q=1,3
      call kron(conjg(A(:,:,q)),D,D,A(:,:,q),D,D,term)
      E=E+term
    end do
    call zgeev('N','N',DD,E,DD,vals,vl,1,vr,1,query,-1,rw,info)
    if (info /= 0) error stop 'Eigenvalue workspace query failed'
    lwork=max(1,nint(real(query(1),dp)));allocate(work(lwork))
    call zgeev('N','N',DD,E,DD,vals,vl,1,vr,1,work,lwork,rw,info)
    if (info /= 0) error stop 'Eigenvalue computation failed'
    mags=abs(vals);k=maxloc(mags,dim=1);lead=mags(k);mags(k)=0;second=maxval(mags)
    if (lead <= 0) error stop 'Zero transfer radius'
    ratio=second/lead;gap=1-ratio
    if (gap <= 1.0d-10) error stop 'Unresolved transfer spectral gap'
    xi=0
    if (ratio>0) xi=-1/log(ratio)
  end subroutine spectrum
end module verified_mps

program run_verified_spin1
  use verified_mps
  implicit none
  complex(dp), allocatable :: A(:,:,:),W(:,:),R(:,:)
  integer :: h,t,reflect,phase,D,rank_now,L
  real(dp) :: ja,jb,err_su2,err_indices,gap,xi
  print *, 'Spin 1, multiplicity 4; actual (Tasaki,time reversal,reflection) indices'
  print *, 'Numerical fixed-length checks; exact certificates: indices_verified.py'
  do h=-1,1,2
    do t=-1,1,2
      do reflect=-1,1,2
        call make_case(h,t,reflect,A,W,R,phase,ja,jb)
        D=size(A,1)
        call check_su2_symmetry(ja,jb,4,1.0_dp,A,D,err_su2)
        call verify_indices(A,W,R,D,t,reflect,phase,err_indices)
        if (max(err_su2,err_indices)>1.0d-9) error stop 'Symmetry or index failure'
        call fixed_length_rank(A,D,9,rank_now,L)
        if (L<1 .or. rank_now/=D*D) error stop 'Normality not established'
        call spectrum(A,D,gap,xi)
        print '(3(I3),A,I3,A,I5,A,I2,A,F12.7,A,F12.5)', &
          h,t,reflect,' D=',D,' rank=',rank_now,' L=',L,' transfer gap=',gap,' xi=',xi
        deallocate(A,W,R)
      end do
    end do
  end do
end program run_verified_spin1
