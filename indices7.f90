! sample_phases_sn.f90
!
! Random sampling of spin-S (integer S >= 1) SO(3) x Z2^T x Z2^R MPS representatives with multiplicity
! dimension n (any n >= 1), generalizing the s = 1, n = 4 ansatz of "Simple Matrix Product State
! Representatives of Symmetry-Protected Phases".
!
! VIRTUAL IRREPS are chosen minimal: (alpha,beta) with beta = alpha + 1, alpha of the parity set by sigma_H
! (half-integer for sigma_H = -1, integer for sigma_H = +1), and alpha the smallest value with
!     S in alpha x alpha,  S in alpha x beta,  S in beta x beta     (i.e. S <= 2 alpha, |beta-alpha| <= S).
! In doubled integers: sigma_H = -1: 2 alpha = S (S odd) or S+1 (S even);  sigma_H = +1: 2 alpha = S (S even) or S+1 (S odd).
!
! Tensor:  A^q = [[C^q_aa (x) M_aa, C^q_ab (x) M_ab],[C^q_ba (x) M_ba, C^q_bb (x) M_bb]],  q = S,...,-S,
! C = SU(2) Clebsch-Gordan blocks, M_gg' = complex n x n multiplicity matrices (4 n^2 complex = 8 n^2 real parameters).
!
! SYMMETRY CONSTRAINTS (solved numerically, so the sign conventions that change with S are not hard-coded):
!     W conj(A^q) W^-1 = (-1)^(S-q) A^(-q),     (A^q)^T = p R^-1 A^q R,
!     W = diag(U_a (x) J, U_b (x) J),  R = diag(U_a (x) K, U_b (x) K),   (U_j)_{m,n} = (-1)^(j-m) delta_{m,-n},
! with J, K real orthogonal n x n matrices, J conj(J) = sigma_J, K^T = sigma_K K, J K = K J, and
! (sigma_H, sigma_T, sigma_R) = (sigma_H, sigma_H sigma_J, sigma_H sigma_K).
! For every triple the program tries ALL of these (J,K) candidates and p = +1,-1:
!     sigma_J = +1: J = I                       sigma_J = -1: J = blockdiag(j,...,j)           (n even)
!     sigma_K = +1: K = I  or  K = diag(I_ceil(n/2), -I_floor(n/2))
!     sigma_K = -1: K = blockdiag(j,-j,j,...)                                                  (n even)
! (candidates with J K /= K J are discarded; n = 1 therefore reaches only sigma_T = sigma_R = sigma_H).
! Random integer matrices (entries in {-3..3}, real and imaginary parts) are projected onto the solution space
! (the paper's "random integer Z, then project onto the constraints" search).  The SAME raw integers are used for
! every triple, candidate and p.
!
! For each sample the transfer operator E(X) = sum_q A^q X A^q^dagger is split EXACTLY into SU(2) spin sectors
! (blocks of dimension <= 4 n^2).  A sample counts as NORMAL (primitive transfer operator) if the leading eigenvalue is
! simple, unique in modulus, lies in the J=0 sector, and its fixed point is faithful.
! xi = -1/ln(|lambda_2|/|lambda_1|) is reported for normal samples only.
!
! Everything is floating point: no prime or ring is needed.  (For an EXACT certificate one needs a prime p for which
! -1 and all primes <= (2 beta + S + 1) are quadratic residues; the program prints a conservative such prime.)
!
! Usage:  ./sample_phases_sn S N NSAMPLES SEED CHECK
!   S          on-site spin (integer >= 1).  Tested for S = 1, 2, 3 (with the previous, non-minimal irreps)
!   N          dimension of the multiplicity space (n = 1, 2, 4, 8, ...)
!   NSAMPLES   number of random samples
!   SEED       seed of the built-in xoshiro256** generator (NOT numpy's generator)
!   CHECK      1 = also compute a numerical injectivity length (SVD, tol 1e-10) for the best sample per triple
!              (very expensive for large D); 0 = skip
!
! Build:  gfortran -O2 -fwrapv -ffree-line-length-none -o sample_phases_sn sample_phases_sn.f90 -llapack -lblas
!
! Cost: D = n (2 alpha + 2 beta + 2).  The constraint solve stores ~ 8 n^2 x (3S+2) x 2 D^2 doubles
! (S=1, n=4: 25 MB; S=1, n=8: ~0.4 GB; S=3, n=4: ~0.2 GB).

module cg_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains
  real(dp) function fact(n)
    integer, intent(in) :: n
    integer :: i
    fact = 1.0_dp
    do i = 2, n
      fact = fact * real(i, dp)
    end do
  end function fact

  ! Clebsch-Gordan <j1 m1; j2 m2 | Jtot Mtot>, Racah formula
  real(dp) function cg(j1, m1, j2, m2, Jtot, Mtot)
    real(dp), intent(in) :: j1, m1, j2, m2, Jtot, Mtot
    real(dp) :: delta, prefac, sumterm, term
    integer :: kmin, kmax, k, a1, a2, a3, a4, a5, a6, a7, a8, a9
    cg = 0.0_dp
    if (abs(Mtot - (m1+m2)) > 1.0d-9) return
    if (Jtot > j1+j2+1.0d-9 .or. Jtot < abs(j1-j2)-1.0d-9) return
    a1 = nint(j1+j2-Jtot); a2 = nint(j1-j2+Jtot); a3 = nint(-j1+j2+Jtot)
    a4 = nint(j1+j2+Jtot+1.0_dp)
    if (a1 < 0 .or. a2 < 0 .or. a3 < 0) return
    delta = sqrt(fact(a1)*fact(a2)*fact(a3) / fact(a4))
    a5 = nint(Jtot+Mtot); a6 = nint(Jtot-Mtot)
    a7 = nint(j1-m1); a8 = nint(j1+m1); a9 = nint(j2-m2)
    if (a5 < 0 .or. a6 < 0 .or. a7 < 0 .or. a8 < 0 .or. a9 < 0) return
    if (nint(j2+m2) < 0) return
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

  ! (C^q)_{a,b} = <jin,mi; s,q | jout,mo>, rows a: jout, columns b: jin
  subroutine cg_matrix(jout, jin, sphys, q, C, dout, din)
    real(dp), intent(in) :: jout, jin
    integer, intent(in) :: sphys, q, dout, din
    real(dp), intent(out) :: C(dout,din)
    integer :: a, b
    do a = 1, dout
      do b = 1, din
        C(a,b) = cg(jin, jin-real(b-1,dp), real(sphys,dp), real(q,dp), jout, jout-real(a-1,dp))
      end do
    end do
  end subroutine cg_matrix

  ! spin-J tensor-operator component (mu=0) from V_d to V_g:  T(m,n) = <d,n; J,0 | g,m>
  subroutine tensor_op(jg, jd, Jt, T, dg, dd)
    real(dp), intent(in) :: jg, jd, Jt
    integer, intent(in) :: dg, dd
    real(dp), intent(out) :: T(dg,dd)
    integer :: a, b
    do a = 1, dg
      do b = 1, dd
        T(a,b) = cg(jd, jd-real(b-1,dp), Jt, 0.0_dp, jg, jg-real(a-1,dp))
      end do
    end do
  end subroutine tensor_op
end module cg_mod

module rng_mod
  implicit none
  integer, parameter :: i8 = selected_int_kind(18)
  integer(i8) :: st(4)
contains
  integer(i8) function rotl(x, k)
    integer(i8), intent(in) :: x
    integer, intent(in) :: k
    rotl = ior(ishft(x, k), ishft(x, k-64))
  end function rotl

  integer(i8) function splitmix(x)
    integer(i8), intent(inout) :: x
    integer(i8) :: z, c1, c2, c3
    c1 = ior(ishft(int(z'9E3779B9', i8), 32), int(z'7F4A7C15', i8))
    c2 = ior(ishft(int(z'BF58476D', i8), 32), int(z'1CE4E5B9', i8))
    c3 = ior(ishft(int(z'94D049BB', i8), 32), int(z'133111EB', i8))
    x = x + c1
    z = x
    z = ieor(z, ishft(z, -30)) * c2
    z = ieor(z, ishft(z, -27)) * c3
    splitmix = ieor(z, ishft(z, -31))
  end function splitmix

  subroutine rng_seed(seed)
    integer, intent(in) :: seed
    integer(i8) :: x
    integer :: i
    x = int(seed, i8)
    do i = 1, 4
      st(i) = splitmix(x)
    end do
  end subroutine rng_seed

  integer(i8) function rng_next()   ! xoshiro256**
    integer(i8) :: t
    rng_next = rotl(st(2)*5_i8, 7) * 9_i8
    t = ishft(st(2), 17)
    st(3) = ieor(st(3), st(1))
    st(4) = ieor(st(4), st(2))
    st(2) = ieor(st(2), st(3))
    st(1) = ieor(st(1), st(4))
    st(3) = ieor(st(3), t)
    st(4) = rotl(st(4), 45)
  end function rng_next

  integer function rand_m3_p3()      ! unbiased integer in {-3,...,3}
    integer(i8) :: u
    do
      u = iand(ishft(rng_next(), -40), 7_i8)
      if (u < 7_i8) exit
    end do
    rand_m3_p3 = int(u) - 3
  end function rand_m3_p3
end module rng_mod

module mps_mod
  use cg_mod
  implicit none
  integer :: NPAR = 128                  ! 8 n^2 real parameters (4 blocks of n x n complex), set in init_family
  ! Clebsch-Gordan blocks of the current family, filled by init_family
  real(dp), allocatable :: Caa(:,:,:), Cab(:,:,:), Cba(:,:,:), Cbb(:,:,:)
  integer :: fam_s, fam_da, fam_db, fam_dim, fam_nph, fam_n
  real(dp) :: fam_ja, fam_jb
contains
  subroutine init_family(s, ja, jb, n)
    integer, intent(in) :: s, n
    real(dp), intent(in) :: ja, jb
    integer :: qi, q
    fam_s = s; fam_ja = ja; fam_jb = jb; fam_nph = 2*s+1; fam_n = n; NPAR = 8*n*n
    fam_da = nint(2*ja+1); fam_db = nint(2*jb+1); fam_dim = n*(fam_da+fam_db)
    if (allocated(Caa)) deallocate(Caa, Cab, Cba, Cbb)
    allocate(Caa(fam_da,fam_da,fam_nph), Cab(fam_da,fam_db,fam_nph), Cba(fam_db,fam_da,fam_nph), Cbb(fam_db,fam_db,fam_nph))
    do qi = 1, fam_nph
      q = s - (qi-1)
      call cg_matrix(ja, ja, s, q, Caa(:,:,qi), fam_da, fam_da)
      call cg_matrix(ja, jb, s, q, Cab(:,:,qi), fam_da, fam_db)
      call cg_matrix(jb, ja, s, q, Cba(:,:,qi), fam_db, fam_da)
      call cg_matrix(jb, jb, s, q, Cbb(:,:,qi), fam_db, fam_db)
    end do
  end subroutine init_family

  subroutine kron_block(C, nr, nc, M, out)       ! out = C (x) M, C real nr x nc, M n x n
    integer, intent(in) :: nr, nc
    real(dp), intent(in) :: C(nr,nc)
    complex(dp), intent(in) :: M(fam_n,fam_n)
    complex(dp), intent(out) :: out(fam_n*nr,fam_n*nc)
    integer :: i, j, n
    n = fam_n
    do i = 1, nr
      do j = 1, nc
        out(n*(i-1)+1:n*i, n*(j-1)+1:n*j) = C(i,j)*M
      end do
    end do
  end subroutine kron_block

  ! parameter vector x -> Kraus matrices A(:,:,q index);  block order: aa, ab, ba, bb
  subroutine build_tensor(x, A)
    real(dp), intent(in) :: x(NPAR)
    complex(dp), intent(out) :: A(fam_dim,fam_dim,fam_nph)
    complex(dp), allocatable :: Mb(:,:,:)
    integer :: blk, i, j, base, qi, na, nb, n, n2
    n = fam_n; n2 = n*n
    allocate(Mb(n,n,4))
    do blk = 1, 4
      base = 2*n2*(blk-1)
      do i = 1, n
        do j = 1, n
          Mb(i,j,blk) = cmplx(x(base+(i-1)*n+j), x(base+n2+(i-1)*n+j), dp)
        end do
      end do
    end do
    na = n*fam_da; nb = n*fam_db
    do qi = 1, fam_nph
      call kron_block(Caa(:,:,qi), fam_da, fam_da, Mb(:,:,1), A(1:na, 1:na, qi))
      call kron_block(Cab(:,:,qi), fam_da, fam_db, Mb(:,:,2), A(1:na, na+1:fam_dim, qi))
      call kron_block(Cba(:,:,qi), fam_db, fam_da, Mb(:,:,3), A(na+1:fam_dim, 1:na, qi))
      call kron_block(Cbb(:,:,qi), fam_db, fam_db, Mb(:,:,4), A(na+1:fam_dim, na+1:fam_dim, qi))
    end do
  end subroutine build_tensor

  ! U_j(a,b) = (-1)^(a-1) delta_{b, 2j+2-a}   (1-based; (U_j)_{m,n} = (-1)^(j-m) delta_{m,-n})
  subroutine umat(j, U)
    real(dp), intent(in) :: j
    real(dp), allocatable, intent(out) :: U(:,:)
    integer :: d, a
    d = nint(2*j+1); allocate(U(d,d)); U = 0.0_dp
    do a = 1, d
      U(a, d+1-a) = real((-1)**(a-1), dp)
    end do
  end subroutine umat

  ! Candidate (J,K).  jtype: 1 = I (sigma_J=+1), 2 = blockdiag(j,..,j) (sigma_J=-1, n even)
  !                    ktype: 1 = I, 2 = diag(I_ceil(n/2),-I_floor(n/2)) (sigma_K=+1), 3 = blockdiag(j,-j,j,..) (sigma_K=-1, n even)
  subroutine make_jk(n, jtype, ktype, Jm, Km, ok)
    integer, intent(in) :: n, jtype, ktype
    real(dp), intent(out) :: Jm(n,n), Km(n,n)
    logical, intent(out) :: ok
    integer :: i, b
    ok = .true.; Jm = 0.0_dp; Km = 0.0_dp
    if (jtype == 1) then
      do i = 1, n
        Jm(i,i) = 1.0_dp
      end do
    else
      if (mod(n,2) /= 0) then; ok = .false.; return; end if
      do b = 1, n/2
        Jm(2*b-1,2*b) = 1.0_dp; Jm(2*b,2*b-1) = -1.0_dp
      end do
    end if
    if (ktype == 1) then
      do i = 1, n
        Km(i,i) = 1.0_dp
      end do
    else if (ktype == 2) then
      if (n == 1) then; ok = .false.; return; end if
      do i = 1, n
        Km(i,i) = 1.0_dp
        if (i > (n+1)/2) Km(i,i) = -1.0_dp
      end do
    else
      if (mod(n,2) /= 0) then; ok = .false.; return; end if
      do b = 1, n/2
        Km(2*b-1,2*b) = 1.0_dp*(-1)**(b-1); Km(2*b,2*b-1) = -1.0_dp*(-1)**(b-1)
      end do
    end if
    if (maxval(abs(matmul(Jm,Km)-matmul(Km,Jm))) > 1.0d-12) ok = .false.
  end subroutine make_jk

  ! minimal irreps: sigma_H = -1 half-integer, +1 integer; beta = alpha + 1; S in alpha x alpha  <=>  2 alpha >= S
  subroutine pick_irreps(S, h, ja, jb)
    integer, intent(in) :: S, h
    real(dp), intent(out) :: ja, jb
    integer :: a2
    if (h == -1) then
      a2 = S; if (mod(S,2) == 0) a2 = S + 1
    else
      a2 = S; if (mod(S,2) == 1) a2 = S + 1
    end if
    ja = 0.5_dp*real(a2,dp); jb = ja + 1.0_dp
  end subroutine pick_irreps

  subroutine build_WR(Jm, Km, W, R)
    real(dp), intent(in) :: Jm(fam_n,fam_n), Km(fam_n,fam_n)
    complex(dp), intent(out) :: W(fam_dim,fam_dim), R(fam_dim,fam_dim)
    real(dp), allocatable :: Ua(:,:), Ub(:,:)
    integer :: i, j, k, l, na, nb, n
    n = fam_n
    call umat(fam_ja, Ua); call umat(fam_jb, Ub)
    na = n*fam_da; nb = n*fam_db
    W = (0.0_dp, 0.0_dp); R = (0.0_dp, 0.0_dp)
    do i = 1, fam_da
      do j = 1, fam_da
        do k = 1, n
          do l = 1, n
            W(n*(i-1)+k, n*(j-1)+l) = Ua(i,j)*Jm(k,l)
            R(n*(i-1)+k, n*(j-1)+l) = Ua(i,j)*Km(k,l)
          end do
        end do
      end do
    end do
    do i = 1, fam_db
      do j = 1, fam_db
        do k = 1, n
          do l = 1, n
            W(na+n*(i-1)+k, na+n*(j-1)+l) = Ub(i,j)*Jm(k,l)
            R(na+n*(i-1)+k, na+n*(j-1)+l) = Ub(i,j)*Km(k,l)
          end do
        end do
      end do
    end do
  end subroutine build_WR

  ! Real residual vector of the symmetry constraints (linear in x).  Time reversal is imposed for q >= 0
  ! only (the q <= 0 equations are equivalent).
  subroutine residual(x, W, R, p, nr, out)
    real(dp), intent(in) :: x(NPAR), p
    complex(dp), intent(in) :: W(fam_dim,fam_dim), R(fam_dim,fam_dim)
    integer, intent(in) :: nr
    real(dp), intent(out) :: out(nr)
    complex(dp), allocatable :: A(:,:,:), tmp(:,:)
    integer :: idx, pos, DD, q
    DD = fam_dim*fam_dim
    allocate(A(fam_dim,fam_dim,fam_nph), tmp(fam_dim,fam_dim))
    call build_tensor(x, A)
    pos = 0
    do idx = 1, fam_nph
      q = fam_s - (idx-1)
      if (q >= 0) then
        tmp = matmul(W, matmul(conjg(A(:,:,idx)), transpose(W))) - real((-1)**(idx-1),dp)*A(:,:,fam_nph+1-idx)
        out(pos+1:pos+DD) = reshape(real(tmp), [DD]); pos = pos + DD
        out(pos+1:pos+DD) = reshape(aimag(tmp), [DD]); pos = pos + DD
      end if
      tmp = transpose(A(:,:,idx)) - p*matmul(transpose(R), matmul(A(:,:,idx), R))
      out(pos+1:pos+DD) = reshape(real(tmp), [DD]); pos = pos + DD
      out(pos+1:pos+DD) = reshape(aimag(tmp), [DD]); pos = pos + DD
    end do
  end subroutine residual

  ! Orthonormal basis N(:,1:nnull) of the real solution space of the constraints
  subroutine constraint_nullspace(W, R, p, N, nnull)
    complex(dp), intent(in) :: W(fam_dim,fam_dim), R(fam_dim,fam_dim)
    real(dp), intent(in) :: p
    real(dp), allocatable, intent(out) :: N(:,:)
    integer, intent(out) :: nnull
    real(dp), allocatable :: Rm(:,:), G(:,:), ev(:), work(:), e(:)
    integer :: nr, k, info, lwork, i
    real(dp) :: tol
    nr = (fam_s+1)*2*fam_dim*fam_dim + fam_nph*2*fam_dim*fam_dim
    allocate(Rm(nr,NPAR), G(NPAR,NPAR), ev(NPAR), e(NPAR))
    do k = 1, NPAR
      e = 0.0_dp; e(k) = 1.0_dp
      call residual(e, W, R, p, nr, Rm(:,k))
    end do
    G = matmul(transpose(Rm), Rm)
    deallocate(Rm)
    lwork = max(64, 4*NPAR); allocate(work(lwork))
    call dsyev('V','U',NPAR,G,NPAR,ev,work,lwork,info)
    if (info /= 0) error stop 'dsyev failed'
    tol = 1.0d-9*max(1.0_dp, ev(NPAR))
    nnull = count(ev < tol)
    allocate(N(NPAR,max(nnull,1)))
    if (nnull > 0) N(:,1:nnull) = G(:,1:nnull)
    i = 0
  end subroutine constraint_nullspace

  subroutine eig_general(M, n, w, vr, want)
    complex(dp), intent(inout) :: M(n,n)
    integer, intent(in) :: n
    complex(dp), intent(out) :: w(n), vr(n,n)
    logical, intent(in) :: want
    complex(dp) :: vl(1,1), query(1)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rw(:)
    character :: jobvr
    integer :: info, lwork
    jobvr = 'N'; if (want) jobvr = 'V'
    allocate(rw(2*n))
    call zgeev('N',jobvr,n,M,n,w,vl,1,vr,n,query,-1,rw,info)
    lwork = max(1, nint(real(query(1),dp))); allocate(work(lwork))
    call zgeev('N',jobvr,n,M,n,w,vl,1,vr,n,work,lwork,rw,info)
    if (info /= 0) error stop 'zgeev failed'
  end subroutine eig_general

  ! Exact SU(2)-sector analysis of the transfer operator.  Outputs: xi, sectors of lambda_1 and lambda_2,
  ! number of eigenvalues of (almost) maximal modulus, and fixed-point diagnostics in the J=0 sector.
  subroutine classify(A, xi, J1o, J2o, nper, herm, fpr)
    complex(dp), intent(in) :: A(fam_dim,fam_dim,fam_nph)
    real(dp), intent(out) :: xi, herm, fpr
    integer, intent(out) :: J1o, J2o, nper
    integer :: nn_, DM, NPH, Jt, ig, id, npair, nb, c, r, k, l, qq, dg, dd, nbasis, Jmax, ntot, i, m, info
    integer :: pg(4), pd(4), dims(2), offs(2), imax
    real(dp) :: js(2), jg, jd, l1, l2, hn, mn, ph_re, ph_im
    real(dp), allocatable :: T(:,:), nrm(:), magv(:), tmpv(:)
    complex(dp), allocatable :: Bs(:,:,:), EB(:,:), Mj(:,:), w(:), vr(:,:), xfix(:)
    complex(dp), allocatable :: mfix(:,:,:), Hm(:,:), wk(:)
    complex(dp) :: ph
    real(dp), allocatable :: ev4(:), rw4(:)
    real(dp) :: emin, emax
    integer, allocatable :: Jv(:)
    integer :: npair0, pg0(4), pd0(4)
    DM = fam_dim; NPH = fam_nph; nn_ = fam_n
    allocate(mfix(nn_,nn_,2), Hm(nn_,nn_), wk(max(8,4*nn_)), ev4(nn_), rw4(max(3,3*nn_)))
    js = [fam_ja, fam_jb]; dims = [fam_da, fam_db]; offs = [0, nn_*fam_da]
    Jmax = nint(fam_ja+fam_jb)
    allocate(magv(DM*DM), Jv(DM*DM)); ntot = 0
    npair0 = 0
    do Jt = 0, Jmax
      npair = 0
      do ig = 1, 2
        do id = 1, 2
          jg = js(ig); jd = js(id)
          if (abs(jg-jd) <= real(Jt,dp)+1.0d-9 .and. real(Jt,dp) <= jg+jd+1.0d-9) then
            dg = dims(ig); dd = dims(id)
            allocate(T(dg,dd))
            call tensor_op(jg, jd, real(Jt,dp), T, dg, dd)
            if (sum(T*T) > 1.0d-12) then
              npair = npair + 1; pg(npair) = ig; pd(npair) = id
            end if
            deallocate(T)
          end if
        end do
      end do
      if (npair == 0) cycle
      nbasis = nn_*nn_*npair
      allocate(Bs(DM,DM,nbasis), EB(DM,DM), Mj(nbasis,nbasis), w(nbasis), vr(nbasis,nbasis), nrm(nbasis))
      Bs = (0.0_dp, 0.0_dp); nb = 0
      do c = 1, npair
        ig = pg(c); id = pd(c); dg = dims(ig); dd = dims(id)
        allocate(T(dg,dd))
        call tensor_op(js(ig), js(id), real(Jt,dp), T, dg, dd)
        do k = 1, nn_
          do l = 1, nn_
            nb = nb + 1
            do r = 1, dg
              do qq = 1, dd
                Bs(offs(ig)+nn_*(r-1)+k, offs(id)+nn_*(qq-1)+l, nb) = T(r,qq)
              end do
            end do
            nrm(nb) = sum(abs(Bs(:,:,nb))**2)
          end do
        end do
        deallocate(T)
      end do
      do c = 1, nbasis
        EB = (0.0_dp, 0.0_dp)
        do qq = 1, NPH
          EB = EB + matmul(A(:,:,qq), matmul(Bs(:,:,c), conjg(transpose(A(:,:,qq)))))
        end do
        do r = 1, nbasis
          Mj(r,c) = sum(conjg(Bs(:,:,r))*EB) / nrm(r)
        end do
      end do
      call eig_general(Mj, nbasis, w, vr, Jt == 0)
      do c = 1, nbasis
        ntot = ntot + 1; magv(ntot) = abs(w(c)); Jv(ntot) = Jt
      end do
      if (Jt == 0) then
        imax = maxloc(abs(w), dim=1)
        allocate(xfix(nbasis)); xfix = vr(:,imax)
        npair0 = npair; pg0 = pg; pd0 = pd
      end if
      deallocate(Bs, EB, Mj, w, vr, nrm)
    end do
    imax = maxloc(magv(1:ntot), dim=1); l1 = magv(imax); J1o = Jv(imax)
    nper = count(magv(1:ntot) > l1*(1.0_dp-1.0d-7))
    allocate(tmpv(ntot)); tmpv = magv(1:ntot); tmpv(imax) = -1.0_dp
    l2 = maxval(tmpv); J2o = Jv(maxloc(tmpv, dim=1))
    xi = -1.0_dp/log(l2/l1)
    ! fixed point from the J=0 eigenvector (valid when lambda_1 lives in J=0)
    mfix = (0.0_dp, 0.0_dp)
    do c = 1, npair0
      if (pg0(c) == pd0(c)) then
        do k = 1, nn_
          do l = 1, nn_
            mfix(k,l,pg0(c)) = xfix((c-1)*nn_*nn_+(k-1)*nn_+l)
          end do
        end do
      end if
    end do
    ph = (0.0_dp, 0.0_dp)
    do i = 1, 2
      do k = 1, nn_
        ph = ph + (2.0_dp*js(i)+1.0_dp)*mfix(k,k,i)
      end do
    end do
    if (abs(ph) > 0.0_dp) mfix = mfix*conjg(ph)/abs(ph)
    herm = 0.0_dp; mn = 0.0_dp; fpr = huge(1.0_dp)
    do i = 1, 2
      hn = sqrt(sum(abs(mfix(:,:,i)-conjg(transpose(mfix(:,:,i))))**2))
      mn = max(mn, sqrt(sum(abs(mfix(:,:,i))**2)))
      herm = max(herm, hn)
      Hm = 0.5_dp*(mfix(:,:,i)+conjg(transpose(mfix(:,:,i))))
      call zheev('N','U',nn_,Hm,nn_,ev4,wk,max(8,4*nn_),rw4,info)
      emin = ev4(1); emax = ev4(nn_)
      if (abs(emax) > 1.0d-14) then
        fpr = min(fpr, emin/emax)
      else
        fpr = min(fpr, 0.0_dp)
      end if
    end do
    if (mn > 0.0_dp) herm = herm/mn
  end subroutine classify

  logical function is_normal(J1o, nper, herm, fpr)
    integer, intent(in) :: J1o, nper
    real(dp), intent(in) :: herm, fpr
    is_normal = (J1o == 0 .and. nper == 1 .and. herm < 1.0d-6 .and. fpr > 1.0d-9)
  end function is_normal

  ! Numerical injectivity length: smallest L <= nmax with span{A^{i1}..A^{iL}} = M_D (SVD rank, tol 1e-10)
  integer function injectivity_length(A, nmax)
    complex(dp), intent(in) :: A(fam_dim,fam_dim,fam_nph)
    integer, intent(in) :: nmax
    complex(dp), allocatable :: basis(:,:), cand(:,:), umat_(:,:), work(:)
    complex(dp) :: vt(1,1), query(1), Mv(fam_dim,fam_dim)
    real(dp), allocatable :: sv(:), rw(:)
    integer :: DD, DM, NPH, nb, ncols, mn, L, k, q, i, info, lwork, rk
    DM = fam_dim; NPH = fam_nph; DD = DM*DM; nb = 1; injectivity_length = -1
    allocate(basis(DD,1)); Mv = (0.0_dp, 0.0_dp)
    do i = 1, DM
      Mv(i,i) = 1.0_dp/sqrt(real(DM,dp))
    end do
    basis(:,1) = reshape(Mv, [DD])
    do L = 1, nmax
      ncols = NPH*nb; mn = min(DD, ncols)
      allocate(cand(DD,ncols), umat_(DD,mn), sv(mn), rw(5*mn))
      do k = 1, nb
        do q = 1, NPH
          cand(:, NPH*(k-1)+q) = reshape(matmul(reshape(basis(:,k),[DM,DM]), A(:,:,q)), [DD])
        end do
      end do
      call zgesvd('S','N',DD,ncols,cand,DD,sv,umat_,DD,vt,1,query,-1,rw,info)
      lwork = max(1, nint(real(query(1),dp))); allocate(work(lwork))
      call zgesvd('S','N',DD,ncols,cand,DD,sv,umat_,DD,vt,1,work,lwork,rw,info)
      if (info /= 0) error stop 'zgesvd failed'
      rk = 0
      if (sv(1) > 0.0_dp) rk = count(sv > 1.0d-10*sv(1))
      deallocate(basis); nb = rk
      allocate(basis(DD,nb)); basis = umat_(:,1:nb)
      deallocate(cand, umat_, sv, rw, work)
      if (rk == DD) then
        injectivity_length = L; exit
      end if
    end do
  end function injectivity_length

  ! ---- conservative prime for an exact certificate: p > B with -1 and all primes q <= B quadratic residues
  integer(8) function modpow(a, e, p)
    integer(8), intent(in) :: a, e, p
    integer(8) :: base, ex, res
    base = mod(a, p); ex = e; res = 1_8
    do while (ex > 0_8)
      if (mod(ex,2_8) == 1_8) res = mod(res*base, p)
      base = mod(base*base, p); ex = ex/2_8
    end do
    modpow = res
  end function modpow

  logical function isprime8(n)
    integer(8), intent(in) :: n
    integer(8) :: i
    isprime8 = .false.
    if (n < 2_8) return
    if (n < 4_8) then; isprime8 = .true.; return; end if
    if (mod(n,2_8) == 0_8) return
    i = 3_8
    do while (i*i <= n)
      if (mod(n,i) == 0_8) return
      i = i + 2_8
    end do
    isprime8 = .true.
  end function isprime8

  integer(8) function find_prime(B)
    integer, intent(in) :: B
    integer(8) :: p, q
    logical :: ok
    find_prime = -1_8
    p = int(B,8) + 1_8
    do while (p < 3000000000_8)
      if (mod(p,8_8) == 1_8 .and. isprime8(p)) then
        ok = (modpow(p-1_8, (p-1_8)/2_8, p) == 1_8)
        q = 2_8
        do while (ok .and. q <= int(B,8))
          if (isprime8(q)) ok = (modpow(q, (p-1_8)/2_8, p) == 1_8)
          q = q + 1_8
        end do
        if (ok) then; find_prime = p; return; end if
      end if
      p = p + 1_8
    end do
  end function find_prime
end module mps_mod

program sample_phases_sn
  use cg_mod
  use rng_mod
  use mps_mod
  implicit none
  integer :: s, nmult, nsamp, seed, docheck, h, t, r, ip, itrip, isamp, i, j, k, nnull, Jl1, Jl2, nper, L
  integer :: jt, kt, sJ, sK, nkt, ik, B, hh
  integer :: ktl(2), bestsamp(8), bestJ2(8), bestjt(8), bestkt(8), bestip(8), nprim
  real(dp) :: ja, jb, xi, herm, fpr, pv, bestxi(8), ja2, jb2
  character(len=32) :: arg
  character(len=4) :: jname, kname
  complex(dp), allocatable :: Wm(:,:), Rr(:,:), A(:,:,:)
  real(dp), allocatable :: N(:,:), zall(:,:), x(:), Jm(:,:), Km(:,:)
  integer(8) :: pr
  logical :: ok

  s = 1; nmult = 4; nsamp = 10; seed = 1; docheck = 0
  if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) s; end if
  if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) nmult; end if
  if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read(arg,*) nsamp; end if
  if (command_argument_count() >= 4) then; call get_command_argument(4, arg); read(arg,*) seed; end if
  if (command_argument_count() >= 5) then; call get_command_argument(5, arg); read(arg,*) docheck; end if
  if (s < 1) error stop 'S must be an integer >= 1'
  if (nmult < 1) error stop 'N must be >= 1'

  ! conservative exact-certificate prime (information only): primes q <= 2 beta + S + 1 for both parities
  B = 0
  do hh = -1, 1, 2
    call pick_irreps(s, hh, ja2, jb2)
    B = max(B, nint(2.0_dp*jb2) + s + 1)
  end do
  pr = find_prime(B)
  write(*,'(A,I0,A,I0,A,I0)') ' S = ', s, '   multiplicity n = ', nmult, '   conservative ring generators: sqrt(q), primes q <= ', B
  if (pr > 0_8) then
    write(*,'(A,I0,A)') '   smallest prime p > B with -1 and all primes <= B squares mod p: ', pr, &
         '   (what an exact mod-p certificate would use; this program does not)'
  else
    write(*,'(A)') '   no such prime below 3e9 (int64-safe range); exact certificate would need big integers'
  end if

  call rng_seed(seed)
  NPAR = 8*nmult*nmult
  allocate(zall(nsamp,NPAR))
  do isamp = 1, nsamp
    do i = 1, NPAR
      zall(isamp,i) = real(rand_m3_p3(), dp)
    end do
  end do

  bestxi = huge(1.0_dp); bestsamp = 0; bestJ2 = -1; bestjt = 0; bestkt = 0; bestip = 0
  allocate(Jm(nmult,nmult), Km(nmult,nmult))
  itrip = 0
  do h = -1, 1, 2
    call pick_irreps(s, h, ja, jb)
    call init_family(s, ja, jb, nmult)
    allocate(Wm(fam_dim,fam_dim), Rr(fam_dim,fam_dim), A(fam_dim,fam_dim,fam_nph), x(NPAR))
    write(*,*)
    write(*,'(A,I3,A,F5.1,A,F5.1,A,I5,A,I5)') ' sigma_H =', h, '   minimal (alpha,beta) = (', ja, ',', jb, ')   D =', fam_dim, &
          '   real parameters =', NPAR
    do t = -1, 1, 2
      do r = -1, 1, 2
        itrip = itrip + 1
        sJ = h*t; sK = h*r
        jt = 1; if (sJ == -1) jt = 2
        if (sK == 1) then
          nkt = 2; ktl = [1, 2]
        else
          nkt = 1; ktl = [3, 3]
        end if
        nprim = 0
        do ik = 1, nkt
          kt = ktl(ik)
          call make_jk(nmult, jt, kt, Jm, Km, ok)
          jname = 'I'; if (jt == 2) jname = 'j'
          kname = 'I'; if (kt == 2) kname = 'mix'; if (kt == 3) kname = 'j'
          if (.not. ok) then
            write(*,'(A,3I3,A,A,A,A,A)') '  triple', h, t, r, '  J=', trim(jname), ' K=', trim(kname), &
                  ' : not admissible at this n (needs even n, or J K /= K J)'
            cycle
          end if
          call build_WR(Jm, Km, Wm, Rr)
          do ip = 1, 2
            pv = 1.0_dp; if (ip == 2) pv = -1.0_dp
            call constraint_nullspace(Wm, Rr, pv, N, nnull)
            nprim = 0
            if (nnull > 0) then
              do isamp = 1, nsamp
                x = matmul(N(:,1:nnull), matmul(transpose(N(:,1:nnull)), zall(isamp,:)))
                call build_tensor(x, A)
                call classify(A, xi, Jl1, Jl2, nper, herm, fpr)
                if (is_normal(Jl1, nper, herm, fpr)) then
                  nprim = nprim + 1
                  if (xi < bestxi(itrip)) then
                    bestxi(itrip) = xi; bestsamp(itrip) = isamp; bestJ2(itrip) = Jl2
                    bestjt(itrip) = jt; bestkt(itrip) = kt; bestip(itrip) = ip
                  end if
                end if
              end do
            end if
            write(*,'(A,3I3,A,A,A,A,A,I3,A,I5,A,I4,A,I4)') '  triple', h, t, r, '  J=', trim(jname), ' K=', trim(kname), &
                  '  p =', nint(pv), '  solution dim', nnull, '  normal samples', nprim, ' /', nsamp
            deallocate(N)
          end do
        end do
      end do
    end do
    deallocate(Wm, Rr, A, x)
  end do

  write(*,*)
  write(*,'(A)') ' ===== lowest correlation length per triple (normal samples only) ====='
  itrip = 0
  do h = -1, 1, 2
    call pick_irreps(s, h, ja, jb)
    call init_family(s, ja, jb, nmult)
    do t = -1, 1, 2
      do r = -1, 1, 2
        itrip = itrip + 1
        if (bestsamp(itrip) == 0) then
          write(*,'(A,3I3,A)') ' (sH,sT,sR) =', h, t, r, &
               '   no normal representative found at this n (more samples, or the ansatz needs larger n)'
          cycle
        end if
        jt = bestjt(itrip); kt = bestkt(itrip); ip = bestip(itrip)
        pv = 1.0_dp; if (ip == 2) pv = -1.0_dp
        jname = 'I'; if (jt == 2) jname = 'j'
        kname = 'I'; if (kt == 2) kname = 'mix'; if (kt == 3) kname = 'j'
        L = -2
        if (docheck == 1) then
          call make_jk(nmult, jt, kt, Jm, Km, ok)
          allocate(Wm(fam_dim,fam_dim), Rr(fam_dim,fam_dim), A(fam_dim,fam_dim,fam_nph), x(NPAR))
          call build_WR(Jm, Km, Wm, Rr)
          call constraint_nullspace(Wm, Rr, pv, N, nnull)
          x = matmul(N(:,1:nnull), matmul(transpose(N(:,1:nnull)), zall(bestsamp(itrip),:)))
          call build_tensor(x, A)
          L = injectivity_length(A, 12)
          deallocate(Wm, Rr, A, x, N)
        end if
        write(*,'(A,3I3,A,F11.4,A,A,A,A,A,I3,A,I5,A,I2)',advance='no') ' (sH,sT,sR) =', h, t, r, '   min xi =', bestxi(itrip), &
             '   J=', trim(jname), ' K=', trim(kname), ' p =', nint(pv), '   sample', bestsamp(itrip), &
             '   lambda2 in J=', bestJ2(itrip)
        if (docheck == 1) then
          write(*,'(A,I3)') '   numerical injectivity length =', L
        else
          write(*,*)
        end if
        write(*,'(A)') '   raw integer data (blocks aa, ab, ba, bb; re | im rows) -- project with the same constraints to reproduce:'
        do k = 1, 4
          do i = 1, nmult
            write(*,'(A,I1,A,*(I3))') '    blk', k, '  re', (nint(zall(bestsamp(itrip), 2*nmult*nmult*(k-1)+(i-1)*nmult+j)), j=1,nmult)
            write(*,'(A,I1,A,*(I3))') '    blk', k, '  im', (nint(zall(bestsamp(itrip), 2*nmult*nmult*(k-1)+nmult*nmult+(i-1)*nmult+j)), j=1,nmult)
          end do
        end do
      end do
    end do
  end do
end program sample_phases_sn
