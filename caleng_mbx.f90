PROGRAM test
  IMPLICIT NONE

  INTERFACE 
      SUBROUTINE mbxeng(coord,Vpot)
          DOUBLE PRECISION, INTENT(OUT) :: Vpot
          DOUBLE PRECISION, DIMENSION(:):: coord
      END SUBROUTINE mbxeng
      SUBROUTINE matpre1(Eulang,rotmat)
          DOUBLE PRECISION, DIMENSION(3), INTENT(IN)::Eulang
          DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT)::rotmat
          DOUBLE PRECISION :: phi, theta, chi, cp, sp, ct, st, ck, sk
      END SUBROUTINE matpre1
      SUBROUTINE rottrn1(rotmat,rwf,rsf,rcom)
          DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rwf,rcom
          DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: rotmat
          DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: rsf
          INTEGER :: i,j
      END SUBROUTINE rottrn1
  END INTERFACE 

  INTEGER :: nmon = 2, tdim, n_at
  DOUBLE PRECISION, ALLOCATABLE :: coord(:)
  CHARACTER(LEN=5), ALLOCATABLE :: at_name(:)
  DOUBLE PRECISION :: Vpot

  DOUBLE PRECISION, DIMENSION(3) :: ROwf
  DOUBLE PRECISION, DIMENSION(3) :: RH1wf
  DOUBLE PRECISION, DIMENSION(3) :: RH2wf
  DOUBLE PRECISION, DIMENSION(3) :: RO_sf
  DOUBLE PRECISION, DIMENSION(3) :: RH1_sf
  DOUBLE PRECISION, DIMENSION(3) :: RH2_sf
  DOUBLE PRECISION, DIMENSION(3) :: Eulang1
  DOUBLE PRECISION, DIMENSION(3) :: com1
  DOUBLE PRECISION, DIMENSION(3,3) :: rotmat
  INTEGER :: i, j, j1, ii, jj, ndim
  INTEGER :: iz, nz
  DOUBLE PRECISION :: zmin, zmax, dz, zz
  CHARACTER(LEN=5), ALLOCATABLE :: monomers(:)
  INTEGER, ALLOCATABLE :: nats(:)
  CHARACTER(LEN=20) :: xyz, json_file

  DOUBLE PRECISION, ALLOCATABLE :: com(:)
  DOUBLE PRECISION, ALLOCATABLE :: Eulang(:)

  DOUBLE PRECISION, PARAMETER :: PI = 4.0D0*ATAN(1.0D0)

  EXTERNAL initialize_system
  EXTERNAL finalize_system

  ndim=3
  n_at = nmon*ndim

  ALLOCATE(nats(nmon),monomers(nmon))

   do i=1,nmon
      nats(i) = 3
      monomers(i) = "h2o"
   enddo

  ! Set json file
  json_file = 'mbx_gas_phase.json'//CHAR(0)

  ALLOCATE(at_name(nmon*3))
  at_name(1) = "O"
  at_name(2) = "H"
  at_name(3) = "H"
  at_name(4) = "O"
  at_name(5) = "H"
  at_name(6) = "H"

  n_at = nmon*3
  do i =1,n_at
    at_name(i)=trim(at_name(i))//CHAR(0)
  enddo
  do i=1,nmon
    monomers(i) = trim(monomers(i))//CHAR(0)
  enddo

  ROwf(1)=0.0
  ROwf(2)=0.0
  ROwf(3)=0.06562        !6.5620209163576804E-002
  RH1wf(1)=-0.757395038  !-0.75739503834488564
  RH1wf(2)=0.0
  RH1wf(3)=-0.5207331156 !-0.52073311562494107
  RH2wf(1)=-RH1wf(1)
  RH2wf(2)=RH1wf(2)
  RH2wf(3)=RH1wf(3)

  tdim = nmon*ndim*3
  ALLOCATE(com(nmon*ndim),Eulang(nmon*ndim))
  ALLOCATE(coord(tdim))

  Eulang(1)=0.0
  Eulang(2)=0.0
  Eulang(3)=0.0
  Eulang(4)=0.0
  Eulang(5)=PI/4.0
  Eulang(6)=0.0

  CALL initialize_system(coord, nats, at_name, monomers, nmon, json_file)
  zmin = 2.5
  zmax = 10.0
  dz = 0.1
  nz = (zmax-zmin-0.5*dz)/0.1

  do iz = 1, nz
    zz = zmin+(iz-1)*dz
      com(1)=0.0
      com(2)=0.0
      com(3)=0.0
      com(4)=0.0
      com(5)=0.0
      com(6)=zz

      do i=1,nmon
        do j=1,ndim
            ii=j+(i-1)*ndim
            com1(j)=com(ii)
            Eulang1(j)=Eulang(ii)
        enddo

        call matpre1(Eulang1, rotmat)
        do j=1,3
            RO_sf(j)=0.d0
        enddo
        call rottrn1(rotmat, ROwf, RO_sf, com1)
        do j=1,3
            jj=j+(i-1)*9
            coord(jj)=RO_sf(j)
        enddo

       do j=1,3
           RH1_sf(j)=0.d0
        enddo
        call rottrn1(rotmat, RH1wf, RH1_sf, com1)
        do j=4,6
            jj=j+(i-1)*9
            j1=j-3
            coord(jj)=RH1_sf(j1)
        enddo

        do j=1,3
            RH2_sf(j)=0.d0
        enddo
        call rottrn1(rotmat, RH2wf, RH2_sf, com1)
        do j=7,9
            jj=j+(i-1)*9
            j1=j-6
            coord(jj)=RH2_sf(j1)
        enddo
      enddo

      CALL mbxeng(coord,Vpot)
      WRITE(7,*) zz, "    ", Vpot
  enddo

  CALL finalize_system()
  DEALLOCATE(coord)
  DEALLOCATE(at_name)
  DEALLOCATE(com, Eulang)

END PROGRAM 

SUBROUTINE mbxeng(coord,Vpot)
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: Vpot
  DOUBLE PRECISION, DIMENSION(:):: coord

  INTEGER :: n_at,i,j,k,nmon

  EXTERNAL get_energy

   i=size(coord)
   nmon = i/(3*3)

  n_at = nmon*3

  ! Energy call no gradients no pbc
  CALL get_energy(coord, n_at, Vpot)

END SUBROUTINE mbxeng

SUBROUTINE matpre1(Eulang,rotmat)
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)::Eulang
    DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT)::rotmat

    DOUBLE PRECISION :: phi, theta, chi, cp, sp, ct, st, ck, sk

    phi=Eulang(1)
    theta=Eulang(2)
    chi=Eulang(3)

    cp=cos(phi)
    sp=sin(phi)
    ct=cos(theta)
    st=sin(theta)
    ck=cos(chi)
    sk=sin(chi)

    rotmat(1,1)=cp*ct*ck-sp*sk
    rotmat(1,2)=-cp*ct*sk-sp*ck
    rotmat(1,3)=cp*st
    rotmat(2,1)=sp*ct*ck+cp*sk
    rotmat(2,2)=-sp*ct*sk+cp*ck
    rotmat(2,3)=sp*st
    rotmat(3,1)=-st*ck
    rotmat(3,2)=st*sk
    rotmat(3,3)=ct

END SUBROUTINE matpre1

SUBROUTINE rottrn1(rotmat,rwf,rsf,rcom)
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rwf,rcom
    DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: rotmat
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: rsf
    INTEGER :: i,j

    DO i=1,3
      rsf(i)=rcom(i)
      DO j=1,3
        rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
      ENDDO
   ENDDO

END SUBROUTINE rottrn1
