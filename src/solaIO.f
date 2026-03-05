c     Copyright Rasmus Munk Larsen, Stanford University, 2003
c
c     Fortran I/O routines for reading and writing binary data files.
c     Replaces the former C implementation in byteswapIO.c.
c

      subroutine readmesh(iunit, filename, rthetw, ldr, N_points)
c
c     Read mesh file (r, theta, weights) written by set-2drls.
c
      implicit none
      integer iunit, ldr, N_points, j
      character*(*) filename
      double precision rthetw(ldr, *)

      open(iunit, file=filename, form='unformatted', status='old')
      read(iunit) N_points
      read(iunit) (rthetw(j,1),j=1,N_points),
     c     (rthetw(j,2),j=1,N_points),
     c     (rthetw(j,3),j=1,N_points)
      close(iunit)
      return
      end


      subroutine readfg(iunit, filename, f1, f2, ldf, g1, g2, ldg,
     c     M_kers, N_points, N_rad, M_nl, N_theta, M_lm, inl, ilm)
c
c     Read kernel matrices (f1, f2, g1, g2) written by set-2drls.
c     The matrices are stored transposed: set-2drls writes f1(mode,rad)
c     and we read into f1(rad,mode).
c
      implicit none
      integer iunit, ldf, ldg, i, j, idummy
      integer M_kers, N_points, N_rad, M_nl, N_theta, M_lm
      character*(*) filename
      double precision f1(ldf, *), f2(ldf, *)
      double precision g1(ldg, *), g2(ldg, *)
      integer inl(*), ilm(*)

      open(iunit, file=filename, form='unformatted', status='old')
c     Record 1: dimensions
      read(iunit) M_kers, N_points
c     Record 2: mode indices (first value is nnlm again, skip it)
      read(iunit) idummy, (inl(i),i=1,M_kers), (ilm(i),i=1,M_kers)
c     Record 3: radial dimensions
      read(iunit) N_rad, M_nl
c     Records 4-5: f1, f2 (transposed read)
      read(iunit) ((f1(j,i), j=1,N_rad), i=1,M_nl)
      read(iunit) ((f2(j,i), j=1,N_rad), i=1,M_nl)
c     Record 6: latitudinal dimensions
      read(iunit) N_theta, M_lm
c     Records 7-8: g1, g2 (transposed read)
      read(iunit) ((g1(j,i), j=1,N_theta), i=1,M_lm)
      read(iunit) ((g2(j,i), j=1,N_theta), i=1,M_lm)
      close(iunit)
      return
      end


      subroutine rdamdl(iunit, filename, x, aa, data, nn, nmod,
     c     ivar, icry, ia)
c
c     Read ADIPLS solar model file (amdl format).
c     Uses stream access since the record has mixed types and
c     the number of elements depends on data read from the record.
c
      implicit none
      integer iunit, nn, nmod, ivar, icry, ia
      character*(*) filename
      double precision x(*), aa(ia, *), data(*)
      integer i, j, recmarker
      double precision d8

      open(iunit, file=filename, access='stream', status='old')
      read(iunit) recmarker
      read(iunit) nmod
      read(iunit) nn
      read(iunit) (data(i), i=1,8)
      d8 = data(8) + 0.1d0
      if (d8.ge.100d0) then
         ivar = 8
      else if (d8.ge.10d0) then
         ivar = 6
      else
         ivar = 5
      endif
      do i=1,nn
         read(iunit) x(i)
         read(iunit) (aa(j,i), j=1,ivar)
      enddo
      read(iunit) recmarker
      close(iunit)
      icry = 0
      return
      end


      subroutine rdmatr(iunit, filename, m, n, a, lda)
c
c     Read a double precision matrix from an unformatted file.
c
      implicit none
      integer iunit, m, n, lda, i, j
      character*(*) filename
      double precision a(lda, *)

      open(iunit, file=filename, form='unformatted', status='old')
      read(iunit) m, n
      read(iunit) ((a(i,j), i=1,m), j=1,n)
      close(iunit)
      return
      end


      subroutine wrmatr(iunit, filename, m, n, a, lda)
c
c     Write a double precision matrix to an unformatted file.
c
      implicit none
      integer iunit, m, n, lda, i, j
      character*(*) filename
      double precision a(lda, *)

      open(iunit, file=filename, form='unformatted', status='unknown')
      write(iunit) m, n
      write(iunit) ((a(i,j), i=1,m), j=1,n)
      close(iunit)
      return
      end


      subroutine writebidiag(iunit, filename, n_points, m_kers,
     c     n_iter, bidiag, ldbidiag, hhvec, ldhhvec,
     c     U, ldu, V, ldv)
c
c     Write Lanczos bidiagonalization data to file.
c
      implicit none
      integer iunit, n_points, m_kers, n_iter, i, j
      integer ldbidiag, ldhhvec, ldu, ldv
      character*(*) filename
      double precision bidiag(ldbidiag, *), hhvec(*)
      double precision U(ldu, *), V(ldv, *)

      open(iunit, file=filename, form='unformatted', status='unknown')
      write(iunit) n_points, m_kers, n_iter
      write(iunit) ((bidiag(i,j), i=1,n_iter), j=1,2)
      write(iunit) (hhvec(i), i=1,m_kers+2)
      write(iunit) ((U(i,j), i=1,n_points), j=1,n_iter)
      write(iunit) ((V(i,j), i=1,m_kers-1), j=1,n_iter)
      close(iunit)
      return
      end


      subroutine readbidiag(iunit, filename, n_points, m_kers,
     c     n_iter, bidiag, ldbidiag, hhvec, ldhhvec,
     c     U, ldu, V, ldv)
c
c     Read Lanczos bidiagonalization data from file.
c
      implicit none
      integer iunit, n_points, m_kers, n_iter, i, j
      integer ldbidiag, ldhhvec, ldu, ldv
      character*(*) filename
      double precision bidiag(ldbidiag, *), hhvec(*)
      double precision U(ldu, *), V(ldv, *)

      open(iunit, file=filename, form='unformatted', status='old')
      read(iunit) n_points, m_kers, n_iter
      read(iunit) ((bidiag(i,j), i=1,n_iter), j=1,2)
      read(iunit) (hhvec(i), i=1,m_kers+2)
      read(iunit) ((U(i,j), i=1,n_points), j=1,n_iter)
      read(iunit) ((V(i,j), i=1,m_kers-1), j=1,n_iter)
      close(iunit)
      return
      end


      subroutine writeavker(iunit, filename, n_targets, n_points,
     c     icase, trdoff, targetparms, ldtp, rot, sigma_rot,
     c     depart, avker, ldavker)
c
c     Write averaging kernels to file.
c
      implicit none
      integer iunit, n_targets, n_points, icase, ldtp, ldavker
      integer i, j
      character*(*) filename
      double precision trdoff(*), targetparms(ldtp, *)
      double precision rot(*), sigma_rot(*), depart(*)
      double precision avker(ldavker, *)
      double precision zero
      parameter(zero = 0d0)

      open(iunit, file=filename, form='unformatted', status='unknown')
      do i=1,n_targets
         write(iunit) icase, trdoff(i),
     c        targetparms(i,1), targetparms(i,2),
     c        targetparms(i,3), targetparms(i,4),
     c        rot(i), sigma_rot(i), depart(i), zero, zero,
     c        n_points, (avker(j,i), j=1,n_points)
      enddo
      close(iunit)
      return
      end


      subroutine appendavker(iunit, filename, n_targets, n_points,
     c     icase, trdoff, targetparms, ldtp, rot, sigma_rot,
     c     depart, avker, ldavker)
c
c     Append averaging kernels to file.
c
      implicit none
      integer iunit, n_targets, n_points, icase, ldtp, ldavker
      integer i, j
      character*(*) filename
      double precision trdoff(*), targetparms(ldtp, *)
      double precision rot(*), sigma_rot(*), depart(*)
      double precision avker(ldavker, *)
      double precision zero
      parameter(zero = 0d0)

      open(iunit, file=filename, form='unformatted', status='old',
     c     position='append')
      do i=1,n_targets
         write(iunit) icase, trdoff(i),
     c        targetparms(i,1), targetparms(i,2),
     c        targetparms(i,3), targetparms(i,4),
     c        rot(i), sigma_rot(i), depart(i), zero, zero,
     c        n_points, (avker(j,i), j=1,n_points)
      enddo
      close(iunit)
      return
      end
