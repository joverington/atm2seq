       program ATM2SEQ
C
C copyright (c) 1988-1991 John Overington
C part of joy
C
C takes a PDB format file and then converts it to one letter code, PIR
C format, in blocks of 75
C
C changed so took input from STDIN, and OUTPUT to STDOUT
C
C problem taking form STDIN is that there is no way of getting file name
C and thus code from stream. Therefore will for now, only use file input
C for any of my stuff.
C
C 23-12-92
C
C Added options -P and -F for PIR and FASTA format resepectively
C       default is PIR format
C
C 14-06-94
C
C Fixed problem with using just CAs to get sequence, will now look for
C differences in NUMRES field.
C
C 5-8-95
C
C Added output of name and species fields. added chack of CA-CA distances
C 4A distance cutoff
C
C 4-9-97
C
C Added correct handling of new species and protein fields
C
C 18-9-97
C
C Improved handling of RFACTOR fields
C
C 30-4-98
C
C Increased size of MAXRES (to cope with 1bgl)
C
C 28-1-00
C
C improved handling and assignment of EXPDTA fields
C
      integer MAXRES, MAXATS, MAXFILELEN
      integer STDIN, STDOUT, STDERR
      parameter (MAXRES=10000, MAXFILELEN=256)
      parameter (STDIN=5, STDOUT=6, STDERR=0)
C
      character*(MAXRES) CARD
      character*(MAXFILELEN) FILE
      character*60 NAME, SOURCE
      character*1 FLAG
      character*5 RESNUM(MAXRES)
      character*3 RESNAM(MAXRES)
      character*4 CODE
      character*1 CHNNAM(MAXRES), SEQ(MAXRES), TSEQ(MAXRES)
      character*1 TOLOWER
      integer NP, NPS
      integer IIN
      integer IARGC, IIO, I, NUMCHN, NTIMES, N, K, NUMRES, LASTCHAR
      integer NUMRESO
      logical PIR, BREAK(MAXRES), DOBREAK
      real RESOL, RFACT
      data PIR /.true./
      data DOBREAK /.true./
C
C check for options
C
      NOPT=0
      if (IARGC() .gt. 0) then
        call GETARG(1,FILE)
        if (FILE(1:1) .eq. '-') then
          NOPT=1
          if      (index(FILE,'F') .gt. 0) then
            PIR=.false.
          else if (index(FILE,'P') .gt. 0) then
            PIR=.true.
          else if (index(FILE,'B') .gt. 0) then
            DOBREAK=.false.
          else
            write (STDERR,'(''atm2seq: unknown option '',A)') FILE(1:10)
          end if
        end if
      end if
C
      if (IARGC()-NOPT .gt. 0) then
        IIN=1
        IIO=STDOUT
        if (NOPT .eq. 1) then
          call GETARG(2,FILE)
        else
          call GETARG(1,FILE)
        end if
        open (file=FILE,unit=IIN,status='old',form='formatted',err=3)
      else
        IIN=5
        IIO=STDOUT
        do I=1,MAXLENFILE
          FILE(I:I)=' '
        end do
      end if
C
      do I=1,MAXRES
        CARD(I:I)=' '
      end do
C
C this is where the problem starts, have no way to sensibly deal with
C residues where CA is missing (it is quite common to get N alone in terminal
C residues
C Can make a number of assumptions though, will have no alternative locations
C since we will have been through stripper
C
      call READCA(IIN,RESNUM,RESNAM,CHNNAM,NUMRES,NUMCHN,NAME,
     +            SOURCE,FLAG,RESOL,RFACT,BREAK)
C
C coPy the original number of resiudes to NUMRESO
C
      NUMRESO=NUMRES
C
      call THRONE(RESNAM,SEQ,NUMRES)
      SEQ(NUMRES+1)='*'
C
C fill in chain breaks with / characters (vars Prefixed with T)
C
      if (DOBREAK) then
        M=0
        do I=1,NUMRES
          M=M+1
          TSEQ(M)=SEQ(I)
          if (BREAK(I)) then
            M=M+1
            TSEQ(M)='/' 
            M=M+1
            TSEQ(M)=CHAR(10)
          end if
        end do
        NUMRES=M
        do I=1,NUMRES
          SEQ(I)=TSEQ(I)
        end do
        SEQ(NUMRES+1)='*'
      end if
C
      if (MOD(NUMRES,75).eq.0) then
        NTIMES=(NUMRES/75)
      else
        NTIMES=(NUMRES/75)+1
      end if
C
C miss out extension of file (if it has one)
C
C catch bug for no . in file name
C
      NP=index(FILE,'.')-1
      if (NP .lt. 1) then
        NP=lastchar(FILE)
      end if
      NN=lastchar(NAME)
      do I=1,NN
        NAME(I:I)=tolower(NAME(I:I))
      end do
      NS=lastchar(SOURCE)
      do I=1,NS
        SOURCE(I:I)=tolower(SOURCE(I:I))
      end do
C
C handle special case of pdb preceeded entries....
C
      if (FILE(1:3) .eq. 'pdb') then
        CODE(1:NP-3)=FILE(4:NP)
        NPS=NP-3
      else
        CODE(1:NP)=FILE(1:NP)
        NPS=NP
      end if
C
      if (PIR) then
        write (IIO,'(''>P1;'',A,/,''structure'',A,'':'',A,
     +               '':'',A,'':'',A,'':'',A,'':'',A,'':'',A,
     +               '':'',A,'':'',F5.2,'':'',F5.2)')
     +         FILE(1:NP),FLAG,CODE(1:NPS),RESNUM(1),CHNNAM(1),RESNUM(NUMRESO),
     +         CHNNAM(NUMRESO),NAME(1:NN),SOURCE(1:NS),
     +         RESOL,RFACT
      else
        write (IIO,'(''>'',A,'' structure'',A,'':'',A,
     +               '':'',A,'':'',A,'':'',A,'':'',A,'':'',A,
     +               '':'',A,'':'',F5.2,'':'',F5.2)')
     +         FILE(1:NP),FLAG,CODE(1:NPS),RESNUM(1),CHNNAM(1),RESNUM(NUMRESO),
     +         CHNNAM(NUMRESO),NAME(1:NN),SOURCE(1:NS),
     +         RESOL,RFACT
      end if
C
      N=1
      do I=1,NTIMES-1
        write (IIO,'(75A)') (SEQ(K),K=N,N+74)
        N=N+75
      end do
C
C if FASTA then omit last *
C
      if (.not. PIR) then
        NUMRES=NUMRES-1
      end if
C
      write (IIO,'(75A)') (SEQ(K),K=N,NUMRES+1)
C
      if (IIN .ne. STDIN) then
        close (UNIT=IIN)
      end if
C
      call exit(0)
C
3     write (STDERR,'(''atm2seq: cannot open '',A)') FILE
C
      end
C
C -----------------------------------------------------------------------------
C
      subroutine THRONE(CHAR,CHARNW,NRES)
C
      integer MAXRES, I, J, NRES
      parameter (MAXRES=10000)
      character*(*) CHAR(MAXRES),CHARNW(MAXRES)
      character*24 ACIDS1
      character*3 ACIDS3(24)
      integer STDERR
      logical FOUND
C
      data STDERR /0/
C
      DATA ACIDS1/'ARNDCQEGHILKMFPSTWYVBZXC'/
      DATA ACIDS3/'ALA','ARG','ASN','ASP','CYS',
     -            'GLN','GLU','GLY','HIS','ILE',
     -            'LEU','LYS','MET','PHE','PRO',
     -            'SER','THR','TRP','TYR','VAL',
     -            'ASX','GLX','UNK','CYH'/
C
C Convert three letter code to one letter code
C
      do I=1,NRES
        FOUND=.false.
        do J=1,23
          if (CHAR(I) .eq. ACIDS3(J)) then
            CHARNW(I)(1:1)=ACIDS1(J:J)
            FOUND=.true.
            goto 1
          end if
        end do
1       continue
        if (.NOT.FOUND) THEN
          write (STDERR,'(''atm2seq: unknown residue code '',A)')
     +           CHAR(I)
          CHARNW(I)='X'
        end if
      end do
C
      return
      end
C
C -----------------------------------------------------------------------------
C
      subroutine READCA(IIN,RESNUM,RESNAM,CHNNAM,NRES,NCHN,NAME,
     +                  SOURCE,FLAG,RESOL,RFACT,BREAK)
C
C subroutine that reads a PDB format file and returns data on a per residue
C basis
C
      integer MAXRES
      parameter (MAXRES=10000)
C
      character*66 CARD(MAXRES), BUFFER
      character*60 NAME, SOURCE
      character*5 RESNUM(MAXRES), OLDNUM
      character*3 RESNAM(MAXRES)
      character*1 CHNNAM(MAXRES), FLAG
      integer NUMRES, IIN, STDERR, STDOUT, STDIN, I
      integer NCHN, NREC, J, M, K, NATS
      real RESOL, RFACT, COORD(3,MAXRES), CUTOFF
      logical BREAK(MAXRES)
C
      data STDIN  /5/
      data STDOUT /6/
      data STDERR /0/
C
C CUTOFF is break cutoff ** 2
C
      CUTOFF=4.5
      CUTOFF=CUTOFF * CUTOFF
C
C make the assumption that most records are X-ray (but signify not sure
C unless EXPDTA tell us otherwise.
C
      FLAG='x'
      NAME=  'undefined                                                   '
      SOURCE='undefined                                                   '
      RFACT=99.99
      RESOL=9.99
C
      NRES=0
      NCHN=0
C
C read all data into CARD buffer, NREC items
C make sure key step point is new subfield (CARD(23:27))
C
      NREC=0
      OLDNUM='@@@@@'
4     read (IIN,'(A)',END=3) BUFFER
        if (BUFFER(1:6) .eq. 'ATOM  ') THEN
          if (BUFFER(23:27) .ne. OLDNUM) then
            NREC=NREC+1
            CARD(NREC)=BUFFER
            OLDNUM=BUFFER(23:27)
            read (BUFFER,'(30X,3F8.3)',err=914) (COORD(K,NREC),K=1,3)
            if (NREC .gt. (MAXRES)) then
              write (STDERR,'(''atm2seq: too many records'')')
              call EXIT(1)
            end if
          end if
        else if (BUFFER(1:6) .eq. 'COMPND') then
          if (BUFFER(10:10) .eq. ' ' .and. BUFFER(11:17) .ne. 'MOL_ID') then
            read (BUFFER,'(10X,A)') NAME
          end if
          if (BUFFER(12:20) .eq. 'MOLECULE:') then
            read (BUFFER,'(21X,A)') NAME
          end if
        else if (BUFFER(1:6) .eq. 'SOURCE') then
          if (BUFFER(10:10) .eq. ' ' .and. BUFFER(11:17) .ne. 'MOL_ID') then
            read (BUFFER,'(10X,A)') SOURCE
          end if
          if (BUFFER(12:31) .eq. 'ORGANISM_SCIENTIFIC:') then
            read (BUFFER,'(32X,A)') SOURCE
          end if
        else if (BUFFER(1:13) .eq. 'EXPDTA    X-R') then
          FLAG='X'
        else if (BUFFER(1:13) .eq. 'EXPDTA    NMR') then
          FLAG='N'
          RESOL=-1.00
          RFACT=-1.00
        else if (BUFFER(1:13) .eq. 'EXPDTA    THE') then
          FLAG='M'
          RESOL=-1.00
          RFACT=-1.00
        else if (BUFFER(1:6) .eq. 'REMARK') then
          if (BUFFER(12:22) .eq. 'RESOLUTION.') then
            read (BUFFER(23:27),fmt='(F5.2)',err=911) RESOL
            FLAG='X'
            goto 913
911         RESOL=9.99
913         continue
          else if (BUFFER(14:47) .eq. 'R VALUE            (WORKING SET) :') then
            read (BUFFER(49:53),fmt='(F5.3)',err=921) RFACT
            FLAG='X'
            RFACT=RFACT*100.0
            goto 923
          else if (BUFFER(14:38) .eq. 'R VALUE                  ') then
            read (BUFFER(41:45),fmt='(F5.3)',err=921) RFACT
            FLAG='X'
            RFACT=RFACT*100.0
            goto 923
921         RFACT=99.9
923         continue
          end if
        else if (BUFFER(1:6) .eq. 'ENDMDL') then
          goto 3
        end if
      GOTO 4
3     continue
cjpo      NREC=NREC-1
C
C Process all these records
C
      do M=1,NREC
        BREAK(M)=.false.
        if (M.gt.1) then
          if (CARD(M)(22:22).NE.CARD(M-1)(22:22)) NCHN=NCHN+1
        end if
        read (CARD(M),'(13X,3X,1X,A3,1X,A1,A5)') RESNAM(M),CHNNAM(M),RESNUM(M)
C
C The chain break algorithm is not really that smart (since there is no
C guarantee that the CAs of a residue will be chosen (it is simPly the first
C which by convention will usually be N
C
        if (M .lt. NREC-1) then
          if (DIST2(COORD(1,M),COORD(1,M+1)) .gt. CUTOFF) then
            BREAK(M)=.true.
          end if
        end if
      end do
C
      NRES=NREC
      NCHN=NCHN+1
C
C remove PDB record terminators
C
      do I=1,60
        if (NAME(I:I) .eq. ';') then
          NAME(I:I)=' '
        end if
        if (SOURCE(I:I) .eq. ';') then
          SOURCE(I:I)=' '
        end if
      end do
C
      return
C
914   write (STDERR,'(''atm2seq: bad read of coordinates'')')
      call EXIT(1)
      end
