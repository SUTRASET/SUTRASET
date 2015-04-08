C     ITERATIVE MATRIX SOLVERS (SLAP) ............ SUTRA VERSION 2.2     SLAP...........100
C                                                                        SLAP...........200
C *** PURPOSE :                                                          SLAP...........300
C ***  TO SOLVE MATRIX PROBLEMS USING ITERATIVE METHODS.  THE COMPUTER   SLAP...........400
C ***  CODE IN THIS FILE IS FROM THE SLATEC SOFTWARE LIBRARY.  IT        SLAP...........500
C ***  INCLUDES A SUBSET OF THE SLAP SOLVER PACKAGE AND ADDITIONAL CODE  SLAP...........600
C ***  UPON WHICH SLAP DEPENDS.  CREDITS ARE LISTED AT THE BEGINNING     SLAP...........700
C ***  OF EACH SUBROUTINE.  MODIFICATIONS TO THE ORIGINAL CODE MADE BY   SLAP...........800
C ***  THE AUTHORS OF SUTRA ARE MARKED AS SUCH; SEARCH FOR THE KEYWORD   SLAP...........900
C ***  "SUTRA".  GENERAL DOCUMENTATION FOR THE SLAP SOLVERS APPEARS IN   SLAP..........1000
C ***  SUBROUTINE DLPDOC BELOW.                                          SLAP..........1100
C ***                                                                    SLAP..........1200
C ***  THE FOLLOWING DISCLAIMER APPEARS IN FILE "aaaaaa.f" OF THE        SLAP..........1300
C ***  SLATEC LIBRARY v4.1:                                              SLAP..........1400
C ***                                                                    SLAP..........1500
C =====================================================================  SLAP..........1600
C                                                                        SLAP..........1700
C   "The SLATEC Common Mathematical Library is issued by the following   SLAP..........1800
C                                                                        SLAP..........1900
C           Air Force Weapons Laboratory, Albuquerque                    SLAP..........2000
C           Lawrence Livermore National Laboratory, Livermore            SLAP..........2100
C           Los Alamos National Laboratory, Los Alamos                   SLAP..........2200
C           National Institute of Standards and Technology, Washington   SLAP..........2300
C           National Energy Research Supercomputer Center, Livermore     SLAP..........2400
C           Oak Ridge National Laboratory, Oak Ridge                     SLAP..........2500
C           Sandia National Laboratories, Albuquerque                    SLAP..........2600
C           Sandia National Laboratories, Livermore                      SLAP..........2700
C                                                                        SLAP..........2800
C   All questions concerning the distribution of the library should be   SLAP..........2900
C   directed to the NATIONAL ENERGY SOFTWARE CENTER, 9700 Cass Ave.,     SLAP..........3000
C   Argonne, Illinois  60439, and not to the authors of the subprograms. SLAP..........3100
C                                                                        SLAP..........3200
C                    * * * * * Notice * * * * *                          SLAP..........3300
C                                                                        SLAP..........3400
C   This material was prepared as an account of work sponsored by the    SLAP..........3500
C   United States Government.  Neither the United States, nor the        SLAP..........3600
C   Department of Energy, nor the Department of Defense, nor any of      SLAP..........3700
C   their employees, nor any of their contractors, subcontractors, or    SLAP..........3800
C   their employees, makes any warranty, expressed or implied, or        SLAP..........3900
C   assumes any legal liability or responsibility for the accuracy,      SLAP..........4000
C   completeness, or usefulness of any information, apparatus, product,  SLAP..........4100
C   or process disclosed, or represents that its use would not infringe  SLAP..........4200
C   upon privately owned rights."                                        SLAP..........4300
C                                                                        SLAP..........4400
C =====================================================================  SLAP..........4500
C ***                                                                    SLAP..........4600
C ***  NOTE:  AS OF THIS WRITING, THE NATIONAL ENERGY SOFTWARE CENTER    SLAP..........4700
C ***  HAD RELOCATED TO OAK RIDGE NATIONAL LABORATORY AND HAD BEEN       SLAP..........4800
C ***  RENAMED THE ENERGY SCIENCE & TECHNOLOGY SOFTWARE CENTER (ESTSC).  SLAP..........4900
C ***                                                                    SLAP..........5000
C ***  BASED ON COMMUNICATIONS WITH ESTSC AND THE PRIMARY AUTHOR OF      SLAP..........5100
C ***  SLAP, THE AUTHORS OF SUTRA UNDERSTAND THAT THE CODE IN THIS FILE  SLAP..........5200
C ***  MAY BE FREELY DISTRIBUTED WITH SUTRA.  HOWEVER, NEITHER THE       SLAP..........5300
C ***  UNITED STATES GOVERNMENT, THE DEPARTMENT OF THE INTERIOR, THE     SLAP..........5400
C ***  U.S. GEOLOGICAL SURVEY, NOR THE AUTHORS OF SUTRA MAKE ANY         SLAP..........5500
C ***  WARRANTY, EXPRESSED OR IMPLIED, OR ASSUME ANY LEGAL LIABILITY OR  SLAP..........5600
C ***  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF   SLAP..........5700
C ***  THIS CODE OR ANY MODIFICATIONS MADE TO THIS CODE, OR GUARANTEE    SLAP..........5800
C ***  THAT ITS UNLIMITED USE WOULD NOT INFRINGE UPON PRIVATELY OWNED    SLAP..........5900
C ***  RIGHTS.                                                           SLAP..........6000
C ***                                                                    SLAP..........6100
C ***  THE CODE LISTING BEGINS BELOW WITH SUBROUTINE DLPDOC, WHICH       SLAP..........6200
C ***  PROVIDES GENERAL DOCUMENTATION FOR THE DOUBLE-PRECISION VERSION   SLAP..........6300
C ***  OF THE SLAP SOLVER PACKAGE.  THE REMAINING SUBROUTINES AND        SLAP..........6400
C ***  FUNCTIONS ARE THEN LISTED IN ALPHABETICAL ORDER.                  SLAP..........6500
C ***                                                                    SLAP..........6600
C =====================================================================  SLAP..........6700
C =====================================================================  SLAP..........6800
C                                                                        SLAP..........6900
*DECK DLPDOC                                                             SLAP..........7000
      SUBROUTINE DLPDOC                                                  SLAP..........7100
C***BEGIN PROLOGUE  DLPDOC                                               SLAP..........7200
C***PURPOSE  Sparse Linear Algebra Package Version 2.0.2 Documentation.  SLAP..........7300
C            Routines to solve large sparse symmetric and nonsymmetric   SLAP..........7400
C            positive definite linear systems, Ax = b, using precondi-   SLAP..........7500
C            tioned iterative methods.                                   SLAP..........7600
C***LIBRARY   SLATEC (SLAP)                                              SLAP..........7700
C***CATEGORY  D2A4, D2B4, Z                                              SLAP..........7800
C***TYPE      DOUBLE PRECISION (SLPDOC-S, DLPDOC-D)                      SLAP..........7900
C***KEYWORDS  BICONJUGATE GRADIENT SQUARED, DOCUMENTATION,               SLAP..........8000
C             GENERALIZED MINIMUM RESIDUAL, ITERATIVE IMPROVEMENT,       SLAP..........8100
C             NORMAL EQUATIONS, ORTHOMIN,                                SLAP..........8200
C             PRECONDITIONED CONJUGATE GRADIENT, SLAP,                   SLAP..........8300
C             SPARSE ITERATIVE METHODS                                   SLAP..........8400
C***AUTHOR  Seager, Mark. K., (LLNL)                                     SLAP..........8500
C             User Systems Division                                      SLAP..........8600
C             Lawrence Livermore National Laboratory                     SLAP..........8700
C             PO BOX 808, L-60                                           SLAP..........8800
C             Livermore, CA 94550                                        SLAP..........8900
C             (FTS) 543-3141, (510) 423-3141                             SLAP..........9000
C             seager@llnl.gov                                            SLAP..........9100
C***DESCRIPTION                                                          SLAP..........9200
C                                 The                                    SLAP..........9300
C                    Sparse Linear Algebra Package                       SLAP..........9400
C                      Double Precision Routines                         SLAP..........9500
C                                                                        SLAP..........9600
C                @@@@@@@  @            @@@    @@@@@@@@                   SLAP..........9700
C               @       @ @           @   @   @       @                  SLAP..........9800
C               @         @          @     @  @       @                  SLAP..........9900
C                @@@@@@@  @         @       @ @@@@@@@@                   SLAP.........10000
C                       @ @         @@@@@@@@@ @                          SLAP.........10100
C               @       @ @         @       @ @                          SLAP.........10200
C                @@@@@@@  @@@@@@@@@ @       @ @                          SLAP.........10300
C                                                                        SLAP.........10400
C      @       @                            @@@@@@@        @@@@@         SLAP.........10500
C      @       @                           @       @      @    @@        SLAP.........10600
C      @       @  @@@@@@@  @ @@                    @     @    @  @       SLAP.........10700
C      @       @ @       @ @@  @             @@@@@@      @   @   @       SLAP.........10800
C       @     @  @@@@@@@@@ @                @            @  @    @       SLAP.........10900
C        @   @   @         @               @         @@@  @@    @        SLAP.........11000
C         @@@     @@@@@@@  @               @@@@@@@@@ @@@   @@@@@         SLAP.........11100
C                                                                        SLAP.........11200
C                                                                        SLAP.........11300
C    =================================================================   SLAP.........11400
C    ========================== Introduction =========================   SLAP.........11500
C    =================================================================   SLAP.........11600
C      This package was  originally derived from a set of  iterative     SLAP.........11700
C      routines written by Anne Greenbaum, as announced in "Routines     SLAP.........11800
C      for Solving Large Sparse Linear Systems",  Tentacle, Lawrence     SLAP.........11900
C      Livermore  National  Laboratory,  Livermore  Computing Center     SLAP.........12000
C      (January 1986), pp 15-21.                                         SLAP.........12100
C                                                                        SLAP.........12200
C    This document  contains the specifications for  the  SLAP Version   SLAP.........12300
C    2.0 package, a Fortran 77  package  for  the  solution  of  large   SLAP.........12400
C    sparse   linear systems, Ax  =  b,  via  preconditioned iterative   SLAP.........12500
C    methods.   Included in  this  package are "core"  routines  to do   SLAP.........12600
C    Iterative   Refinement  (Jacobi's  method),  Conjugate  Gradient,   SLAP.........12700
C    Conjugate Gradient on the normal equations, AA'y = b,  (where x =   SLAP.........12800
C    A'y and  A' denotes the  transpose of   A), BiConjugate Gradient,   SLAP.........12900
C    BiConjugate  Gradient  Squared, Orthomin and  Generalized Minimum   SLAP.........13000
C    Residual Iteration.    These "core" routines   do  not  require a   SLAP.........13100
C    "fixed"   data  structure   for storing  the   matrix  A  and the   SLAP.........13200
C    preconditioning   matrix  M.   The  user  is free  to  choose any   SLAP.........13300
C    structure that facilitates  efficient solution  of the problem at   SLAP.........13400
C    hand.  The drawback  to this approach  is that the user must also   SLAP.........13500
C    supply at least two routines  (MATVEC and MSOLVE,  say).   MATVEC   SLAP.........13600
C    must calculate, y = Ax, given x and the user's data structure for   SLAP.........13700
C    A.  MSOLVE must solve,  r = Mz, for z (*NOT*  r) given r  and the   SLAP.........13800
C    user's data  structure for  M (or its  inverse).  The user should   SLAP.........13900
C    choose M so that  inv(M)*A  is approximately the identity and the   SLAP.........14000
C    solution step r = Mz is "easy" to  solve.  For some of the "core"   SLAP.........14100
C    routines (Orthomin,  BiConjugate Gradient and  Conjugate Gradient   SLAP.........14200
C    on the  normal equations)   the user must  also  supply  a matrix   SLAP.........14300
C    transpose times   vector  routine  (MTTVEC,  say)  and (possibly,   SLAP.........14400
C    depending    on the "core"  method)   a  routine  that solves the   SLAP.........14500
C    transpose  of   the   preconditioning    step     (MTSOLV,  say).   SLAP.........14600
C    Specifically, MTTVEC is a routine which calculates y = A'x, given   SLAP.........14700
C    x and the user's data structure for A (A' is the transpose of A).   SLAP.........14800
C    MTSOLV is a routine which solves the system r = M'z for z given r   SLAP.........14900
C    and the user's data structure for M.                                SLAP.........15000
C                                                                        SLAP.........15100
C    This process of writing the matrix vector operations  can be time   SLAP.........15200
C    consuming and error  prone.  To alleviate  these problems we have   SLAP.........15300
C    written drivers   for  the  "core" methods  that  assume the user   SLAP.........15400
C    supplies one of two specific data structures (SLAP Triad and SLAP   SLAP.........15500
C    Column format), see  below.  Utilizing these  data structures  we   SLAP.........15600
C    have augmented   each  "core" method  with   two preconditioners:   SLAP.........15700
C    Diagonal  Scaling and Incomplete Factorization.  Diagonal scaling   SLAP.........15800
C    is easy to implement, vectorizes very  well and for problems that   SLAP.........15900
C    are  not too  ill-conditioned  reduces the  number  of iterations   SLAP.........16000
C    enough   to warrant its use.  On   the other  hand, an Incomplete   SLAP.........16100
C    factorization  (Incomplete  Cholesky for  symmetric systems   and   SLAP.........16200
C    Incomplete LU for nonsymmetric  systems) may  take much longer to   SLAP.........16300
C    calculate, but it reduces the iteration count (for most problems)   SLAP.........16400
C    significantly.  Our implementations  of IC and ILU  vectorize for   SLAP.........16500
C    machines with hardware gather scatter, but the vector lengths can   SLAP.........16600
C    be quite short if  the  number  of non-zeros  in a column is  not   SLAP.........16700
C    large.                                                              SLAP.........16800
C                                                                        SLAP.........16900
C    =================================================================   SLAP.........17000
C    ==================== Supplied Data Structures ===================   SLAP.........17100
C    =================================================================   SLAP.........17200
C    The following describes the data   structures supplied  with  the   SLAP.........17300
C    package: SLAP Triad and Column formats.                             SLAP.........17400
C                                                                        SLAP.........17500
C    ====================== S L A P Triad format =====================   SLAP.........17600
C                                                                        SLAP.........17700
C    In the SLAP Triad format only the non-zeros are stored.  They may   SLAP.........17800
C    appear in *ANY* order.  The user supplies three  arrays of length   SLAP.........17900
C    NELT, where NELT  is the   number of  non-zeros  in the   matrix:   SLAP.........18000
C    (IA(NELT),  JA(NELT), A(NELT)).  If  the matrix is symmetric then   SLAP.........18100
C    one need only store the lower triangle (including  the  diagonal)   SLAP.........18200
C    and NELT would be the corresponding  number  of non-zeros stored.   SLAP.........18300
C    For each non-zero the user puts the row and column  index of that   SLAP.........18400
C    matrix  element   in the  IA  and JA  arrays.  The  value  of the   SLAP.........18500
C    non-zero matrix element is placed  in  the corresponding location   SLAP.........18600
C    of  the A array.   This  is an extremely  easy  data structure to   SLAP.........18700
C    generate.  On the other hand, it is not very  efficient on vector   SLAP.........18800
C    computers for the iterative  solution of  linear systems.  Hence,   SLAP.........18900
C    SLAP changes this input data structure to  the SLAP Column format   SLAP.........19000
C    for the iteration (but does not change it back).                    SLAP.........19100
C                                                                        SLAP.........19200
C    Here  is an example   of  the  SLAP  Triad storage  format  for a   SLAP.........19300
C    nonsymmetric 5x5 Matrix.  NELT=11.   Recall that the  entries may   SLAP.........19400
C    appear in any order.                                                SLAP.........19500
C                                                                        SLAP.........19600
C     5x5 Matrix       SLAP Triad format for 5x5 matrix on left.         SLAP.........19700
C                           1  2  3  4  5  6  7  8  9 10 11              SLAP.........19800
C    |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21              SLAP.........19900
C    |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2              SLAP.........20000
C    | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1              SLAP.........20100
C    | 0  0  0 44  0|                                                    SLAP.........20200
C    |51  0 53  0 55|                                                    SLAP.........20300
C                                                                        SLAP.........20400
C    ====================== S L A P Column format ====================   SLAP.........20500
C                                                                        SLAP.........20600
C    In the SLAP Column format  the non-zeros are stored counting down   SLAP.........20700
C    columns (except for the  diagonal entry,  which must appear first   SLAP.........20800
C    in each "column") and are stored in the double precision array A.   SLAP.........20900
C    In  other words,  for each  column  in the matrix  first put  the   SLAP.........21000
C    diagonal  entry  in A.  Then  put in the other  non-zero elements   SLAP.........21100
C    going  down the column  (except the  diagonal)  in order.  The IA   SLAP.........21200
C    array holds the row index  for each non-zero.  The JA array holds   SLAP.........21300
C    the  offsets  into the  IA, A  arrays for the  beginning of  each   SLAP.........21400
C    column. That is, IA(JA(ICOL)), A(JA(ICOL)) are the first elements   SLAP.........21500
C    of  the  ICOL-th  column  in  IA  and  A,  and  IA(JA(ICOL+1)-1),   SLAP.........21600
C    A(JA(ICOL+1)-1) are the last elements of the ICOL-th column. Note   SLAP.........21700
C    that we  always have  JA(N+1) = NELT+1, where  N is the number of   SLAP.........21800
C    columns in the matrix  and NELT is the number of non-zeros in the   SLAP.........21900
C    matrix.  If the matrix is symmetric one need only store the lower   SLAP.........22000
C    triangle  (including the diagonal)  and NELT would be the  corre-   SLAP.........22100
C    sponding number of non-zeros stored.                                SLAP.........22200
C                                                                        SLAP.........22300
C    Here is  an  example of the  SLAP   Column storage format  for  a   SLAP.........22400
C    nonsymmetric 5x5 Matrix (in the  A and  IA arrays '|' denotes the   SLAP.........22500
C    end of a column):                                                   SLAP.........22600
C                                                                        SLAP.........22700
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.       SLAP.........22800
C                           1  2  3    4  5    6  7    8    9 10 11      SLAP.........22900
C    |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35      SLAP.........23000
C    |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3      SLAP.........23100
C    | 0  0 33  0 35|  JA:  1  4  6    8  9   12                         SLAP.........23200
C    | 0  0  0 44  0|                                                    SLAP.........23300
C    |51  0 53  0 55|                                                    SLAP.........23400
C                                                                        SLAP.........23500
C    =================================================================   SLAP.........23600
C    ====================== Which Method To Use ======================   SLAP.........23700
C    =================================================================   SLAP.........23800
C                                                                        SLAP.........23900
C                          BACKGROUND                                    SLAP.........24000
C    In solving a large sparse linear system Ax = b using an iterative   SLAP.........24100
C    method, it   is  not necessary to actually   store  the matrix A.   SLAP.........24200
C    Rather, what is needed is a procedure  for multiplying the matrix   SLAP.........24300
C    A times a given vector y to obtain the matrix-vector product, Ay.   SLAP.........24400
C    SLAP has been written to take advantage of this fact.  The higher   SLAP.........24500
C    level routines in the package require storage only of the non-zero  SLAP.........24600
C    elements of   A (and  their  positions), and  even this   can  be   SLAP.........24700
C    avoided, if the  user  writes his own subroutine for  multiplying   SLAP.........24800
C    the matrix times a vector  and   calls the lower-level  iterative   SLAP.........24900
C    routines in the package.                                            SLAP.........25000
C                                                                        SLAP.........25100
C    If  the matrix A is ill-conditioned,  then most iterative methods   SLAP.........25200
C    will be slow to converge (if they converge  at all!).  To improve   SLAP.........25300
C    the  convergence  rate,  one  may use  a "matrix  splitting," or,   SLAP.........25400
C    "preconditioning matrix," say, M.  It is then necessary to solve,   SLAP.........25500
C    at each iteration, a linear system  with coefficient matrix M.  A   SLAP.........25600
C    good preconditioner  M should have  two  properties: (1) M should   SLAP.........25700
C    "approximate" A, in the sense that the  matrix inv(M)*A  (or some   SLAP.........25800
C    variant  thereof) is better conditioned  than the original matrix   SLAP.........25900
C    A; and  (2) linear  systems with coefficient  matrix M should  be   SLAP.........26000
C    much easier  to solve  than  the original system with coefficient   SLAP.........26100
C    matrix   A.   Preconditioning routines  in the   SLAP package are   SLAP.........26200
C    separate from the  iterative   routines,  so   that any of    the   SLAP.........26300
C    preconditioners provided in the package,   or one that the   user   SLAP.........26400
C    codes himself, can be used with any of the iterative routines.      SLAP.........26500
C                                                                        SLAP.........26600
C                        CHOICE OF PRECONDITIONER                        SLAP.........26700
C    If you  willing   to live with   either the SLAP Triad or  Column   SLAP.........26800
C    matrix data structure  you  can then  choose one  of two types of   SLAP.........26900
C    preconditioners   to   use:   diagonal  scaling    or  incomplete   SLAP.........27000
C    factorization.  To  choose   between these two   methods requires   SLAP.........27100
C    knowing  something  about the computer you're going  to run these   SLAP.........27200
C    codes on  and how well incomplete factorization  approximates the   SLAP.........27300
C    inverse of your matrix.                                             SLAP.........27400
C                                                                        SLAP.........27500
C    Let us  suppose you have   a scalar  machine.   Then,  unless the   SLAP.........27600
C    incomplete factorization is very,  very poor this  is *GENERALLY*   SLAP.........27700
C    the method to choose.  It  will reduce the  number of  iterations   SLAP.........27800
C    significantly and is not all  that expensive  to compute.  So  if   SLAP.........27900
C    you have just one  linear system to solve  and  "just want to get   SLAP.........28000
C    the job  done" then try  incomplete factorization first.   If you   SLAP.........28100
C    are thinking of integrating some SLAP  iterative method into your   SLAP.........28200
C    favorite   "production  code" then  try incomplete  factorization   SLAP.........28300
C    first,  but  also check  to see that  diagonal  scaling is indeed   SLAP.........28400
C    slower for a large sample of test problems.                         SLAP.........28500
C                                                                        SLAP.........28600
C    Let us now suppose  you have  a  vector  computer  with  hardware   SLAP.........28700
C    gather/scatter support (Cray X-MP, Y-MP, SCS-40 or Cyber 205, ETA   SLAP.........28800
C    10,  ETA Piper,  Convex C-1,  etc.).   Then  it is much harder to   SLAP.........28900
C    choose  between the  two  methods.   The  versions  of incomplete   SLAP.........29000
C    factorization in SLAP do in fact vectorize, but have short vector   SLAP.........29100
C    lengths and the factorization step is relatively  more expensive.   SLAP.........29200
C    Hence,  for  most problems (i.e.,  unless  your  problem  is  ill   SLAP.........29300
C    conditioned,  sic!)  diagonal  scaling is  faster,  with its very   SLAP.........29400
C    fast    set up  time    and  vectorized  (with   long    vectors)   SLAP.........29500
C    preconditioning step (even though  it  may take more iterations).   SLAP.........29600
C    If you have several systems (or  right hand sides) to  solve that   SLAP.........29700
C    can  utilize  the  same  preconditioner  then the   cost   of the   SLAP.........29800
C    incomplete factorization can   be  amortized over these  several    SLAP.........29900
C    solutions.  This situation gives more advantage to the incomplete   SLAP.........30000
C    factorization methods.  If  you have  a  vector  machine  without   SLAP.........30100
C    hardware  gather/scatter (Cray  1,  Cray  2  &  Cray 3) then  the   SLAP.........30200
C    advantages for incomplete factorization are even less.              SLAP.........30300
C                                                                        SLAP.........30400
C    If you're trying to shoehorn SLAP into your  favorite "production   SLAP.........30500
C    code" and can not easily generate either the SLAP Triad or Column   SLAP.........30600
C    format  then  you are  left  to   your  own  devices in terms  of   SLAP.........30700
C    preconditioning.  Also,  you may  find that the   preconditioners   SLAP.........30800
C    supplied with SLAP are not sufficient  for your problem.  In this   SLAP.........30900
C    situation we would  recommend  that you   talk  with a  numerical   SLAP.........31000
C    analyst  versed in   iterative   methods   about   writing  other   SLAP.........31100
C    preconditioning  subroutines (e.g.,  polynomial  preconditioning,   SLAP.........31200
C    shifted incomplete factorization,  SOR  or SSOR  iteration).  You   SLAP.........31300
C    can always "roll your own"  by using the "core" iterative methods   SLAP.........31400
C    and supplying your own MSOLVE and MATVEC (and possibly MTSOLV and   SLAP.........31500
C    MTTVEC) routines.                                                   SLAP.........31600
C                                                                        SLAP.........31700
C                          SYMMETRIC SYSTEMS                             SLAP.........31800
C    If your matrix is symmetric then you would want to use one of the   SLAP.........31900
C    symmetric system  solvers.    If  your  system  is  also positive   SLAP.........32000
C    definite,   (Ax,x) (Ax dot  product  with x) is  positive for all   SLAP.........32100
C    non-zero  vectors x,  then use   Conjugate Gradient (DCG,  DSDCG,   SLAP.........32200
C    DSICSG).  If you're  not sure it's SPD   (symmetric and  Positive   SLAP.........32300
C    Definite)  then try DCG anyway and  if it works, fine.  If you're   SLAP.........32400
C    sure your matrix is not  positive definite  then you  may want to   SLAP.........32500
C    try the iterative refinement   methods  (DIR)  or the  GMRES code   SLAP.........32600
C    (DGMRES) if DIR converges too slowly.                               SLAP.........32700
C                                                                        SLAP.........32800
C                         NONSYMMETRIC SYSTEMS                           SLAP.........32900
C    This   is currently  an  area  of  active research  in  numerical   SLAP.........33000
C    analysis  and   there   are   new  strategies  being   developed.   SLAP.........33100
C    Consequently take the following advice with a grain of salt.   If   SLAP.........33200
C    you matrix is positive definite, (Ax,x)  (Ax  dot product  with x   SLAP.........33300
C    is positive for all non-zero  vectors x), then you can use any of   SLAP.........33400
C    the    methods   for   nonsymmetric   systems (Orthomin,   GMRES,   SLAP.........33500
C    BiConjugate Gradient, BiConjugate Gradient  Squared and Conjugate   SLAP.........33600
C    Gradient applied to the normal equations).  If your system is not   SLAP.........33700
C    too ill conditioned then try  BiConjugate Gradient Squared (BCGS)   SLAP.........33800
C    or GMRES (DGMRES).  Both  of  these methods converge very quickly   SLAP.........33900
C    and do  not require A'  or M' ('  denotes transpose) information.   SLAP.........34000
C    DGMRES  does require  some  additional storage,  though.  If  the   SLAP.........34100
C    system is very  ill conditioned  or   nearly positive  indefinite   SLAP.........34200
C    ((Ax,x) is positive,  but may be  very small),  then GMRES should   SLAP.........34300
C    be the first choice,  but try the  other  methods  if you have to   SLAP.........34400
C    fine tune  the solution process for a  "production code".  If you   SLAP.........34500
C    have a great preconditioner for the normal  equations (i.e., M is   SLAP.........34600
C    an approximation to the inverse of AA' rather than  just  A) then   SLAP.........34700
C    this is not a bad route to travel.  Old wisdom would say that the   SLAP.........34800
C    normal equations are a disaster  (since it squares the  condition   SLAP.........34900
C    number of the system and DCG convergence is linked to this number   SLAP.........35000
C    of    infamy), but   some     preconditioners    (like incomplete   SLAP.........35100
C    factorization) can reduce the condition number back below that of   SLAP.........35200
C    the original system.                                                SLAP.........35300
C                                                                        SLAP.........35400
C    =================================================================   SLAP.........35500
C    ======================= Naming Conventions ======================   SLAP.........35600
C    =================================================================   SLAP.........35700
C    SLAP  iterative  methods,    matrix vector    and  preconditioner   SLAP.........35800
C    calculation  routines   follow a naming   convention  which, when   SLAP.........35900
C    understood, allows one to determine the iterative method and data   SLAP.........36000
C    structure(s) used.  The  subroutine  naming convention  takes the   SLAP.........36100
C    following form:                                                     SLAP.........36200
C                          P[S][M]DESC                                   SLAP.........36300
C    where                                                               SLAP.........36400
C        P  stands for the precision (or data type) of the routine and   SLAP.........36500
C           is required in all names,                                    SLAP.........36600
C        S  denotes whether or not the routine requires the SLAP Triad   SLAP.........36700
C           or Column format (it does if the second letter of the name   SLAP.........36800
C           is S and does not otherwise),                                SLAP.........36900
C        M  stands for the type of preconditioner used (only appears     SLAP.........37000
C           in drivers for "core" routines), and                         SLAP.........37100
C     DESC  is some number of letters describing the method or purpose   SLAP.........37200
C           of the routine.  The following is a list of the "DESC"       SLAP.........37300
C           fields for iterative methods and their meaning:              SLAP.........37400
C             BCG,BC:       BiConjugate Gradient                         SLAP.........37500
C             CG:           Conjugate Gradient                           SLAP.........37600
C             CGN,CN:       Conjugate Gradient on the Normal equations   SLAP.........37700
C             CGS,CS:       biConjugate Gradient Squared                 SLAP.........37800
C             GMRES,GMR,GM: Generalized Minimum RESidual                 SLAP.........37900
C             IR,R:         Iterative Refinement                         SLAP.........38000
C             JAC:          JACobi's method                              SLAP.........38100
C             GS:           Gauss-Seidel                                 SLAP.........38200
C             OMN,OM:       OrthoMiN                                     SLAP.........38300
C                                                                        SLAP.........38400
C    In the double precision version of SLAP, all routine names start    SLAP.........38500
C    with a D. The brackets around the S and M designate that these      SLAP.........38600
C    fields are optional.                                                SLAP.........38700
C                                                                        SLAP.........38800
C    Here are some examples of the routines:                             SLAP.........38900
C    1) DBCG: Double precision BiConjugate Gradient "core" routine.      SLAP.........39000
C       One can deduce that this is a "core" routine, because the S and  SLAP.........39100
C       M fields are missing and BiConjugate Gradient is an iterative    SLAP.........39200
C       method.                                                          SLAP.........39300
C    2) DSDBCG: Double precision, SLAP data structure BCG with Diagonal  SLAP.........39400
C       scaling.                                                         SLAP.........39500
C    3) DSLUBC: Double precision, SLAP data structure BCG with incom-    SLAP.........39600
C       plete LU factorization as the preconditioning.                   SLAP.........39700
C    4) DCG: Double precision Conjugate Gradient "core" routine.         SLAP.........39800
C    5) DSDCG: Double precision, SLAP data structure Conjugate Gradient  SLAP.........39900
C       with Diagonal scaling.                                           SLAP.........40000
C    6) DSICCG: Double precision, SLAP data structure Conjugate Gra-     SLAP.........40100
C       dient with Incomplete Cholesky factorization preconditioning.    SLAP.........40200
C                                                                        SLAP.........40300
C                                                                        SLAP.........40400
C    =================================================================   SLAP.........40500
C    ===================== USER CALLABLE ROUTINES ====================   SLAP.........40600
C    =================================================================   SLAP.........40700
C    The following is a list of  the "user callable" SLAP routines and   SLAP.........40800
C    their one line descriptions.  The headers denote  the  file names   SLAP.........40900
C    where the routines can be found, as distributed for UNIX systems.   SLAP.........41000
C                                                                        SLAP.........41100
C    Note:  Each core routine, DXXX, has a corresponding stop routine,   SLAP.........41200
C         ISDXXX.  If the stop routine does not have the specific stop   SLAP.........41300
C         test the user requires (e.g., weighted infinity norm),  then   SLAP.........41400
C         the user should modify the source for ISDXXX accordingly.      SLAP.........41500
C                                                                        SLAP.........41600
C    ============================= dir.f =============================   SLAP.........41700
C    DIR: Preconditioned Iterative Refinement Sparse Ax = b Solver.      SLAP.........41800
C    DSJAC: Jacobi's Method Iterative Sparse Ax = b Solver.              SLAP.........41900
C    DSGS: Gauss-Seidel Method Iterative Sparse Ax = b Solver.           SLAP.........42000
C    DSILUR: Incomplete LU Iterative Refinement Sparse Ax = b Solver.    SLAP.........42100
C                                                                        SLAP.........42200
C    ============================= dcg.f =============================   SLAP.........42300
C    DCG: Preconditioned Conjugate Gradient Sparse Ax=b Solver.          SLAP.........42400
C    DSDCG: Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver.     SLAP.........42500
C    DSICCG: Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver.  SLAP.........42600
C                                                                        SLAP.........42700
C    ============================= dcgn.f ============================   SLAP.........42800
C    DCGN: Preconditioned CG Sparse Ax=b Solver for Normal Equations.    SLAP.........42900
C    DSDCGN: Diagonally Scaled CG Sparse Ax=b Solver for Normal Eqn's.   SLAP.........43000
C    DSLUCN: Incomplete LU CG Sparse Ax=b Solver for Normal Equations.   SLAP.........43100
C                                                                        SLAP.........43200
C    ============================= dbcg.f ============================   SLAP.........43300
C    DBCG: Preconditioned BiConjugate Gradient Sparse Ax = b Solver.     SLAP.........43400
C    DSDBCG: Diagonally Scaled BiConjugate Gradient Sparse Ax=b Solver.  SLAP.........43500
C    DSLUBC: Incomplete LU BiConjugate Gradient Sparse Ax=b Solver.      SLAP.........43600
C                                                                        SLAP.........43700
C    ============================= dcgs.f ============================   SLAP.........43800
C    DCGS: Preconditioned BiConjugate Gradient Squared Ax=b Solver.      SLAP.........43900
C    DSDCGS: Diagonally Scaled CGS Sparse Ax=b Solver.                   SLAP.........44000
C    DSLUCS: Incomplete LU BiConjugate Gradient Squared Ax=b Solver.     SLAP.........44100
C                                                                        SLAP.........44200
C    ============================= domn.f ============================   SLAP.........44300
C    DOMN: Preconditioned Orthomin Sparse Iterative Ax=b Solver.         SLAP.........44400
C    DSDOMN: Diagonally Scaled Orthomin Sparse Iterative Ax=b Solver.    SLAP.........44500
C    DSLUOM: Incomplete LU Orthomin Sparse Iterative Ax=b Solver.        SLAP.........44600
C                                                                        SLAP.........44700
C    ============================ dgmres.f ===========================   SLAP.........44800
C    DGMRES: Preconditioned GMRES Iterative Sparse Ax=b Solver.          SLAP.........44900
C    DSDGMR: Diagonally Scaled GMRES Iterative Sparse Ax=b Solver.       SLAP.........45000
C    DSLUGM: Incomplete LU GMRES Iterative Sparse Ax=b Solver.           SLAP.........45100
C                                                                        SLAP.........45200
C    ============================ dmset.f ============================   SLAP.........45300
C       The following routines are used to set up preconditioners.       SLAP.........45400
C                                                                        SLAP.........45500
C    DSDS: Diagonal Scaling Preconditioner SLAP Set Up.                  SLAP.........45600
C    DSDSCL: Diagonally Scales/Unscales a SLAP Column Matrix.            SLAP.........45700
C    DSD2S: Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up.     SLAP.........45800
C    DS2LT: Lower Triangle Preconditioner SLAP Set Up.                   SLAP.........45900
C    DSICS: Incomplete Cholesky Decomp. Preconditioner SLAP Set Up.      SLAP.........46000
C    DSILUS: Incomplete LU Decomposition Preconditioner SLAP Set Up.     SLAP.........46100
C                                                                        SLAP.........46200
C    ============================ dmvops.f ===========================   SLAP.........46300
C       Most of the incomplete  factorization  (LL' and LDU) solvers     SLAP.........46400
C       in this  file require an  intermediate routine  to translate     SLAP.........46500
C       from the SLAP MSOLVE(N, R, Z, NELT, IA,  JA, A, ISYM, RWORK,     SLAP.........46600
C       IWORK) calling  convention to the calling  sequence required     SLAP.........46700
C       by  the solve routine.   This generally  is  accomplished by     SLAP.........46800
C       fishing out pointers to the preconditioner (stored in RWORK)     SLAP.........46900
C       from the  IWORK  array and then making a call to the routine     SLAP.........47000
C       that actually does the backsolve.                                SLAP.........47100
C                                                                        SLAP.........47200
C    DSMV: SLAP Column Format Sparse Matrix Vector Product.              SLAP.........47300
C    DSMTV: SLAP Column Format Sparse Matrix (transpose) Vector Prod.    SLAP.........47400
C    DSDI: Diagonal Matrix Vector Multiply.                              SLAP.........47500
C    DSLI: SLAP MSOLVE for Lower Triangle Matrix (set up for DSLI2).     SLAP.........47600
C    DSLI2: Lower Triangle Matrix Backsolve.                             SLAP.........47700
C    DSLLTI: SLAP MSOLVE for LDL' (IC) Fact. (set up for DLLTI2).        SLAP.........47800
C    DLLTI2: Backsolve routine for LDL' Factorization.                   SLAP.........47900
C    DSLUI: SLAP MSOLVE for LDU Factorization (set up for DSLUI2).       SLAP.........48000
C    DSLUI2: SLAP Backsolve for LDU Factorization.                       SLAP.........48100
C    DSLUTI: SLAP MTSOLV for LDU Factorization (set up for DSLUI4).      SLAP.........48200
C    DSLUI4: SLAP Backsolve for LDU Factorization.                       SLAP.........48300
C    DSMMTI: SLAP MSOLVE for LDU Fact of Normal Eq (set up for DSMMI2).  SLAP.........48400
C    DSMMI2: SLAP Backsolve for LDU Factorization of Normal Equations.   SLAP.........48500
C                                                                        SLAP.........48600
C    =========================== dlaputil.f ==========================   SLAP.........48700
C       The following utility routines are useful additions to SLAP.     SLAP.........48800
C                                                                        SLAP.........48900
C    DBHIN: Read Sparse Linear System in the Boeing/Harwell Format.      SLAP.........49000
C    DCHKW: SLAP WORK/IWORK Array Bounds Checker.                        SLAP.........49100
C    DCPPLT: Printer Plot of SLAP Column Format Matrix.                  SLAP.........49200
C    DS2Y: SLAP Triad to SLAP Column Format Converter.                   SLAP.........49300
C    QS2I1D: Quick Sort Integer array, moving integer and DP arrays.     SLAP.........49400
C            (Used by DS2Y.)                                             SLAP.........49500
C    DTIN: Read in SLAP Triad Format Linear System.                      SLAP.........49600
C    DTOUT: Write out SLAP Triad Format Linear System.                   SLAP.........49700
C                                                                        SLAP.........49800
C                                                                        SLAP.........49900
C***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in             SLAP.........50000
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,   SLAP.........50100
C                  Algorithms and Applications, Wiley, 1989, pp.135-155. SLAP.........50200
C***ROUTINES CALLED  (NONE)                                              SLAP.........50300
C***REVISION HISTORY  (YYMMDD)                                           SLAP.........50400
C   890404  DATE WRITTEN                                                 SLAP.........50500
C   890404  Previous REVISION DATE                                       SLAP.........50600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP.........50700
C   890921  Removed TeX from comments.  (FNF)                            SLAP.........50800
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP.........50900
C           standard.  (FNF)                                             SLAP.........51000
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP.........51100
C           -----( This produced Version 2.0.1. )-----                   SLAP.........51200
C   891003  Rearranged list of user callable routines to agree with      SLAP.........51300
C           order in source deck.  (FNF)                                 SLAP.........51400
C   891004  Updated reference.                                           SLAP.........51500
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP.........51600
C           -----( This produced Version 2.0.2. )-----                   SLAP.........51700
C   910506  Minor improvements to prologue.  (FNF)                       SLAP.........51800
C   920511  Added complete declaration section.  (WRB)                   SLAP.........51900
C   920929  Corrected format of reference.  (FNF)                        SLAP.........52000
C   921019  Improved one-line descriptions, reordering some.  (FNF)      SLAP.........52100
C***END PROLOGUE  DLPDOC                                                 SLAP.........52200
C***FIRST EXECUTABLE STATEMENT  DLPDOC                                   SLAP.........52300
C                                                                        SLAP.........52400
C     This is a *DUMMY* subroutine and should never be called.           SLAP.........52500
C                                                                        SLAP.........52600
      RETURN                                                             SLAP.........52700
C------------- LAST LINE OF DLPDOC FOLLOWS ----------------------------- SLAP.........52800
      END                                                                SLAP.........52900
*DECK D1MACH                                                             SLAP.........53000
      DOUBLE PRECISION FUNCTION D1MACH (I)                               SLAP.........53100
C***BEGIN PROLOGUE  D1MACH                                               SLAP.........53200
C***PURPOSE  Return floating point machine dependent constants.          SLAP.........53300
C***LIBRARY   SLATEC                                                     SLAP.........53400
C***CATEGORY  R1                                                         SLAP.........53500
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)                      SLAP.........53600
C***KEYWORDS  MACHINE CONSTANTS                                          SLAP.........53700
C***AUTHOR  Fox, P. A., (Bell Labs)                                      SLAP.........53800
C           Hall, A. D., (Bell Labs)                                     SLAP.........53900
C           Schryer, N. L., (Bell Labs)                                  SLAP.........54000
C***DESCRIPTION                                                          SLAP.........54100
C                                                                        SLAP.........54200
C   D1MACH can be used to obtain machine-dependent parameters for the    SLAP.........54300
C   local machine environment.  It is a function subprogram with one     SLAP.........54400
C   (input) argument, and can be referenced as follows:                  SLAP.........54500
C                                                                        SLAP.........54600
C        D = D1MACH(I)                                                   SLAP.........54700
C                                                                        SLAP.........54800
C   where I=1,...,5.  The (output) value of D above is determined by     SLAP.........54900
C   the (input) value of I.  The results for various values of I are     SLAP.........55000
C   discussed below.                                                     SLAP.........55100
C                                                                        SLAP.........55200
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.           SLAP.........55300
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.           SLAP.........55400
C   D1MACH( 3) = B**(-T), the smallest relative spacing.                 SLAP.........55500
C   D1MACH( 4) = B**(1-T), the largest relative spacing.                 SLAP.........55600
C   D1MACH( 5) = LOG10(B)                                                SLAP.........55700
C                                                                        SLAP.........55800
C   Assume double precision numbers are represented in the T-digit,      SLAP.........55900
C   base-B form                                                          SLAP.........56000
C                                                                        SLAP.........56100
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )              SLAP.........56200
C                                                                        SLAP.........56300
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and             SLAP.........56400
C   EMIN .LE. E .LE. EMAX.                                               SLAP.........56500
C                                                                        SLAP.........56600
C   The values of B, T, EMIN and EMAX are provided in I1MACH as          SLAP.........56700
C   follows:                                                             SLAP.........56800
C   I1MACH(10) = B, the base.                                            SLAP.........56900
C   I1MACH(14) = T, the number of base-B digits.                         SLAP.........57000
C   I1MACH(15) = EMIN, the smallest exponent E.                          SLAP.........57100
C   I1MACH(16) = EMAX, the largest exponent E.                           SLAP.........57200
C                                                                        SLAP.........57300
C   To alter this function for a particular environment, the desired     SLAP.........57400
C   set of DATA statements should be activated by removing the C from    SLAP.........57500
C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be       SLAP.........57600
C   checked for consistency with the local operating system.             SLAP.........57700
C                                                                        SLAP.........57800
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for   SLAP.........57900
C                 a portable library, ACM Transactions on Mathematical   SLAP.........58000
C                 Software 4, 2 (June 1978), pp. 177-188.                SLAP.........58100
C***ROUTINES CALLED  XERMSG                                              SLAP.........58200
C***REVISION HISTORY  (YYMMDD)                                           SLAP.........58300
C   750101  DATE WRITTEN                                                 SLAP.........58400
C   890213  REVISION DATE from Version 3.2                               SLAP.........58500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP.........58600
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)           SLAP.........58700
C   900618  Added DEC RISC constants.  (WRB)                             SLAP.........58800
C   900723  Added IBM RS 6000 constants.  (WRB)                          SLAP.........58900
C   900911  Added SUN 386i constants.  (WRB)                             SLAP.........59000
C   910710  Added HP 730 constants.  (SMR)                               SLAP.........59100
C   911114  Added Convex IEEE constants.  (WRB)                          SLAP.........59200
C   920121  Added SUN -r8 compiler option constants.  (WRB)              SLAP.........59300
C   920229  Added Touchstone Delta i860 constants.  (WRB)                SLAP.........59400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP.........59500
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.         SLAP.........59600
C           (BKS, WRB)                                                   SLAP.........59700
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)            SLAP.........59800
C***END PROLOGUE  D1MACH                                                 SLAP.........59900
C                                                                        SLAP.........60000
      INTEGER SMALL(4)                                                   SLAP.........60100
      INTEGER LARGE(4)                                                   SLAP.........60200
      INTEGER RIGHT(4)                                                   SLAP.........60300
      INTEGER DIVER(4)                                                   SLAP.........60400
      INTEGER LOG10(4)                                                   SLAP.........60500
C                                                                        SLAP.........60600
      DOUBLE PRECISION DMACH(5)                                          SLAP.........60700
      SAVE DMACH                                                         SLAP.........60800
C                                                                        SLAP.........60900
      EQUIVALENCE (DMACH(1),SMALL(1))                                    SLAP.........61000
      EQUIVALENCE (DMACH(2),LARGE(1))                                    SLAP.........61100
      EQUIVALENCE (DMACH(3),RIGHT(1))                                    SLAP.........61200
      EQUIVALENCE (DMACH(4),DIVER(1))                                    SLAP.........61300
      EQUIVALENCE (DMACH(5),LOG10(1))                                    SLAP.........61400
C                                                                        SLAP.........61500
C     MACHINE CONSTANTS FOR THE AMIGA                                    SLAP.........61600
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION      SLAP.........61700
C                                                                        SLAP.........61800
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /               SLAP.........61900
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /               SLAP.........62000
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /               SLAP.........62100
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /               SLAP.........62200
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /               SLAP.........62300
C                                                                        SLAP.........62400
C     MACHINE CONSTANTS FOR THE AMIGA                                    SLAP.........62500
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT              SLAP.........62600
C                                                                        SLAP.........62700
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /               SLAP.........62800
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /               SLAP.........62900
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /               SLAP.........63000
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /               SLAP.........63100
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /               SLAP.........63200
C                                                                        SLAP.........63300
C     MACHINE CONSTANTS FOR THE APOLLO                                   SLAP.........63400
C                                                                        SLAP.........63500
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /               SLAP.........63600
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /               SLAP.........63700
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /               SLAP.........63800
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /               SLAP.........63900
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /               SLAP.........64000
C                                                                        SLAP.........64100
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM                    SLAP.........64200
C                                                                        SLAP.........64300
C     DATA SMALL(1) / ZC00800000 /                                       SLAP.........64400
C     DATA SMALL(2) / Z000000000 /                                       SLAP.........64500
C     DATA LARGE(1) / ZDFFFFFFFF /                                       SLAP.........64600
C     DATA LARGE(2) / ZFFFFFFFFF /                                       SLAP.........64700
C     DATA RIGHT(1) / ZCC5800000 /                                       SLAP.........64800
C     DATA RIGHT(2) / Z000000000 /                                       SLAP.........64900
C     DATA DIVER(1) / ZCC6800000 /                                       SLAP.........65000
C     DATA DIVER(2) / Z000000000 /                                       SLAP.........65100
C     DATA LOG10(1) / ZD00E730E7 /                                       SLAP.........65200
C     DATA LOG10(2) / ZC77800DC0 /                                       SLAP.........65300
C                                                                        SLAP.........65400
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM                    SLAP.........65500
C                                                                        SLAP.........65600
C     DATA SMALL(1) / O1771000000000000 /                                SLAP.........65700
C     DATA SMALL(2) / O0000000000000000 /                                SLAP.........65800
C     DATA LARGE(1) / O0777777777777777 /                                SLAP.........65900
C     DATA LARGE(2) / O0007777777777777 /                                SLAP.........66000
C     DATA RIGHT(1) / O1461000000000000 /                                SLAP.........66100
C     DATA RIGHT(2) / O0000000000000000 /                                SLAP.........66200
C     DATA DIVER(1) / O1451000000000000 /                                SLAP.........66300
C     DATA DIVER(2) / O0000000000000000 /                                SLAP.........66400
C     DATA LOG10(1) / O1157163034761674 /                                SLAP.........66500
C     DATA LOG10(2) / O0006677466732724 /                                SLAP.........66600
C                                                                        SLAP.........66700
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS              SLAP.........66800
C                                                                        SLAP.........66900
C     DATA SMALL(1) / O1771000000000000 /                                SLAP.........67000
C     DATA SMALL(2) / O7770000000000000 /                                SLAP.........67100
C     DATA LARGE(1) / O0777777777777777 /                                SLAP.........67200
C     DATA LARGE(2) / O7777777777777777 /                                SLAP.........67300
C     DATA RIGHT(1) / O1461000000000000 /                                SLAP.........67400
C     DATA RIGHT(2) / O0000000000000000 /                                SLAP.........67500
C     DATA DIVER(1) / O1451000000000000 /                                SLAP.........67600
C     DATA DIVER(2) / O0000000000000000 /                                SLAP.........67700
C     DATA LOG10(1) / O1157163034761674 /                                SLAP.........67800
C     DATA LOG10(2) / O0006677466732724 /                                SLAP.........67900
C                                                                        SLAP.........68000
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE          SLAP.........68100
C                                                                        SLAP.........68200
C     DATA SMALL(1) / Z"3001800000000000" /                              SLAP.........68300
C     DATA SMALL(2) / Z"3001000000000000" /                              SLAP.........68400
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /                              SLAP.........68500
C     DATA LARGE(2) / Z"4FFE000000000000" /                              SLAP.........68600
C     DATA RIGHT(1) / Z"3FD2800000000000" /                              SLAP.........68700
C     DATA RIGHT(2) / Z"3FD2000000000000" /                              SLAP.........68800
C     DATA DIVER(1) / Z"3FD3800000000000" /                              SLAP.........68900
C     DATA DIVER(2) / Z"3FD3000000000000" /                              SLAP.........69000
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /                              SLAP.........69100
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /                              SLAP.........69200
C                                                                        SLAP.........69300
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES                     SLAP.........69400
C                                                                        SLAP.........69500
C     DATA SMALL(1) / 00564000000000000000B /                            SLAP.........69600
C     DATA SMALL(2) / 00000000000000000000B /                            SLAP.........69700
C     DATA LARGE(1) / 37757777777777777777B /                            SLAP.........69800
C     DATA LARGE(2) / 37157777777777777777B /                            SLAP.........69900
C     DATA RIGHT(1) / 15624000000000000000B /                            SLAP.........70000
C     DATA RIGHT(2) / 00000000000000000000B /                            SLAP.........70100
C     DATA DIVER(1) / 15634000000000000000B /                            SLAP.........70200
C     DATA DIVER(2) / 00000000000000000000B /                            SLAP.........70300
C     DATA LOG10(1) / 17164642023241175717B /                            SLAP.........70400
C     DATA LOG10(2) / 16367571421742254654B /                            SLAP.........70500
C                                                                        SLAP.........70600
C     MACHINE CONSTANTS FOR THE CELERITY C1260                           SLAP.........70700
C                                                                        SLAP.........70800
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /               SLAP.........70900
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /               SLAP.........71000
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /               SLAP.........71100
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /               SLAP.........71200
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /               SLAP.........71300
C                                                                        SLAP.........71400
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP.........71500
C     USING THE -fn OR -pd8 COMPILER OPTION                              SLAP.........71600
C                                                                        SLAP.........71700
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........71800
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /                              SLAP.........71900
C     DATA DMACH(3) / Z'3CC0000000000000' /                              SLAP.........72000
C     DATA DMACH(4) / Z'3CD0000000000000' /                              SLAP.........72100
C     DATA DMACH(5) / Z'3FF34413509F79FF' /                              SLAP.........72200
C                                                                        SLAP.........72300
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP.........72400
C     USING THE -fi COMPILER OPTION                                      SLAP.........72500
C                                                                        SLAP.........72600
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........72700
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /                              SLAP.........72800
C     DATA DMACH(3) / Z'3CA0000000000000' /                              SLAP.........72900
C     DATA DMACH(4) / Z'3CB0000000000000' /                              SLAP.........73000
C     DATA DMACH(5) / Z'3FD34413509F79FF' /                              SLAP.........73100
C                                                                        SLAP.........73200
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP.........73300
C     USING THE -p8 COMPILER OPTION                                      SLAP.........73400
C                                                                        SLAP.........73500
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /              SLAP.........73600
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /              SLAP.........73700
C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /              SLAP.........73800
C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /              SLAP.........73900
C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /              SLAP.........74000
C                                                                        SLAP.........74100
C     MACHINE CONSTANTS FOR THE CRAY                                     SLAP.........74200
C                                                                        SLAP.........74300
C     DATA SMALL(1) / 201354000000000000000B /                           SLAP.........74400
C     DATA SMALL(2) / 000000000000000000000B /                           SLAP.........74500
C     DATA LARGE(1) / 577767777777777777777B /                           SLAP.........74600
C     DATA LARGE(2) / 000007777777777777774B /                           SLAP.........74700
C     DATA RIGHT(1) / 376434000000000000000B /                           SLAP.........74800
C     DATA RIGHT(2) / 000000000000000000000B /                           SLAP.........74900
C     DATA DIVER(1) / 376444000000000000000B /                           SLAP.........75000
C     DATA DIVER(2) / 000000000000000000000B /                           SLAP.........75100
C     DATA LOG10(1) / 377774642023241175717B /                           SLAP.........75200
C     DATA LOG10(2) / 000007571421742254654B /                           SLAP.........75300
C                                                                        SLAP.........75400
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200               SLAP.........75500
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -       SLAP.........75600
C     STATIC DMACH(5)                                                    SLAP.........75700
C                                                                        SLAP.........75800
C     DATA SMALL /    20K, 3*0 /                                         SLAP.........75900
C     DATA LARGE / 77777K, 3*177777K /                                   SLAP.........76000
C     DATA RIGHT / 31420K, 3*0 /                                         SLAP.........76100
C     DATA DIVER / 32020K, 3*0 /                                         SLAP.........76200
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /                      SLAP.........76300
C                                                                        SLAP.........76400
C     MACHINE CONSTANTS FOR THE DEC ALPHA                                SLAP.........76500
C     USING G_FLOAT                                                      SLAP.........76600
C                                                                        SLAP.........76700
C     DATA DMACH(1) / '0000000000000010'X /                              SLAP.........76800
C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /                              SLAP.........76900
C     DATA DMACH(3) / '0000000000003CC0'X /                              SLAP.........77000
C     DATA DMACH(4) / '0000000000003CD0'X /                              SLAP.........77100
C     DATA DMACH(5) / '79FF509F44133FF3'X /                              SLAP.........77200
C                                                                        SLAP.........77300
C     MACHINE CONSTANTS FOR THE DEC ALPHA                                SLAP.........77400
C     USING IEEE_FORMAT                                                  SLAP.........77500
C                                                                        SLAP.........77600
C     DATA DMACH(1) / '0010000000000000'X /                              SLAP.........77700
C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /                              SLAP.........77800
C     DATA DMACH(3) / '3CA0000000000000'X /                              SLAP.........77900
C     DATA DMACH(4) / '3CB0000000000000'X /                              SLAP.........78000
C     DATA DMACH(5) / '3FD34413509F79FF'X /                              SLAP.........78100
C                                                                        SLAP.........78200
C     MACHINE CONSTANTS FOR THE DEC RISC                                 SLAP.........78300
C                                                                        SLAP.........78400
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/                SLAP.........78500
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/                SLAP.........78600
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/                SLAP.........78700
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/                SLAP.........78800
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/                SLAP.........78900
C                                                                        SLAP.........79000
C     MACHINE CONSTANTS FOR THE DEC VAX                                  SLAP.........79100
C     USING D_FLOATING                                                   SLAP.........79200
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)                             SLAP.........79300
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS          SLAP.........79400
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS                   SLAP.........79500
C                                                                        SLAP.........79600
C     DATA SMALL(1), SMALL(2) /        128,           0 /                SLAP.........79700
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /                SLAP.........79800
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /                SLAP.........79900
C     DATA DIVER(1), DIVER(2) /       9472,           0 /                SLAP.........80000
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /                SLAP.........80100
C                                                                        SLAP.........80200
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /                   SLAP.........80300
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /                   SLAP.........80400
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /                   SLAP.........80500
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /                   SLAP.........80600
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /                   SLAP.........80700
C                                                                        SLAP.........80800
C     MACHINE CONSTANTS FOR THE DEC VAX                                  SLAP.........80900
C     USING G_FLOATING                                                   SLAP.........81000
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)                             SLAP.........81100
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS          SLAP.........81200
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS                   SLAP.........81300
C                                                                        SLAP.........81400
C     DATA SMALL(1), SMALL(2) /         16,           0 /                SLAP.........81500
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /                SLAP.........81600
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /                SLAP.........81700
C     DATA DIVER(1), DIVER(2) /      15568,           0 /                SLAP.........81800
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /                SLAP.........81900
C                                                                        SLAP.........82000
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /                   SLAP.........82100
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /                   SLAP.........82200
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /                   SLAP.........82300
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /                   SLAP.........82400
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /                   SLAP.........82500
C                                                                        SLAP.........82600
C     MACHINE CONSTANTS FOR THE ELXSI 6400                               SLAP.........82700
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)                  SLAP.........82800
C                                                                        SLAP.........82900
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /                SLAP.........83000
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /                SLAP.........83100
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /                SLAP.........83200
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /                SLAP.........83300
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /                SLAP.........83400
C                                                                        SLAP.........83500
C     MACHINE CONSTANTS FOR THE HARRIS 220                               SLAP.........83600
C                                                                        SLAP.........83700
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /                   SLAP.........83800
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /                   SLAP.........83900
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /                   SLAP.........84000
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /                   SLAP.........84100
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /                   SLAP.........84200
C                                                                        SLAP.........84300
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES                SLAP.........84400
C                                                                        SLAP.........84500
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /           SLAP.........84600
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /           SLAP.........84700
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /           SLAP.........84800
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /           SLAP.........84900
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /           SLAP.........85000
C                                                                        SLAP.........85100
C     MACHINE CONSTANTS FOR THE HP 730                                   SLAP.........85200
C                                                                        SLAP.........85300
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........85400
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /                              SLAP.........85500
C     DATA DMACH(3) / Z'3CA0000000000000' /                              SLAP.........85600
C     DATA DMACH(4) / Z'3CB0000000000000' /                              SLAP.........85700
C     DATA DMACH(5) / Z'3FD34413509F79FF' /                              SLAP.........85800
C                                                                        SLAP.........85900
C     MACHINE CONSTANTS FOR THE HP 2100                                  SLAP.........86000
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4                       SLAP.........86100
C                                                                        SLAP.........86200
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /     SLAP.........86300
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /     SLAP.........86400
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /     SLAP.........86500
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /     SLAP.........86600
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /     SLAP.........86700
C                                                                        SLAP.........86800
C     MACHINE CONSTANTS FOR THE HP 2100                                  SLAP.........86900
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4                        SLAP.........87000
C                                                                        SLAP.........87100
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /                       SLAP.........87200
C     DATA SMALL(3), SMALL(4) /       0,       1 /                       SLAP.........87300
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /                       SLAP.........87400
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /                       SLAP.........87500
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /                       SLAP.........87600
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /                       SLAP.........87700
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /                       SLAP.........87800
C     DATA DIVER(3), DIVER(4) /       0,    227B /                       SLAP.........87900
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /                       SLAP.........88000
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /                       SLAP.........88100
C                                                                        SLAP.........88200
C     MACHINE CONSTANTS FOR THE HP 9000                                  SLAP.........88300
C                                                                        SLAP.........88400
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /             SLAP.........88500
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /             SLAP.........88600
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /             SLAP.........88700
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /             SLAP.........88800
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /             SLAP.........88900
C                                                                        SLAP.........89000
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                      SLAP.........89100
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND                  SLAP.........89200
C     THE PERKIN ELMER (INTERDATA) 7/32.                                 SLAP.........89300
C                                                                        SLAP.........89400
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /                   SLAP.........89500
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /                   SLAP.........89600
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /                   SLAP.........89700
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /                   SLAP.........89800
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /                   SLAP.........89900
C                                                                        SLAP.........90000
C     MACHINE CONSTANTS FOR THE IBM PC                                   SLAP.........90100
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION            SLAP.........90200
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.                    SLAP.........90300
C                                                                        SLAP.........90400
C     DATA SMALL(1) / 2.23D-308  /                                       SLAP.........90500
C     DATA LARGE(1) / 1.79D+308  /                                       SLAP.........90600
C     DATA RIGHT(1) / 1.11D-16   /                                       SLAP.........90700
C     DATA DIVER(1) / 2.22D-16   /                                       SLAP.........90800
C     DATA LOG10(1) / 0.301029995663981195D0 /                           SLAP.........90900
C                                                                        SLAP.........91000
C     MACHINE CONSTANTS FOR THE IBM RS 6000                              SLAP.........91100
C                                                                        SLAP.........91200
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........91300
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /                              SLAP.........91400
C     DATA DMACH(3) / Z'3CA0000000000000' /                              SLAP.........91500
C     DATA DMACH(4) / Z'3CB0000000000000' /                              SLAP.........91600
C     DATA DMACH(5) / Z'3FD34413509F79FF' /                              SLAP.........91700
C                                                                        SLAP.........91800
C     MACHINE CONSTANTS FOR THE INTEL i860                               SLAP.........91900
C                                                                        SLAP.........92000
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........92100
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /                              SLAP.........92200
C     DATA DMACH(3) / Z'3CA0000000000000' /                              SLAP.........92300
C     DATA DMACH(4) / Z'3CB0000000000000' /                              SLAP.........92400
C     DATA DMACH(5) / Z'3FD34413509F79FF' /                              SLAP.........92500
C                                                                        SLAP.........92600
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)                    SLAP.........92700
C                                                                        SLAP.........92800
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /           SLAP.........92900
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /           SLAP.........93000
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /           SLAP.........93100
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /           SLAP.........93200
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /           SLAP.........93300
C                                                                        SLAP.........93400
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)                    SLAP.........93500
C                                                                        SLAP.........93600
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /           SLAP.........93700
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /           SLAP.........93800
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /           SLAP.........93900
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /           SLAP.........94000
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /           SLAP.........94100
C                                                                        SLAP.........94200
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                    SLAP.........94300
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                  SLAP.........94400
C                                                                        SLAP.........94500
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /                SLAP.........94600
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /                SLAP.........94700
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /                SLAP.........94800
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /                SLAP.........94900
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /                SLAP.........95000
C                                                                        SLAP.........95100
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /             SLAP.........95200
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /             SLAP.........95300
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /             SLAP.........95400
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /             SLAP.........95500
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /             SLAP.........95600
C                                                                        SLAP.........95700
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                    SLAP.........95800
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                  SLAP.........95900
C                                                                        SLAP.........96000
C     DATA SMALL(1), SMALL(2) /    128,      0 /                         SLAP.........96100
C     DATA SMALL(3), SMALL(4) /      0,      0 /                         SLAP.........96200
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /                         SLAP.........96300
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /                         SLAP.........96400
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /                         SLAP.........96500
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /                         SLAP.........96600
C     DATA DIVER(1), DIVER(2) /   9472,      0 /                         SLAP.........96700
C     DATA DIVER(3), DIVER(4) /      0,      0 /                         SLAP.........96800
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /                         SLAP.........96900
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /                         SLAP.........97000
C                                                                        SLAP.........97100
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /                       SLAP.........97200
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /                       SLAP.........97300
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /                       SLAP.........97400
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /                       SLAP.........97500
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /                       SLAP.........97600
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /                       SLAP.........97700
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /                       SLAP.........97800
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /                       SLAP.........97900
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /                       SLAP.........98000
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /                       SLAP.........98100
C                                                                        SLAP.........98200
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS                         SLAP.........98300
C                                                                        SLAP.........98400
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /               SLAP.........98500
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /               SLAP.........98600
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /               SLAP.........98700
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /               SLAP.........98800
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /               SLAP.........98900
C                                                                        SLAP.........99000
C     MACHINE CONSTANTS FOR THE SUN                                      SLAP.........99100
C                                                                        SLAP.........99200
C     DATA DMACH(1) / Z'0010000000000000' /                              SLAP.........99300
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /                              SLAP.........99400
C     DATA DMACH(3) / Z'3CA0000000000000' /                              SLAP.........99500
C     DATA DMACH(4) / Z'3CB0000000000000' /                              SLAP.........99600
C     DATA DMACH(5) / Z'3FD34413509F79FF' /                              SLAP.........99700
C                                                                        SLAP.........99800
C     MACHINE CONSTANTS FOR THE SUN                                      SLAP.........99900
C     USING THE -r8 COMPILER OPTION                                      SLAP........100000
C                                                                        SLAP........100100
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /              SLAP........100200
C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /              SLAP........100300
C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /              SLAP........100400
C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /              SLAP........100500
C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /              SLAP........100600
C                                                                        SLAP........100700
C     MACHINE CONSTANTS FOR THE SUN 386i                                 SLAP........100800
C                                                                        SLAP........100900
C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /               SLAP........101000
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /               SLAP........101100
C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /               SLAP........101200
C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'                 SLAP........101300
C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /               SLAP........101400
C                                                                        SLAP........101500
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER          SLAP........101600
C                                                                        SLAP........101700
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /           SLAP........101800
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /           SLAP........101900
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /           SLAP........102000
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /           SLAP........102100
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /           SLAP........102200
C                                                                        SLAP........102300
C.....GENERAL-PURPOSE VALUES INSERTED DURING INTEGRATION OF SLAP         SLAP........102400
C        WITH SUTRA. (ONLY SPECIFIED VALUES USED BY SLAP ROUTINES.       SLAP........102500
C        DMACH(2) IS NOT CRITICAL.)                                      SLAP........102600
C                                                                        SLAP........102700
      DATA DMACH(1) / 1.00D-200  /                                       SLAP........102800
      DATA DMACH(2) / 1.00D+200  /                                       SLAP........102900
      DATA DMACH(3) / 1.00D-16   /                                       SLAP........103000
C                                                                        SLAP........103100
C***FIRST EXECUTABLE STATEMENT  D1MACH                                   SLAP........103200
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH',       SLAP........103300
     +   'I OUT OF BOUNDS', 1, 2)                                        SLAP........103400
C                                                                        SLAP........103500
      D1MACH = DMACH(I)                                                  SLAP........103600
      RETURN                                                             SLAP........103700
C                                                                        SLAP........103800
      END                                                                SLAP........103900
*DECK DAXPY                                                              SLAP........104000
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)                       SLAP........104100
C***BEGIN PROLOGUE  DAXPY                                                SLAP........104200
C***PURPOSE  Compute a constant times a vector plus a vector.            SLAP........104300
C***LIBRARY   SLATEC (BLAS)                                              SLAP........104400
C***CATEGORY  D1A7                                                       SLAP........104500
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)               SLAP........104600
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR                        SLAP........104700
C***AUTHOR  Lawson, C. L., (JPL)                                         SLAP........104800
C           Hanson, R. J., (SNLA)                                        SLAP........104900
C           Kincaid, D. R., (U. of Texas)                                SLAP........105000
C           Krogh, F. T., (JPL)                                          SLAP........105100
C***DESCRIPTION                                                          SLAP........105200
C                                                                        SLAP........105300
C                B L A S  Subprogram                                     SLAP........105400
C    Description of Parameters                                           SLAP........105500
C                                                                        SLAP........105600
C     --Input--                                                          SLAP........105700
C        N  number of elements in input vector(s)                        SLAP........105800
C       DA  double precision scalar multiplier                           SLAP........105900
C       DX  double precision vector with N elements                      SLAP........106000
C     INCX  storage spacing between elements of DX                       SLAP........106100
C       DY  double precision vector with N elements                      SLAP........106200
C     INCY  storage spacing between elements of DY                       SLAP........106300
C                                                                        SLAP........106400
C     --Output--                                                         SLAP........106500
C       DY  double precision result (unchanged if N .LE. 0)              SLAP........106600
C                                                                        SLAP........106700
C     Overwrite double precision DY with double precision DA*DX + DY.    SLAP........106800
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +   SLAP........106900
C       DY(LY+I*INCY),                                                   SLAP........107000
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is     SLAP........107100
C     defined in a similar way using INCY.                               SLAP........107200
C                                                                        SLAP........107300
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.      SLAP........107400
C                 Krogh, Basic linear algebra subprograms for Fortran    SLAP........107500
C                 usage, Algorithm No. 539, Transactions on Mathematical SLAP........107600
C                 Software 5, 3 (September 1979), pp. 308-323.           SLAP........107700
C***ROUTINES CALLED  (NONE)                                              SLAP........107800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........107900
C   791001  DATE WRITTEN                                                 SLAP........108000
C   890831  Modified array declarations.  (WRB)                          SLAP........108100
C   890831  REVISION DATE from Version 3.2                               SLAP........108200
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........108300
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)            SLAP........108400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........108500
C***END PROLOGUE  DAXPY                                                  SLAP........108600
      DOUBLE PRECISION DX(*), DY(*), DA                                  SLAP........108700
C***FIRST EXECUTABLE STATEMENT  DAXPY                                    SLAP........108800
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN                                SLAP........108900
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60                            SLAP........109000
C                                                                        SLAP........109100
C     Code for unequal or nonpositive increments.                        SLAP........109200
C                                                                        SLAP........109300
    5 IX = 1                                                             SLAP........109400
      IY = 1                                                             SLAP........109500
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1                              SLAP........109600
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1                              SLAP........109700
      DO 10 I = 1,N                                                      SLAP........109800
        DY(IY) = DY(IY) + DA*DX(IX)                                      SLAP........109900
        IX = IX + INCX                                                   SLAP........110000
        IY = IY + INCY                                                   SLAP........110100
   10 CONTINUE                                                           SLAP........110200
      RETURN                                                             SLAP........110300
C                                                                        SLAP........110400
C     Code for both increments equal to 1.                               SLAP........110500
C                                                                        SLAP........110600
C     Clean-up loop so remaining vector length is a multiple of 4.       SLAP........110700
C                                                                        SLAP........110800
   20 M = MOD(N,4)                                                       SLAP........110900
      IF (M .EQ. 0) GO TO 40                                             SLAP........111000
      DO 30 I = 1,M                                                      SLAP........111100
        DY(I) = DY(I) + DA*DX(I)                                         SLAP........111200
   30 CONTINUE                                                           SLAP........111300
      IF (N .LT. 4) RETURN                                               SLAP........111400
   40 MP1 = M + 1                                                        SLAP........111500
      DO 50 I = MP1,N,4                                                  SLAP........111600
        DY(I) = DY(I) + DA*DX(I)                                         SLAP........111700
        DY(I+1) = DY(I+1) + DA*DX(I+1)                                   SLAP........111800
        DY(I+2) = DY(I+2) + DA*DX(I+2)                                   SLAP........111900
        DY(I+3) = DY(I+3) + DA*DX(I+3)                                   SLAP........112000
   50 CONTINUE                                                           SLAP........112100
      RETURN                                                             SLAP........112200
C                                                                        SLAP........112300
C     Code for equal, positive, non-unit increments.                     SLAP........112400
C                                                                        SLAP........112500
   60 NS = N*INCX                                                        SLAP........112600
      DO 70 I = 1,NS,INCX                                                SLAP........112700
        DY(I) = DA*DX(I) + DY(I)                                         SLAP........112800
   70 CONTINUE                                                           SLAP........112900
      RETURN                                                             SLAP........113000
      END                                                                SLAP........113100
*DECK DCG                                                                SLAP........113200
      SUBROUTINE DCG (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,    SLAP........113300
     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,   SLAP........113400
     +   IWORK)                                                          SLAP........113500
C***BEGIN PROLOGUE  DCG                                                  SLAP........113600
C***PURPOSE  Preconditioned Conjugate Gradient Sparse Ax=b Solver.       SLAP........113700
C            Routine to solve a symmetric positive definite linear       SLAP........113800
C            system  Ax = b  using the Preconditioned Conjugate          SLAP........113900
C            Gradient method.                                            SLAP........114000
C***LIBRARY   SLATEC (SLAP)                                              SLAP........114100
C***CATEGORY  D2B4                                                       SLAP........114200
C***TYPE      DOUBLE PRECISION (SCG-S, DCG-D)                            SLAP........114300
C***KEYWORDS  ITERATIVE PRECONDITION, SLAP, SPARSE,                      SLAP........114400
C             SYMMETRIC LINEAR SYSTEM                                    SLAP........114500
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........114600
C           Seager, Mark K., (LLNL)                                      SLAP........114700
C             Lawrence Livermore National Laboratory                     SLAP........114800
C             PO BOX 808, L-60                                           SLAP........114900
C             Livermore, CA 94550 (510) 423-3141                         SLAP........115000
C             seager@llnl.gov                                            SLAP........115100
C***DESCRIPTION                                                          SLAP........115200
C                                                                        SLAP........115300
C *Usage:                                                                SLAP........115400
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX            SLAP........115500
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)                    SLAP........115600
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)         SLAP........115700
C     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINED)                  SLAP........115800
C     EXTERNAL MATVEC, MSOLVE                                            SLAP........115900
C                                                                        SLAP........116000
C     CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,           SLAP........116100
C    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,        SLAP........116200
C    $     RWORK, IWORK )                                                SLAP........116300
C                                                                        SLAP........116400
C *Arguments:                                                            SLAP........116500
C N      :IN       Integer.                                              SLAP........116600
C         Order of the Matrix.                                           SLAP........116700
C B      :IN       Double Precision B(N).                                SLAP........116800
C         Right-hand side vector.                                        SLAP........116900
C X      :INOUT    Double Precision X(N).                                SLAP........117000
C         On input X is your initial guess for solution vector.          SLAP........117100
C         On output X is the final approximate solution.                 SLAP........117200
C NELT   :IN       Integer.                                              SLAP........117300
C         Number of Non-Zeros stored in A.                               SLAP........117400
C IA     :IN       Integer IA(NELT).                                     SLAP........117500
C JA     :IN       Integer JA(NELT).                                     SLAP........117600
C A      :IN       Double Precision A(NELT).                             SLAP........117700
C         These arrays contain the matrix data structure for A.          SLAP........117800
C         It could take any form.  See "Description", below,             SLAP........117900
C         for more details.                                              SLAP........118000
C ISYM   :IN       Integer.                                              SLAP........118100
C         Flag to indicate symmetric storage format.                     SLAP........118200
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........118300
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........118400
C         or lower triangle of the matrix is stored.                     SLAP........118500
C MATVEC :EXT      External.                                             SLAP........118600
C         Name of a routine which performs the matrix vector multiply    SLAP........118700
C         Y = A*X given A and X.  The name of the MATVEC routine must    SLAP........118800
C         be declared external in the calling program.  The calling      SLAP........118900
C         sequence to MATVEC is:                                         SLAP........119000
C                                                                        SLAP........119100
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )              SLAP........119200
C                                                                        SLAP........119300
C         Where N is the number of unknowns, Y is the product A*X        SLAP........119400
C         upon return X is an input vector, NELT is the number of        SLAP........119500
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.      SLAP........119600
C         ISYM is a flag which, if non-zero, denotest that A is          SLAP........119700
C         symmetric and only the lower or upper triangle is stored.      SLAP........119800
C MSOLVE :EXT      External.                                             SLAP........119900
C         Name of a routine which solves a linear system MZ = R for      SLAP........120000
C         Z given R with the preconditioning matrix M (M is supplied via SLAP........120100
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  SLAP........120200
C         be declared external in the calling program.  The calling      SLAP........120300
C         sequence to MSOLVE is:                                         SLAP........120400
C                                                                        SLAP........120500
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........120600
C                                                                        SLAP........120700
C         Where N is the number of unknowns, R is the right-hand side    SLAP........120800
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........120900
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........121000
C         that can be used to pass necessary preconditioning information SLAP........121100
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........121200
C         for the same purpose as RWORK.                                 SLAP........121300
C ITOL   :IN       Integer.                                              SLAP........121400
C         Flag to indicate type of convergence criterion.                SLAP........121500
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........121600
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........121700
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........121800
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........121900
C         side is less than TOL, where M-inv is the inverse of the       SLAP........122000
C         diagonal of A.                                                 SLAP........122100
C         ITOL=11 is often useful for checking and comparing different   SLAP........122200
C         routines.  For this case, the user must supply the "exact"     SLAP........122300
C         solution or a very accurate approximation (one with an error   SLAP........122400
C         much less than TOL) through a common block,                    SLAP........122500
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........122600
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........122700
C         between the iterative approximation and the user-supplied      SLAP........122800
C         solution divided by the 2-norm of the user-supplied solution   SLAP........122900
C         is less than TOL.  Note that this requires the user to set up  SLAP........123000
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........123100
C         The routine with this declaration should be loaded before the  SLAP........123200
C         stop test so that the correct length is used by the loader.    SLAP........123300
C         This procedure is not standard Fortran and may not work        SLAP........123400
C         correctly on your system (although it has worked on every      SLAP........123500
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........123600
C         common block is indeed standard Fortran.                       SLAP........123700
C TOL    :INOUT    Double Precision.                                     SLAP........123800
C         Convergence criterion, as described above.  (Reset if IERR=4.) SLAP........123900
C ITMAX  :IN       Integer.                                              SLAP........124000
C         Maximum number of iterations.                                  SLAP........124100
C ITER   :OUT      Integer.                                              SLAP........124200
C         Number of iterations required to reach convergence, or         SLAP........124300
C         ITMAX+1 if convergence criterion could not be achieved in      SLAP........124400
C         ITMAX iterations.                                              SLAP........124500
C ERR    :OUT      Double Precision.                                     SLAP........124600
C         Error estimate of error in final approximate solution, as      SLAP........124700
C         defined by ITOL.                                               SLAP........124800
C IERR   :OUT      Integer.                                              SLAP........124900
C         Return error flag.                                             SLAP........125000
C           IERR = 0 => All went well.                                   SLAP........125100
C           IERR = 1 => Insufficient space allocated for WORK or IWORK.  SLAP........125200
C           IERR = 2 => Method failed to converge in ITMAX steps.        SLAP........125300
C           IERR = 3 => Error in user input.                             SLAP........125400
C                       Check input values of N, ITOL.                   SLAP........125500
C           IERR = 4 => User error tolerance set too tight.              SLAP........125600
C                       Reset to 500*D1MACH(3).  Iteration proceeded.    SLAP........125700
C           IERR = 5 => Preconditioning matrix, M, is not positive       SLAP........125800
C                       definite.  (r,z) < 0.                            SLAP........125900
C           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.  SLAP........126000
C IUNIT  :IN       Integer.                                              SLAP........126100
C         Unit number on which to write the error at each iteration,     SLAP........126200
C         if this is desired for monitoring convergence.  If unit        SLAP........126300
C         number is 0, no writing will occur.                            SLAP........126400
C R      :WORK     Double Precision R(N).                                SLAP........126500
C Z      :WORK     Double Precision Z(N).                                SLAP........126600
C P      :WORK     Double Precision P(N).                                SLAP........126700
C DZ     :WORK     Double Precision DZ(N).                               SLAP........126800
C         Double Precision arrays used for workspace.                    SLAP........126900
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........127000
C         Double Precision array that can be used by  MSOLVE.            SLAP........127100
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........127200
C         Integer array that can be used by  MSOLVE.                     SLAP........127300
C                                                                        SLAP........127400
C *Description                                                           SLAP........127500
C       This routine does  not care  what matrix data   structure is     SLAP........127600
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE     SLAP........127700
C       routines, with  the arguments as  described above.  The user     SLAP........127800
C       could write any type of structure and the appropriate MATVEC     SLAP........127900
C       and MSOLVE routines.  It is assumed  that A is stored in the     SLAP........128000
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is     SLAP........128100
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP     SLAP........128200
C       routines DSDCG and DSICCG are examples of this procedure.        SLAP........128300
C                                                                        SLAP........128400
C       Two  examples  of  matrix  data structures  are the: 1) SLAP     SLAP........128500
C       Triad  format and 2) SLAP Column format.                         SLAP........128600
C                                                                        SLAP........128700
C       =================== S L A P Triad format ===================     SLAP........128800
C                                                                        SLAP........128900
C       In  this   format only the  non-zeros are  stored.  They may     SLAP........129000
C       appear  in *ANY* order.   The user  supplies three arrays of     SLAP........129100
C       length NELT, where  NELT  is the number  of non-zeros in the     SLAP........129200
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero     SLAP........129300
C       the  user puts   the row  and  column index   of that matrix     SLAP........129400
C       element in the IA and JA arrays.  The  value of the non-zero     SLAP........129500
C       matrix  element is  placed in  the corresponding location of     SLAP........129600
C       the A  array.  This is  an extremely easy data  structure to     SLAP........129700
C       generate.  On  the other hand it  is  not too  efficient  on     SLAP........129800
C       vector  computers   for the  iterative  solution  of  linear     SLAP........129900
C       systems.  Hence, SLAP  changes this input  data structure to     SLAP........130000
C       the SLAP   Column  format for the  iteration (but   does not     SLAP........130100
C       change it back).                                                 SLAP........130200
C                                                                        SLAP........130300
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........130400
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........130500
C                                                                        SLAP........130600
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........130700
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........130800
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........130900
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........131000
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........131100
C       | 0  0  0 44  0|                                                 SLAP........131200
C       |51  0 53  0 55|                                                 SLAP........131300
C                                                                        SLAP........131400
C       =================== S L A P Column format ==================     SLAP........131500
C                                                                        SLAP........131600
C       In  this format   the non-zeros are    stored counting  down     SLAP........131700
C       columns (except  for the diagonal  entry, which must  appear     SLAP........131800
C       first  in each "column") and are  stored in the  double pre-     SLAP........131900
C       cision array  A. In  other  words,  for each  column  in the     SLAP........132000
C       matrix  first put  the diagonal entry in A.  Then put in the     SLAP........132100
C       other non-zero  elements going  down the column  (except the     SLAP........132200
C       diagonal)  in order.  The IA array  holds the  row index for     SLAP........132300
C       each non-zero.  The JA array  holds the offsets into the IA,     SLAP........132400
C       A  arrays  for  the  beginning  of  each  column.  That  is,     SLAP........132500
C       IA(JA(ICOL)),A(JA(ICOL)) are the first elements of the ICOL-     SLAP........132600
C       th column in IA and A, and IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)     SLAP........132700
C       are  the last elements of the ICOL-th column.   Note that we     SLAP........132800
C       always have JA(N+1)=NELT+1, where N is the number of columns     SLAP........132900
C       in the matrix  and NELT  is the number  of non-zeros  in the     SLAP........133000
C       matrix.                                                          SLAP........133100
C                                                                        SLAP........133200
C       Here is an example of the  SLAP Column  storage format for a     SLAP........133300
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........133400
C       column):                                                         SLAP........133500
C                                                                        SLAP........133600
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........133700
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........133800
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........133900
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........134000
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........134100
C       | 0  0  0 44  0|                                                 SLAP........134200
C       |51  0 53  0 55|                                                 SLAP........134300
C                                                                        SLAP........134400
C *Cautions:                                                             SLAP........134500
C     This routine will attempt to write to the Fortran logical output   SLAP........134600
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........134700
C     this logical unit is attached to a file or terminal before calling SLAP........134800
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........134900
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........135000
C                                                                        SLAP........135100
C***SEE ALSO  DSDCG, DSICCG                                              SLAP........135200
C***REFERENCES  1. Louis Hageman and David Young, Applied Iterative      SLAP........135300
C                  Methods, Academic Press, New York, 1981.              SLAP........135400
C               2. Concus, Golub and O'Leary, A Generalized Conjugate    SLAP........135500
C                  Gradient Method for the Numerical Solution of         SLAP........135600
C                  Elliptic Partial Differential Equations, in Sparse    SLAP........135700
C                  Matrix Computations, Bunch and Rose, Eds., Academic   SLAP........135800
C                  Press, New York, 1979.                                SLAP........135900
C               3. Mark K. Seager, A SLAP for the Masses, in             SLAP........136000
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,   SLAP........136100
C                  Algorithms and Applications, Wiley, 1989, pp.135-155. SLAP........136200
C***ROUTINES CALLED  D1MACH, DAXPY, DCOPY, DDOT, ISDCG                   SLAP........136300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........136400
C   890404  DATE WRITTEN                                                 SLAP........136500
C   890404  Previous REVISION DATE                                       SLAP........136600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........136700
C   890921  Removed TeX from comments.  (FNF)                            SLAP........136800
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........136900
C           standard.  (FNF)                                             SLAP........137000
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........137100
C   891004  Added new reference.                                         SLAP........137200
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........137300
C   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)  SLAP........137400
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........137500
C   920511  Added complete declaration section.  (WRB)                   SLAP........137600
C   920929  Corrected format of references.  (FNF)                       SLAP........137700
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)     SLAP........137800
C***END PROLOGUE  DCG                                                    SLAP........137900
C     .. Scalar Arguments ..                                             SLAP........138000
      DOUBLE PRECISION ERR, TOL                                          SLAP........138100
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT              SLAP........138200
C     .. Array Arguments ..                                              SLAP........138300
      DOUBLE PRECISION A(NELT), B(N), DZ(N), P(N), R(N), RWORK(*), X(N), SLAP........138400
     +                 Z(N)                                              SLAP........138500
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........138600
C     .. Subroutine Arguments ..                                         SLAP........138700
      EXTERNAL MATVEC, MSOLVE                                            SLAP........138800
C     .. Local Scalars ..                                                SLAP........138900
      DOUBLE PRECISION AK, AKDEN, BK, BKDEN, BKNUM, BNRM, SOLNRM, TOLMIN SLAP........139000
      INTEGER I, K                                                       SLAP........139100
C     .. External Functions ..                                           SLAP........139200
      DOUBLE PRECISION D1MACH, DDOT                                      SLAP........139300
      INTEGER ISDCG                                                      SLAP........139400
      EXTERNAL D1MACH, DDOT, ISDCG                                       SLAP........139500
C     .. External Subroutines ..                                         SLAP........139600
      EXTERNAL DAXPY, DCOPY                                              SLAP........139700
C***FIRST EXECUTABLE STATEMENT  DCG                                      SLAP........139800
C                                                                        SLAP........139900
C         Check some of the input data.                                  SLAP........140000
C                                                                        SLAP........140100
      ITER = 0                                                           SLAP........140200
      IERR = 0                                                           SLAP........140300
      IF( N.LT.1 ) THEN                                                  SLAP........140400
         IERR = 3                                                        SLAP........140500
         RETURN                                                          SLAP........140600
      ENDIF                                                              SLAP........140700
      TOLMIN = 500*D1MACH(3)                                             SLAP........140800
      IF( TOL.LT.TOLMIN ) THEN                                           SLAP........140900
         TOL = TOLMIN                                                    SLAP........141000
         IERR = 4                                                        SLAP........141100
      ENDIF                                                              SLAP........141200
C                                                                        SLAP........141300
C         Calculate initial residual and pseudo-residual, and check      SLAP........141400
C         stopping criterion.                                            SLAP........141500
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)                        SLAP........141600
      DO 10 I = 1, N                                                     SLAP........141700
         R(I) = B(I) - R(I)                                              SLAP........141800
 10   CONTINUE                                                           SLAP........141900
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)          SLAP........142000
C                                                                        SLAP........142100
      IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,       SLAP........142200
     $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,                   SLAP........142300
     $     RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200        SLAP........142400
      IF( IERR.NE.0 ) RETURN                                             SLAP........142500
C                                                                        SLAP........142600
C         ***** Iteration loop *****                                     SLAP........142700
C                                                                        SLAP........142800
      DO 100 K=1,ITMAX                                                   SLAP........142900
         ITER = K                                                        SLAP........143000
C                                                                        SLAP........143100
C         Calculate coefficient bk and direction vector p.               SLAP........143200
         BKNUM = DDOT(N, Z, 1, R, 1)                                     SLAP........143300
         IF( BKNUM.LE.0.0D0 ) THEN                                       SLAP........143400
            IERR = 5                                                     SLAP........143500
            RETURN                                                       SLAP........143600
         ENDIF                                                           SLAP........143700
         IF(ITER .EQ. 1) THEN                                            SLAP........143800
            CALL DCOPY(N, Z, 1, P, 1)                                    SLAP........143900
         ELSE                                                            SLAP........144000
            BK = BKNUM/BKDEN                                             SLAP........144100
            DO 20 I = 1, N                                               SLAP........144200
               P(I) = Z(I) + BK*P(I)                                     SLAP........144300
 20         CONTINUE                                                     SLAP........144400
         ENDIF                                                           SLAP........144500
         BKDEN = BKNUM                                                   SLAP........144600
C                                                                        SLAP........144700
C         Calculate coefficient ak, new iterate x, new residual r,       SLAP........144800
C         and new pseudo-residual z.                                     SLAP........144900
         CALL MATVEC(N, P, Z, NELT, IA, JA, A, ISYM)                     SLAP........145000
         AKDEN = DDOT(N, P, 1, Z, 1)                                     SLAP........145100
         IF( AKDEN.LE.0.0D0 ) THEN                                       SLAP........145200
            IERR = 6                                                     SLAP........145300
            RETURN                                                       SLAP........145400
         ENDIF                                                           SLAP........145500
         AK = BKNUM/AKDEN                                                SLAP........145600
         CALL DAXPY(N, AK, P, 1, X, 1)                                   SLAP........145700
         CALL DAXPY(N, -AK, Z, 1, R, 1)                                  SLAP........145800
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)       SLAP........145900
C                                                                        SLAP........146000
C         check stopping criterion.                                      SLAP........146100
         IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,    SLAP........146200
     $        ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,         SLAP........146300
     $        IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200            SLAP........146400
C                                                                        SLAP........146500
 100  CONTINUE                                                           SLAP........146600
C                                                                        SLAP........146700
C         *****   end of loop  *****                                     SLAP........146800
C                                                                        SLAP........146900
C         stopping criterion not satisfied.                              SLAP........147000
      ITER = ITMAX + 1                                                   SLAP........147100
      IERR = 2                                                           SLAP........147200
C                                                                        SLAP........147300
 200  RETURN                                                             SLAP........147400
C------------- LAST LINE OF DCG FOLLOWS -----------------------------    SLAP........147500
      END                                                                SLAP........147600
*DECK DCHKW                                                              SLAP........147700
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR) SLAP........147800
C***BEGIN PROLOGUE  DCHKW                                                SLAP........147900
C***SUBSIDIARY                                                           SLAP........148000
C***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.                       SLAP........148100
C            This routine checks the work array lengths and interfaces   SLAP........148200
C            to the SLATEC error handler if a problem is found.          SLAP........148300
C***LIBRARY   SLATEC (SLAP)                                              SLAP........148400
C***CATEGORY  R2                                                         SLAP........148500
C***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D)                        SLAP........148600
C***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING                   SLAP........148700
C***AUTHOR  Seager, Mark K., (LLNL)                                      SLAP........148800
C             Lawrence Livermore National Laboratory                     SLAP........148900
C             PO BOX 808, L-60                                           SLAP........149000
C             Livermore, CA 94550 (510) 423-3141                         SLAP........149100
C             seager@llnl.gov                                            SLAP........149200
C***DESCRIPTION                                                          SLAP........149300
C                                                                        SLAP........149400
C *Usage:                                                                SLAP........149500
C     CHARACTER*(*) NAME                                                 SLAP........149600
C     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER                       SLAP........149700
C     DOUBLE PRECISION ERR                                               SLAP........149800
C                                                                        SLAP........149900
C     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )      SLAP........150000
C                                                                        SLAP........150100
C *Arguments:                                                            SLAP........150200
C NAME   :IN       Character*(*).                                        SLAP........150300
C         Name of the calling routine.  This is used in the output       SLAP........150400
C         message, if an error is detected.                              SLAP........150500
C LOCIW  :IN       Integer.                                              SLAP........150600
C         Location of the first free element in the integer workspace    SLAP........150700
C         array.                                                         SLAP........150800
C LENIW  :IN       Integer.                                              SLAP........150900
C         Length of the integer workspace array.                         SLAP........151000
C LOCW   :IN       Integer.                                              SLAP........151100
C         Location of the first free element in the double precision     SLAP........151200
C         workspace array.                                               SLAP........151300
C LENRW  :IN       Integer.                                              SLAP........151400
C         Length of the double precision workspace array.                SLAP........151500
C IERR   :OUT      Integer.                                              SLAP........151600
C         Return error flag.                                             SLAP........151700
C               IERR = 0 => All went well.                               SLAP........151800
C               IERR = 1 => Insufficient storage allocated for           SLAP........151900
C                           WORK or IWORK.                               SLAP........152000
C ITER   :OUT      Integer.                                              SLAP........152100
C         Set to zero on return.                                         SLAP........152200
C ERR    :OUT      Double Precision.                                     SLAP........152300
C         Set to the smallest positive magnitude if all went well.       SLAP........152400
C         Set to a very large number if an error is detected.            SLAP........152500
C                                                                        SLAP........152600
C***REFERENCES  (NONE)                                                   SLAP........152700
C***ROUTINES CALLED  D1MACH, XERMSG                                      SLAP........152800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........152900
C   880225  DATE WRITTEN                                                 SLAP........153000
C   881213  Previous REVISION DATE                                       SLAP........153100
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........153200
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........153300
C           standard.  (FNF)                                             SLAP........153400
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........153500
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)              SLAP........153600
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........153700
C   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI      SLAP........153800
C           X3.9-1978.  (FNF)                                            SLAP........153900
C   910506  Made subsidiary.  (FNF)                                      SLAP........154000
C   920511  Added complete declaration section.  (WRB)                   SLAP........154100
C   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)    SLAP........154200
C***END PROLOGUE  DCHKW                                                  SLAP........154300
C     .. Scalar Arguments ..                                             SLAP........154400
      DOUBLE PRECISION ERR                                               SLAP........154500
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW                       SLAP........154600
      CHARACTER NAME*(*)                                                 SLAP........154700
C     .. Local Scalars ..                                                SLAP........154800
      CHARACTER XERN1*8, XERN2*8, XERNAM*8                               SLAP........154900
C     .. External Functions ..                                           SLAP........155000
      DOUBLE PRECISION D1MACH                                            SLAP........155100
      EXTERNAL D1MACH                                                    SLAP........155200
C     .. External Subroutines ..                                         SLAP........155300
      EXTERNAL XERMSG                                                    SLAP........155400
C***FIRST EXECUTABLE STATEMENT  DCHKW                                    SLAP........155500
C                                                                        SLAP........155600
C         Check the Integer workspace situation.                         SLAP........155700
C                                                                        SLAP........155800
      IERR = 0                                                           SLAP........155900
      ITER = 0                                                           SLAP........156000
      ERR = D1MACH(1)                                                    SLAP........156100
      IF( LOCIW.GT.LENIW ) THEN                                          SLAP........156200
         IERR = 1                                                        SLAP........156300
         ERR = D1MACH(2)                                                 SLAP........156400
         XERNAM = NAME                                                   SLAP........156500
         WRITE (XERN1, '(I8)') LOCIW                                     SLAP........156600
         WRITE (XERN2, '(I8)') LENIW                                     SLAP........156700
         CALL XERMSG ('SLATEC', 'DCHKW',                                 SLAP........156800
     $      'In ' // XERNAM // ', INTEGER work array too short.  ' //    SLAP........156900
     $      'IWORK needs ' // XERN1 // '; have allocated ' // XERN2,     SLAP........157000
     $      1, 1)                                                        SLAP........157100
      ENDIF                                                              SLAP........157200
C                                                                        SLAP........157300
C         Check the Double Precision workspace situation.                SLAP........157400
      IF( LOCW.GT.LENW ) THEN                                            SLAP........157500
         IERR = 1                                                        SLAP........157600
         ERR = D1MACH(2)                                                 SLAP........157700
         XERNAM = NAME                                                   SLAP........157800
         WRITE (XERN1, '(I8)') LOCW                                      SLAP........157900
         WRITE (XERN2, '(I8)') LENW                                      SLAP........158000
         CALL XERMSG ('SLATEC', 'DCHKW',                                 SLAP........158100
     $      'In ' // XERNAM // ', DOUBLE PRECISION work array too ' //   SLAP........158200
     $      'short.  RWORK needs ' // XERN1 // '; have allocated ' //    SLAP........158300
     $      XERN2, 1, 1)                                                 SLAP........158400
      ENDIF                                                              SLAP........158500
      RETURN                                                             SLAP........158600
C------------- LAST LINE OF DCHKW FOLLOWS ----------------------------   SLAP........158700
      END                                                                SLAP........158800
*DECK DCOPY                                                              SLAP........158900
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)                           SLAP........159000
C***BEGIN PROLOGUE  DCOPY                                                SLAP........159100
C***PURPOSE  Copy a vector.                                              SLAP........159200
C***LIBRARY   SLATEC (BLAS)                                              SLAP........159300
C***CATEGORY  D1A5                                                       SLAP........159400
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)      SLAP........159500
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR                         SLAP........159600
C***AUTHOR  Lawson, C. L., (JPL)                                         SLAP........159700
C           Hanson, R. J., (SNLA)                                        SLAP........159800
C           Kincaid, D. R., (U. of Texas)                                SLAP........159900
C           Krogh, F. T., (JPL)                                          SLAP........160000
C***DESCRIPTION                                                          SLAP........160100
C                                                                        SLAP........160200
C                B L A S  Subprogram                                     SLAP........160300
C    Description of Parameters                                           SLAP........160400
C                                                                        SLAP........160500
C     --Input--                                                          SLAP........160600
C        N  number of elements in input vector(s)                        SLAP........160700
C       DX  double precision vector with N elements                      SLAP........160800
C     INCX  storage spacing between elements of DX                       SLAP........160900
C       DY  double precision vector with N elements                      SLAP........161000
C     INCY  storage spacing between elements of DY                       SLAP........161100
C                                                                        SLAP........161200
C     --Output--                                                         SLAP........161300
C       DY  copy of vector DX (unchanged if N .LE. 0)                    SLAP........161400
C                                                                        SLAP........161500
C     Copy double precision DX to double precision DY.                   SLAP........161600
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),             SLAP........161700
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is     SLAP........161800
C     defined in a similar way using INCY.                               SLAP........161900
C                                                                        SLAP........162000
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.      SLAP........162100
C                 Krogh, Basic linear algebra subprograms for Fortran    SLAP........162200
C                 usage, Algorithm No. 539, Transactions on Mathematical SLAP........162300
C                 Software 5, 3 (September 1979), pp. 308-323.           SLAP........162400
C***ROUTINES CALLED  (NONE)                                              SLAP........162500
C***REVISION HISTORY  (YYMMDD)                                           SLAP........162600
C   791001  DATE WRITTEN                                                 SLAP........162700
C   890831  Modified array declarations.  (WRB)                          SLAP........162800
C   890831  REVISION DATE from Version 3.2                               SLAP........162900
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........163000
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)            SLAP........163100
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........163200
C***END PROLOGUE  DCOPY                                                  SLAP........163300
      DOUBLE PRECISION DX(*), DY(*)                                      SLAP........163400
C***FIRST EXECUTABLE STATEMENT  DCOPY                                    SLAP........163500
      IF (N .LE. 0) RETURN                                               SLAP........163600
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60                            SLAP........163700
C                                                                        SLAP........163800
C     Code for unequal or nonpositive increments.                        SLAP........163900
C                                                                        SLAP........164000
    5 IX = 1                                                             SLAP........164100
      IY = 1                                                             SLAP........164200
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1                              SLAP........164300
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1                              SLAP........164400
      DO 10 I = 1,N                                                      SLAP........164500
        DY(IY) = DX(IX)                                                  SLAP........164600
        IX = IX + INCX                                                   SLAP........164700
        IY = IY + INCY                                                   SLAP........164800
   10 CONTINUE                                                           SLAP........164900
      RETURN                                                             SLAP........165000
C                                                                        SLAP........165100
C     Code for both increments equal to 1.                               SLAP........165200
C                                                                        SLAP........165300
C     Clean-up loop so remaining vector length is a multiple of 7.       SLAP........165400
C                                                                        SLAP........165500
   20 M = MOD(N,7)                                                       SLAP........165600
      IF (M .EQ. 0) GO TO 40                                             SLAP........165700
      DO 30 I = 1,M                                                      SLAP........165800
        DY(I) = DX(I)                                                    SLAP........165900
   30 CONTINUE                                                           SLAP........166000
      IF (N .LT. 7) RETURN                                               SLAP........166100
   40 MP1 = M + 1                                                        SLAP........166200
      DO 50 I = MP1,N,7                                                  SLAP........166300
        DY(I) = DX(I)                                                    SLAP........166400
        DY(I+1) = DX(I+1)                                                SLAP........166500
        DY(I+2) = DX(I+2)                                                SLAP........166600
        DY(I+3) = DX(I+3)                                                SLAP........166700
        DY(I+4) = DX(I+4)                                                SLAP........166800
        DY(I+5) = DX(I+5)                                                SLAP........166900
        DY(I+6) = DX(I+6)                                                SLAP........167000
   50 CONTINUE                                                           SLAP........167100
      RETURN                                                             SLAP........167200
C                                                                        SLAP........167300
C     Code for equal, positive, non-unit increments.                     SLAP........167400
C                                                                        SLAP........167500
   60 NS = N*INCX                                                        SLAP........167600
      DO 70 I = 1,NS,INCX                                                SLAP........167700
        DY(I) = DX(I)                                                    SLAP........167800
   70 CONTINUE                                                           SLAP........167900
      RETURN                                                             SLAP........168000
      END                                                                SLAP........168100
*DECK DDOT                                                               SLAP........168200
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)             SLAP........168300
C***BEGIN PROLOGUE  DDOT                                                 SLAP........168400
C***PURPOSE  Compute the inner product of two vectors.                   SLAP........168500
C***LIBRARY   SLATEC (BLAS)                                              SLAP........168600
C***CATEGORY  D1A4                                                       SLAP........168700
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)                 SLAP........168800
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR                SLAP........168900
C***AUTHOR  Lawson, C. L., (JPL)                                         SLAP........169000
C           Hanson, R. J., (SNLA)                                        SLAP........169100
C           Kincaid, D. R., (U. of Texas)                                SLAP........169200
C           Krogh, F. T., (JPL)                                          SLAP........169300
C***DESCRIPTION                                                          SLAP........169400
C                                                                        SLAP........169500
C                B L A S  Subprogram                                     SLAP........169600
C    Description of Parameters                                           SLAP........169700
C                                                                        SLAP........169800
C     --Input--                                                          SLAP........169900
C        N  number of elements in input vector(s)                        SLAP........170000
C       DX  double precision vector with N elements                      SLAP........170100
C     INCX  storage spacing between elements of DX                       SLAP........170200
C       DY  double precision vector with N elements                      SLAP........170300
C     INCY  storage spacing between elements of DY                       SLAP........170400
C                                                                        SLAP........170500
C     --Output--                                                         SLAP........170600
C     DDOT  double precision dot product (zero if N .LE. 0)              SLAP........170700
C                                                                        SLAP........170800
C     Returns the dot product of double precision DX and DY.             SLAP........170900
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),     SLAP........171000
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is     SLAP........171100
C     defined in a similar way using INCY.                               SLAP........171200
C                                                                        SLAP........171300
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.      SLAP........171400
C                 Krogh, Basic linear algebra subprograms for Fortran    SLAP........171500
C                 usage, Algorithm No. 539, Transactions on Mathematical SLAP........171600
C                 Software 5, 3 (September 1979), pp. 308-323.           SLAP........171700
C***ROUTINES CALLED  (NONE)                                              SLAP........171800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........171900
C   791001  DATE WRITTEN                                                 SLAP........172000
C   890831  Modified array declarations.  (WRB)                          SLAP........172100
C   890831  REVISION DATE from Version 3.2                               SLAP........172200
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........172300
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)            SLAP........172400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........172500
C***END PROLOGUE  DDOT                                                   SLAP........172600
      DOUBLE PRECISION DX(*), DY(*)                                      SLAP........172700
C***FIRST EXECUTABLE STATEMENT  DDOT                                     SLAP........172800
      DDOT = 0.0D0                                                       SLAP........172900
      IF (N .LE. 0) RETURN                                               SLAP........173000
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60                            SLAP........173100
C                                                                        SLAP........173200
C     Code for unequal or nonpositive increments.                        SLAP........173300
C                                                                        SLAP........173400
    5 IX = 1                                                             SLAP........173500
      IY = 1                                                             SLAP........173600
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1                              SLAP........173700
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1                              SLAP........173800
      DO 10 I = 1,N                                                      SLAP........173900
        DDOT = DDOT + DX(IX)*DY(IY)                                      SLAP........174000
        IX = IX + INCX                                                   SLAP........174100
        IY = IY + INCY                                                   SLAP........174200
   10 CONTINUE                                                           SLAP........174300
      RETURN                                                             SLAP........174400
C                                                                        SLAP........174500
C     Code for both increments equal to 1.                               SLAP........174600
C                                                                        SLAP........174700
C     Clean-up loop so remaining vector length is a multiple of 5.       SLAP........174800
C                                                                        SLAP........174900
   20 M = MOD(N,5)                                                       SLAP........175000
      IF (M .EQ. 0) GO TO 40                                             SLAP........175100
      DO 30 I = 1,M                                                      SLAP........175200
         DDOT = DDOT + DX(I)*DY(I)                                       SLAP........175300
   30 CONTINUE                                                           SLAP........175400
      IF (N .LT. 5) RETURN                                               SLAP........175500
   40 MP1 = M + 1                                                        SLAP........175600
      DO 50 I = MP1,N,5                                                  SLAP........175700
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +    SLAP........175800
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)                    SLAP........175900
   50 CONTINUE                                                           SLAP........176000
      RETURN                                                             SLAP........176100
C                                                                        SLAP........176200
C     Code for equal, positive, non-unit increments.                     SLAP........176300
C                                                                        SLAP........176400
   60 NS = N*INCX                                                        SLAP........176500
      DO 70 I = 1,NS,INCX                                                SLAP........176600
        DDOT = DDOT + DX(I)*DY(I)                                        SLAP........176700
   70 CONTINUE                                                           SLAP........176800
      RETURN                                                             SLAP........176900
      END                                                                SLAP........177000
*DECK DGMRES                                                             SLAP........177100
      SUBROUTINE DGMRES (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, SLAP........177200
     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, RGWK, LRGW,   SLAP........177300
     +   IGWK, LIGW, RWORK, IWORK)                                       SLAP........177400
C***BEGIN PROLOGUE  DGMRES                                               SLAP........177500
C***PURPOSE  Preconditioned GMRES iterative sparse Ax=b solver.          SLAP........177600
C            This routine uses the generalized minimum residual          SLAP........177700
C            (GMRES) method with preconditioning to solve                SLAP........177800
C            non-symmetric linear systems of the form: Ax = b.           SLAP........177900
C***LIBRARY   SLATEC (SLAP)                                              SLAP........178000
C***CATEGORY  D2A4, D2B4                                                 SLAP........178100
C***TYPE      DOUBLE PRECISION (SGMRES-S, DGMRES-D)                      SLAP........178200
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........178300
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........178400
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........178500
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........178600
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........178700
C             Lawrence Livermore National Laboratory                     SLAP........178800
C             PO Box 808, L-60                                           SLAP........178900
C             Livermore, CA 94550 (510) 423-3141                         SLAP........179000
C***DESCRIPTION                                                          SLAP........179100
C                                                                        SLAP........179200
C *Usage:                                                                SLAP........179300
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX          SLAP........179400
C      INTEGER   ITER, IERR, IUNIT, LRGW, IGWK(LIGW), LIGW               SLAP........179500
C      INTEGER   IWORK(USER DEFINED)                                     SLAP........179600
C      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N)      SLAP........179700
C      DOUBLE PRECISION RGWK(LRGW), RWORK(USER DEFINED)                  SLAP........179800
C      EXTERNAL  MATVEC, MSOLVE                                          SLAP........179900
C                                                                        SLAP........180000
C      CALL DGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,       SLAP........180100
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX,            SLAP........180200
C     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK)                        SLAP........180300
C                                                                        SLAP........180400
C *Arguments:                                                            SLAP........180500
C N      :IN       Integer.                                              SLAP........180600
C         Order of the Matrix.                                           SLAP........180700
C B      :IN       Double Precision B(N).                                SLAP........180800
C         Right-hand side vector.                                        SLAP........180900
C X      :INOUT    Double Precision X(N).                                SLAP........181000
C         On input X is your initial guess for the solution vector.      SLAP........181100
C         On output X is the final approximate solution.                 SLAP........181200
C NELT   :IN       Integer.                                              SLAP........181300
C         Number of Non-Zeros stored in A.                               SLAP........181400
C IA     :IN       Integer IA(NELT).                                     SLAP........181500
C JA     :IN       Integer JA(NELT).                                     SLAP........181600
C A      :IN       Double Precision A(NELT).                             SLAP........181700
C         These arrays contain the matrix data structure for A.          SLAP........181800
C         It could take any form.  See "Description", below,             SLAP........181900
C         for more details.                                              SLAP........182000
C ISYM   :IN       Integer.                                              SLAP........182100
C         Flag to indicate symmetric storage format.                     SLAP........182200
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........182300
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........182400
C         or lower triangle of the matrix is stored.                     SLAP........182500
C MATVEC :EXT      External.                                             SLAP........182600
C         Name of a routine which performs the matrix vector multiply    SLAP........182700
C         Y = A*X given A and X.  The name of the MATVEC routine must    SLAP........182800
C         be declared external in the calling program.  The calling      SLAP........182900
C         sequence to MATVEC is:                                         SLAP........183000
C             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)                SLAP........183100
C         where N is the number of unknowns, Y is the product A*X        SLAP........183200
C         upon return, X is an input vector, and NELT is the number of   SLAP........183300
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.      SLAP........183400
C         ISYM is a flag which, if non-zero, denotes that A is           SLAP........183500
C         symmetric and only the lower or upper triangle is stored.      SLAP........183600
C MSOLVE :EXT      External.                                             SLAP........183700
C         Name of the routine which solves a linear system Mz = r for    SLAP........183800
C         z given r with the preconditioning matrix M (M is supplied via SLAP........183900
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must   SLAP........184000
C         be declared external in the calling program.  The calling      SLAP........184100
C         sequence to MSOLVE is:                                         SLAP........184200
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........184300
C         Where N is the number of unknowns, R is the right-hand side    SLAP........184400
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........184500
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........184600
C         that can be used to pass necessary preconditioning information SLAP........184700
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........184800
C         for the same purpose as RWORK.                                 SLAP........184900
C ITOL   :IN       Integer.                                              SLAP........185000
C         Flag to indicate the type of convergence criterion used.       SLAP........185100
C         ITOL=0  Means the  iteration stops when the test described     SLAP........185200
C                 below on  the  residual RL  is satisfied.  This is     SLAP........185300
C                 the  "Natural Stopping Criteria" for this routine.     SLAP........185400
C                 Other values  of   ITOL  cause  extra,   otherwise     SLAP........185500
C                 unnecessary, computation per iteration and     are     SLAP........185600
C                 therefore  much less  efficient.  See  ISDGMR (the     SLAP........185700
C                 stop test routine) for more information.               SLAP........185800
C         ITOL=1  Means   the  iteration stops   when the first test     SLAP........185900
C                 described below on  the residual RL  is satisfied,     SLAP........186000
C                 and there  is either right  or  no preconditioning     SLAP........186100
C                 being used.                                            SLAP........186200
C         ITOL=2  Implies     that   the  user    is   using    left     SLAP........186300
C                 preconditioning, and the second stopping criterion     SLAP........186400
C                 below is used.                                         SLAP........186500
C         ITOL=3  Means the  iteration stops   when  the  third test     SLAP........186600
C                 described below on Minv*Residual is satisfied, and     SLAP........186700
C                 there is either left  or no  preconditioning being     SLAP........186800
C                 used.                                                  SLAP........186900
C         ITOL=11 is    often  useful  for   checking  and comparing     SLAP........187000
C                 different routines.  For this case, the  user must     SLAP........187100
C                 supply  the  "exact" solution or  a  very accurate     SLAP........187200
C                 approximation (one with  an  error much less  than     SLAP........187300
C                 TOL) through a common block,                           SLAP........187400
C                     COMMON /DSLBLK/ SOLN( )                            SLAP........187500
C                 If ITOL=11, iteration stops when the 2-norm of the     SLAP........187600
C                 difference between the iterative approximation and     SLAP........187700
C                 the user-supplied solution  divided by the  2-norm     SLAP........187800
C                 of the  user-supplied solution  is  less than TOL.     SLAP........187900
C                 Note that this requires  the  user to  set up  the     SLAP........188000
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling     SLAP........188100
C                 routine.  The routine with this declaration should     SLAP........188200
C                 be loaded before the stop test so that the correct     SLAP........188300
C                 length is used by  the loader.  This procedure  is     SLAP........188400
C                 not standard Fortran and may not work correctly on     SLAP........188500
C                 your   system (although  it  has  worked  on every     SLAP........188600
C                 system the authors have tried).  If ITOL is not 11     SLAP........188700
C                 then this common block is indeed standard Fortran.     SLAP........188800
C TOL    :INOUT    Double Precision.                                     SLAP........188900
C         Convergence criterion, as described below.  If TOL is set      SLAP........189000
C         to zero on input, then a default value of 500*(the smallest    SLAP........189100
C         positive magnitude, machine epsilon) is used.                  SLAP........189200
C ITMAX  :DUMMY    Integer.                                              SLAP........189300
C         Maximum number of iterations in most SLAP routines.  In        SLAP........189400
C         this routine this does not make sense.  The maximum number     SLAP........189500
C         of iterations here is given by ITMAX = MAXL*(NRMAX+1).         SLAP........189600
C         See IGWK for definitions of MAXL and NRMAX.                    SLAP........189700
C ITER   :OUT      Integer.                                              SLAP........189800
C         Number of iterations required to reach convergence, or         SLAP........189900
C         ITMAX if convergence criterion could not be achieved in        SLAP........190000
C         ITMAX iterations.                                              SLAP........190100
C ERR    :OUT      Double Precision.                                     SLAP........190200
C         Error estimate of error in final approximate solution, as      SLAP........190300
C         defined by ITOL.  Letting norm() denote the Euclidean          SLAP........190400
C         norm, ERR is defined as follows..                              SLAP........190500
C                                                                        SLAP........190600
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........190700
C                               for right or no preconditioning, and     SLAP........190800
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........190900
C                                norm(SB*(M-inverse)*B),                 SLAP........191000
C                               for left preconditioning.                SLAP........191100
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........191200
C                               since right or no preconditioning        SLAP........191300
C                               being used.                              SLAP........191400
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........191500
C                                norm(SB*(M-inverse)*B),                 SLAP........191600
C                               since left preconditioning is being      SLAP........191700
C                               used.                                    SLAP........191800
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|        SLAP........191900
C                               i=1,n                                    SLAP........192000
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).     SLAP........192100
C IERR   :OUT      Integer.                                              SLAP........192200
C         Return error flag.                                             SLAP........192300
C               IERR = 0 => All went well.                               SLAP........192400
C               IERR = 1 => Insufficient storage allocated for           SLAP........192500
C                           RGWK or IGWK.                                SLAP........192600
C               IERR = 2 => Routine DGMRES failed to reduce the norm     SLAP........192700
C                           of the current residual on its last call,    SLAP........192800
C                           and so the iteration has stalled.  In        SLAP........192900
C                           this case, X equals the last computed        SLAP........193000
C                           approximation.  The user must either         SLAP........193100
C                           increase MAXL, or choose a different         SLAP........193200
C                           initial guess.                               SLAP........193300
C               IERR =-1 => Insufficient length for RGWK array.          SLAP........193400
C                           IGWK(6) contains the required minimum        SLAP........193500
C                           length of the RGWK array.                    SLAP........193600
C               IERR =-2 => Illegal value of ITOL, or ITOL and JPRE      SLAP........193700
C                           values are inconsistent.                     SLAP........193800
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the        SLAP........193900
C         left-hand-side of the relevant stopping test defined           SLAP........194000
C         below associated with the residual for the current             SLAP........194100
C         approximation X(L).                                            SLAP........194200
C IUNIT  :IN       Integer.                                              SLAP........194300
C         Unit number on which to write the error at each iteration,     SLAP........194400
C         if this is desired for monitoring convergence.  If unit        SLAP........194500
C         number is 0, no writing will occur.                            SLAP........194600
C SB     :IN       Double Precision SB(N).                               SLAP........194700
C         Array of length N containing scale factors for the right       SLAP........194800
C         hand side vector B.  If JSCAL.eq.0 (see below), SB need        SLAP........194900
C         not be supplied.                                               SLAP........195000
C SX     :IN       Double Precision SX(N).                               SLAP........195100
C         Array of length N containing scale factors for the solution    SLAP........195200
C         vector X.  If JSCAL.eq.0 (see below), SX need not be           SLAP........195300
C         supplied.  SB and SX can be the same array in the calling      SLAP........195400
C         program if desired.                                            SLAP........195500
C RGWK   :INOUT    Double Precision RGWK(LRGW).                          SLAP........195600
C         Double Precision array used for workspace by DGMRES.           SLAP........195700
C         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL.   SLAP........195800
C LRGW   :IN       Integer.                                              SLAP........195900
C         Length of the double precision workspace, RGWK.                SLAP........196000
C         LRGW >= 1 + N*(MAXL+6) + MAXL*(MAXL+3).                        SLAP........196100
C         See below for definition of MAXL.                              SLAP........196200
C         For the default values, RGWK has size at least 131 + 16*N.     SLAP........196300
C IGWK   :INOUT    Integer IGWK(LIGW).                                   SLAP........196400
C         The following IGWK parameters should be set by the user        SLAP........196500
C         before calling this routine.                                   SLAP........196600
C         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in       SLAP........196700
C            which X - X0 is to be found (where, X0 is the initial       SLAP........196800
C            guess).  The default value of MAXL is 10.                   SLAP........196900
C         IGWK(2) = KMP.  Maximum number of previous Krylov basis        SLAP........197000
C            vectors to which each new basis vector is made orthogonal.  SLAP........197100
C            The default value of KMP is MAXL.                           SLAP........197200
C         IGWK(3) = JSCAL.  Flag indicating whether the scaling          SLAP........197300
C            arrays SB and SX are to be used.                            SLAP........197400
C            JSCAL = 0 => SB and SX are not used and the algorithm       SLAP........197500
C               will perform as if all SB(I) = 1 and SX(I) = 1.          SLAP........197600
C            JSCAL = 1 =>  Only SX is used, and the algorithm            SLAP........197700
C               performs as if all SB(I) = 1.                            SLAP........197800
C            JSCAL = 2 =>  Only SB is used, and the algorithm            SLAP........197900
C               performs as if all SX(I) = 1.                            SLAP........198000
C            JSCAL = 3 =>  Both SB and SX are used.                      SLAP........198100
C         IGWK(4) = JPRE.  Flag indicating whether preconditioning       SLAP........198200
C            is being used.                                              SLAP........198300
C            JPRE = 0  =>  There is no preconditioning.                  SLAP........198400
C            JPRE > 0  =>  There is preconditioning on the right         SLAP........198500
C               only, and the solver will call routine MSOLVE.           SLAP........198600
C            JPRE < 0  =>  There is preconditioning on the left          SLAP........198700
C               only, and the solver will call routine MSOLVE.           SLAP........198800
C         IGWK(5) = NRMAX.  Maximum number of restarts of the            SLAP........198900
C            Krylov iteration.  The default value of NRMAX = 10.         SLAP........199000
C            if IWORK(5) = -1,  then no restarts are performed (in       SLAP........199100
C            this case, NRMAX is set to zero internally).                SLAP........199200
C         The following IWORK parameters are diagnostic information      SLAP........199300
C         made available to the user after this routine completes.       SLAP........199400
C         IGWK(6) = MLWK.  Required minimum length of RGWK array.        SLAP........199500
C         IGWK(7) = NMS.  The total number of calls to MSOLVE.           SLAP........199600
C LIGW   :IN       Integer.                                              SLAP........199700
C         Length of the integer workspace, IGWK.  LIGW >= 20.            SLAP........199800
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........199900
C         Double Precision array that can be used for workspace in       SLAP........200000
C         MSOLVE.                                                        SLAP........200100
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........200200
C         Integer array that can be used for workspace in MSOLVE.        SLAP........200300
C                                                                        SLAP........200400
C *Description:                                                          SLAP........200500
C       DGMRES solves a linear system A*X = B rewritten in the form:     SLAP........200600
C                                                                        SLAP........200700
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,                SLAP........200800
C                                                                        SLAP........200900
C       with right preconditioning, or                                   SLAP........201000
C                                                                        SLAP........201100
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,      SLAP........201200
C                                                                        SLAP........201300
C       with left preconditioning, where A is an N-by-N double precision SLAP........201400
C       matrix, X and B are N-vectors,  SB and SX  are diagonal scaling  SLAP........201500
C       matrices,   and M is  a preconditioning    matrix.   It uses     SLAP........201600
C       preconditioned  Krylov   subpace  methods  based     on  the     SLAP........201700
C       generalized minimum residual  method (GMRES).   This routine     SLAP........201800
C       optionally performs  either  the  full     orthogonalization     SLAP........201900
C       version of the  GMRES  algorithm or an incomplete variant of     SLAP........202000
C       it.  Both versions use restarting of the linear iteration by     SLAP........202100
C       default, although the user can disable this feature.             SLAP........202200
C                                                                        SLAP........202300
C       The GMRES  algorithm generates a sequence  of approximations     SLAP........202400
C       X(L) to the  true solution of the above  linear system.  The     SLAP........202500
C       convergence criteria for stopping the  iteration is based on     SLAP........202600
C       the size  of the  scaled norm of  the residual  R(L)  =  B -     SLAP........202700
C       A*X(L).  The actual stopping test is either:                     SLAP........202800
C                                                                        SLAP........202900
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),                 SLAP........203000
C                                                                        SLAP........203100
C       for right preconditioning, or                                    SLAP........203200
C                                                                        SLAP........203300
C               norm(SB*(M-inverse)*(B-A*X(L))) .le.                     SLAP........203400
C                       TOL*norm(SB*(M-inverse)*B),                      SLAP........203500
C                                                                        SLAP........203600
C       for left preconditioning, where norm() denotes the Euclidean     SLAP........203700
C       norm, and TOL is  a positive scalar less  than one  input by     SLAP........203800
C       the user.  If TOL equals zero  when DGMRES is called, then a     SLAP........203900
C       default  value  of 500*(the   smallest  positive  magnitude,     SLAP........204000
C       machine epsilon) is used.  If the  scaling arrays SB  and SX     SLAP........204100
C       are used, then  ideally they  should be chosen  so  that the     SLAP........204200
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components     SLAP........204300
C       approximately equal  to  one in  magnitude.  If one wants to     SLAP........204400
C       use the same scaling in X  and B, then  SB and SX can be the     SLAP........204500
C       same array in the calling program.                               SLAP........204600
C                                                                        SLAP........204700
C       The following is a list of the other routines and their          SLAP........204800
C       functions used by DGMRES:                                        SLAP........204900
C       DPIGMR  Contains the main iteration loop for GMRES.              SLAP........205000
C       DORTH   Orthogonalizes a new vector against older basis vectors. SLAP........205100
C       DHEQR   Computes a QR decomposition of a Hessenberg matrix.      SLAP........205200
C       DHELS   Solves a Hessenberg least-squares system, using QR       SLAP........205300
C               factors.                                                 SLAP........205400
C       DRLCAL  Computes the scaled residual RL.                         SLAP........205500
C       DXLCAL  Computes the solution XL.                                SLAP........205600
C       ISDGMR  User-replaceable stopping routine.                       SLAP........205700
C                                                                        SLAP........205800
C       This routine does  not care  what matrix data   structure is     SLAP........205900
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE     SLAP........206000
C       routines, with  the arguments as  described above.  The user     SLAP........206100
C       could write any type of structure and the appropriate MATVEC     SLAP........206200
C       and MSOLVE routines.  It is assumed  that A is stored in the     SLAP........206300
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is     SLAP........206400
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP     SLAP........206500
C       routines DSDCG and DSICCG are examples of this procedure.        SLAP........206600
C                                                                        SLAP........206700
C       Two  examples  of  matrix  data structures  are the: 1) SLAP     SLAP........206800
C       Triad  format and 2) SLAP Column format.                         SLAP........206900
C                                                                        SLAP........207000
C       =================== S L A P Triad format ===================     SLAP........207100
C       This routine requires that the  matrix A be   stored in  the     SLAP........207200
C       SLAP  Triad format.  In  this format only the non-zeros  are     SLAP........207300
C       stored.  They may appear in  *ANY* order.  The user supplies     SLAP........207400
C       three arrays of  length NELT, where  NELT is  the number  of     SLAP........207500
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For     SLAP........207600
C       each non-zero the user puts the row and column index of that     SLAP........207700
C       matrix element  in the IA and  JA arrays.  The  value of the     SLAP........207800
C       non-zero   matrix  element is  placed  in  the corresponding     SLAP........207900
C       location of the A array.   This is  an  extremely  easy data     SLAP........208000
C       structure to generate.  On  the  other hand it   is  not too     SLAP........208100
C       efficient on vector computers for  the iterative solution of     SLAP........208200
C       linear systems.  Hence,   SLAP changes   this  input    data     SLAP........208300
C       structure to the SLAP Column format  for  the iteration (but     SLAP........208400
C       does not change it back).                                        SLAP........208500
C                                                                        SLAP........208600
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........208700
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........208800
C                                                                        SLAP........208900
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........209000
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........209100
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........209200
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........209300
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........209400
C       | 0  0  0 44  0|                                                 SLAP........209500
C       |51  0 53  0 55|                                                 SLAP........209600
C                                                                        SLAP........209700
C       =================== S L A P Column format ==================     SLAP........209800
C                                                                        SLAP........209900
C       This routine  requires that  the matrix A  be stored in  the     SLAP........210000
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........210100
C       counting down columns (except for  the diagonal entry, which     SLAP........210200
C       must appear first in each  "column")  and are stored  in the     SLAP........210300
C       double precision array A.   In other words,  for each column     SLAP........210400
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........210500
C       other non-zero  elements going down  the column (except  the     SLAP........210600
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........210700
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........210800
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........210900
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........211000
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........211100
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........211200
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........211300
C       number of columns in  the matrix and NELT  is the number  of     SLAP........211400
C       non-zeros in the matrix.                                         SLAP........211500
C                                                                        SLAP........211600
C       Here is an example of the  SLAP Column  storage format for a     SLAP........211700
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........211800
C       column):                                                         SLAP........211900
C                                                                        SLAP........212000
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........212100
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........212200
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........212300
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........212400
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........212500
C       | 0  0  0 44  0|                                                 SLAP........212600
C       |51  0 53  0 55|                                                 SLAP........212700
C                                                                        SLAP........212800
C *Cautions:                                                             SLAP........212900
C     This routine will attempt to write to the Fortran logical output   SLAP........213000
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........213100
C     this logical unit is attached to a file or terminal before calling SLAP........213200
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........213300
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........213400
C                                                                        SLAP........213500
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage   SLAP........213600
C                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-  SLAP........213700
C                  more National Laboratory Report UCRL-95088, Rev. 1,   SLAP........213800
C                  Livermore, California, June 1987.                     SLAP........213900
C               2. Mark K. Seager, A SLAP for the Masses, in             SLAP........214000
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,   SLAP........214100
C                  Algorithms and Applications, Wiley, 1989, pp.135-155. SLAP........214200
C***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DPIGMR                        SLAP........214300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........214400
C   890404  DATE WRITTEN                                                 SLAP........214500
C   890404  Previous REVISION DATE                                       SLAP........214600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........214700
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........214800
C           standard.  (FNF)                                             SLAP........214900
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........215000
C   891004  Added new reference.                                         SLAP........215100
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........215200
C   910506  Corrected errors in C***ROUTINES CALLED list.  (FNF)         SLAP........215300
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........215400
C   920511  Added complete declaration section.  (WRB)                   SLAP........215500
C   920929  Corrected format of references.  (FNF)                       SLAP........215600
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)     SLAP........215700
C   921026  Added check for valid value of ITOL.  (FNF)                  SLAP........215800
C***END PROLOGUE  DGMRES                                                 SLAP........215900
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........216000
CLLL. OPTIMIZE                                                           SLAP........216100
C     .. Scalar Arguments ..                                             SLAP........216200
      DOUBLE PRECISION ERR, TOL                                          SLAP........216300
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LIGW, LRGW, N, NELT  SLAP........216400
C     .. Array Arguments ..                                              SLAP........216500
      DOUBLE PRECISION A(NELT), B(N), RGWK(LRGW), RWORK(*), SB(N),       SLAP........216600
     +                 SX(N), X(N)                                       SLAP........216700
      INTEGER IA(NELT), IGWK(LIGW), IWORK(*), JA(NELT)                   SLAP........216800
C     .. Subroutine Arguments ..                                         SLAP........216900
      EXTERNAL MATVEC, MSOLVE                                            SLAP........217000
C     .. Local Scalars ..                                                SLAP........217100
      DOUBLE PRECISION BNRM, RHOL, SUM                                   SLAP........217200
      INTEGER I, IFLAG, JPRE, JSCAL, KMP, LDL, LGMR, LHES, LQ, LR, LV,   SLAP........217300
     +        LW, LXL, LZ, LZM1, MAXL, MAXLP1, NMS, NMSL, NRMAX, NRSTS   SLAP........217400
C     .. External Functions ..                                           SLAP........217500
      DOUBLE PRECISION D1MACH, DNRM2                                     SLAP........217600
      EXTERNAL D1MACH, DNRM2                                             SLAP........217700
C     .. External Subroutines ..                                         SLAP........217800
      EXTERNAL DCOPY, DPIGMR                                             SLAP........217900
C     .. Intrinsic Functions ..                                          SLAP........218000
      INTRINSIC SQRT                                                     SLAP........218100
C***FIRST EXECUTABLE STATEMENT  DGMRES                                   SLAP........218200
      IERR = 0                                                           SLAP........218300
C   ------------------------------------------------------------------   SLAP........218400
C         Load method parameters with user values or defaults.           SLAP........218500
C   ------------------------------------------------------------------   SLAP........218600
      MAXL = IGWK(1)                                                     SLAP........218700
      IF (MAXL .EQ. 0) MAXL = 10                                         SLAP........218800
      IF (MAXL .GT. N) MAXL = N                                          SLAP........218900
      KMP = IGWK(2)                                                      SLAP........219000
      IF (KMP .EQ. 0) KMP = MAXL                                         SLAP........219100
      IF (KMP .GT. MAXL) KMP = MAXL                                      SLAP........219200
      JSCAL = IGWK(3)                                                    SLAP........219300
      JPRE = IGWK(4)                                                     SLAP........219400
C         Check for valid value of ITOL.                                 SLAP........219500
      IF( (ITOL.LT.0) .OR. ((ITOL.GT.3).AND.(ITOL.NE.11)) ) GOTO 650     SLAP........219600
C         Check for consistent values of ITOL and JPRE.                  SLAP........219700
      IF( ITOL.EQ.1 .AND. JPRE.LT.0 ) GOTO 650                           SLAP........219800
      IF( ITOL.EQ.2 .AND. JPRE.GE.0 ) GOTO 650                           SLAP........219900
      NRMAX = IGWK(5)                                                    SLAP........220000
C.....THE COMMENT LINE THAT BEGINS WITH "CCC" IMMEDIATELY BELOW          SLAP........220100
C        THIS COMMENT WAS ORIGINALLY ACTIVE CODE AND WAS COMMENTED       SLAP........220200
C        OUT DURING INTEGRATION OF SLAP WITH SUTRA.                      SLAP........220300
CCC   IF( NRMAX.EQ.0 ) NRMAX = 10                                        SLAP........220400
C         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting.   SLAP........220500
      IF( NRMAX.EQ.-1 ) NRMAX = 0                                        SLAP........220600
C         If input value of TOL is zero, set it to its default value.    SLAP........220700
      IF( TOL.EQ.0.0D0 ) TOL = 500*D1MACH(3)                             SLAP........220800
C                                                                        SLAP........220900
C         Initialize counters.                                           SLAP........221000
      ITER = 0                                                           SLAP........221100
      NMS = 0                                                            SLAP........221200
      NRSTS = 0                                                          SLAP........221300
C   ------------------------------------------------------------------   SLAP........221400
C         Form work array segment pointers.                              SLAP........221500
C   ------------------------------------------------------------------   SLAP........221600
      MAXLP1 = MAXL + 1                                                  SLAP........221700
      LV = 1                                                             SLAP........221800
      LR = LV + N*MAXLP1                                                 SLAP........221900
      LHES = LR + N + 1                                                  SLAP........222000
      LQ = LHES + MAXL*MAXLP1                                            SLAP........222100
      LDL = LQ + 2*MAXL                                                  SLAP........222200
      LW = LDL + N                                                       SLAP........222300
      LXL = LW + N                                                       SLAP........222400
      LZ = LXL + N                                                       SLAP........222500
C                                                                        SLAP........222600
C         Load IGWK(6) with required minimum length of the RGWK array.   SLAP........222700
      IGWK(6) = LZ + N - 1                                               SLAP........222800
      IF( LZ+N-1.GT.LRGW ) GOTO 640                                      SLAP........222900
C   ------------------------------------------------------------------   SLAP........223000
C         Calculate scaled-preconditioned norm of RHS vector b.          SLAP........223100
C   ------------------------------------------------------------------   SLAP........223200
      IF (JPRE .LT. 0) THEN                                              SLAP........223300
         CALL MSOLVE(N, B, RGWK(LR), NELT, IA, JA, A, ISYM,              SLAP........223400
     $        RWORK, IWORK)                                              SLAP........223500
         NMS = NMS + 1                                                   SLAP........223600
      ELSE                                                               SLAP........223700
         CALL DCOPY(N, B, 1, RGWK(LR), 1)                                SLAP........223800
      ENDIF                                                              SLAP........223900
      IF( JSCAL.EQ.2 .OR. JSCAL.EQ.3 ) THEN                              SLAP........224000
         SUM = 0                                                         SLAP........224100
         DO 10 I = 1,N                                                   SLAP........224200
            SUM = SUM + (RGWK(LR-1+I)*SB(I))**2                          SLAP........224300
 10      CONTINUE                                                        SLAP........224400
         BNRM = SQRT(SUM)                                                SLAP........224500
      ELSE                                                               SLAP........224600
         BNRM = DNRM2(N,RGWK(LR),1)                                      SLAP........224700
      ENDIF                                                              SLAP........224800
C   ------------------------------------------------------------------   SLAP........224900
C         Calculate initial residual.                                    SLAP........225000
C   ------------------------------------------------------------------   SLAP........225100
      CALL MATVEC(N, X, RGWK(LR), NELT, IA, JA, A, ISYM)                 SLAP........225200
      DO 50 I = 1,N                                                      SLAP........225300
         RGWK(LR-1+I) = B(I) - RGWK(LR-1+I)                              SLAP........225400
 50   CONTINUE                                                           SLAP........225500
C   ------------------------------------------------------------------   SLAP........225600
C         If performing restarting, then load the residual into the      SLAP........225700
C         correct location in the RGWK array.                            SLAP........225800
C   ------------------------------------------------------------------   SLAP........225900
 100  CONTINUE                                                           SLAP........226000
      IF( NRSTS.GT.NRMAX ) GOTO 610                                      SLAP........226100
      IF( NRSTS.GT.0 ) THEN                                              SLAP........226200
C         Copy the current residual to a different location in the RGWK  SLAP........226300
C         array.                                                         SLAP........226400
         CALL DCOPY(N, RGWK(LDL), 1, RGWK(LR), 1)                        SLAP........226500
      ENDIF                                                              SLAP........226600
C   ------------------------------------------------------------------   SLAP........226700
C         Use the DPIGMR algorithm to solve the linear system A*Z = R.   SLAP........226800
C   ------------------------------------------------------------------   SLAP........226900
      CALL DPIGMR(N, RGWK(LR), SB, SX, JSCAL, MAXL, MAXLP1, KMP,         SLAP........227000
     $       NRSTS, JPRE, MATVEC, MSOLVE, NMSL, RGWK(LZ), RGWK(LV),      SLAP........227100
     $       RGWK(LHES), RGWK(LQ), LGMR, RWORK, IWORK, RGWK(LW),         SLAP........227200
     $       RGWK(LDL), RHOL, NRMAX, B, BNRM, X, RGWK(LXL), ITOL,        SLAP........227300
C............THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION OF    SLAP........227400
C               SLAP WITH SUTRA.  IT ORIGINALLY READ:                    SLAP........227500
C               $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)   SLAP........227600
     $       TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR, ITMAX)       SLAP........227700
      ITER = ITER + LGMR                                                 SLAP........227800
      NMS = NMS + NMSL                                                   SLAP........227900
C                                                                        SLAP........228000
C         Increment X by the current approximate solution Z of A*Z = R.  SLAP........228100
C                                                                        SLAP........228200
      LZM1 = LZ - 1                                                      SLAP........228300
      DO 110 I = 1,N                                                     SLAP........228400
         X(I) = X(I) + RGWK(LZM1+I)                                      SLAP........228500
 110  CONTINUE                                                           SLAP........228600
      IF( IFLAG.EQ.0 ) GOTO 600                                          SLAP........228700
      IF( IFLAG.EQ.1 ) THEN                                              SLAP........228800
         NRSTS = NRSTS + 1                                               SLAP........228900
         GOTO 100                                                        SLAP........229000
      ENDIF                                                              SLAP........229100
      IF( IFLAG.EQ.2 ) GOTO 620                                          SLAP........229200
C   ------------------------------------------------------------------   SLAP........229300
C         All returns are made through this section.                     SLAP........229400
C   ------------------------------------------------------------------   SLAP........229500
C         The iteration has converged.                                   SLAP........229600
C                                                                        SLAP........229700
 600  CONTINUE                                                           SLAP........229800
      IGWK(7) = NMS                                                      SLAP........229900
      RGWK(1) = RHOL                                                     SLAP........230000
      IERR = 0                                                           SLAP........230100
      RETURN                                                             SLAP........230200
C                                                                        SLAP........230300
C         Max number((NRMAX+1)*MAXL) of linear iterations performed.     SLAP........230400
 610  CONTINUE                                                           SLAP........230500
      IGWK(7) = NMS                                                      SLAP........230600
      RGWK(1) = RHOL                                                     SLAP........230700
      IERR = 1                                                           SLAP........230800
      RETURN                                                             SLAP........230900
C                                                                        SLAP........231000
C         GMRES failed to reduce last residual in MAXL iterations.       SLAP........231100
C         The iteration has stalled.                                     SLAP........231200
 620  CONTINUE                                                           SLAP........231300
      IGWK(7) = NMS                                                      SLAP........231400
      RGWK(1) = RHOL                                                     SLAP........231500
      IERR = 2                                                           SLAP........231600
      RETURN                                                             SLAP........231700
C         Error return.  Insufficient length for RGWK array.             SLAP........231800
 640  CONTINUE                                                           SLAP........231900
      ERR = TOL                                                          SLAP........232000
      IERR = -1                                                          SLAP........232100
      RETURN                                                             SLAP........232200
C         Error return.  Inconsistent ITOL and JPRE values.              SLAP........232300
 650  CONTINUE                                                           SLAP........232400
      ERR = TOL                                                          SLAP........232500
      IERR = -2                                                          SLAP........232600
      RETURN                                                             SLAP........232700
C------------- LAST LINE OF DGMRES FOLLOWS ----------------------------  SLAP........232800
      END                                                                SLAP........232900
*DECK DHELS                                                              SLAP........233000
      SUBROUTINE DHELS (A, LDA, N, Q, B)                                 SLAP........233100
C***BEGIN PROLOGUE  DHELS                                                SLAP........233200
C***SUBSIDIARY                                                           SLAP........233300
C***PURPOSE  Internal routine for DGMRES.                                SLAP........233400
C***LIBRARY   SLATEC (SLAP)                                              SLAP........233500
C***CATEGORY  D2A4, D2B4                                                 SLAP........233600
C***TYPE      DOUBLE PRECISION (SHELS-S, DHELS-D)                        SLAP........233700
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........233800
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........233900
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........234000
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........234100
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........234200
C             Lawrence Livermore National Laboratory                     SLAP........234300
C             PO Box 808, L-60                                           SLAP........234400
C             Livermore, CA 94550 (510) 423-3141                         SLAP........234500
C***DESCRIPTION                                                          SLAP........234600
C        This routine is extracted from the LINPACK routine SGESL with   SLAP........234700
C        changes due to the fact that A is an upper Hessenberg matrix.   SLAP........234800
C                                                                        SLAP........234900
C        DHELS solves the least squares problem:                         SLAP........235000
C                                                                        SLAP........235100
C                   MIN(B-A*X,B-A*X)                                     SLAP........235200
C                                                                        SLAP........235300
C        using the factors computed by DHEQR.                            SLAP........235400
C                                                                        SLAP........235500
C *Usage:                                                                SLAP........235600
C      INTEGER LDA, N                                                    SLAP........235700
C      DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)                         SLAP........235800
C                                                                        SLAP........235900
C      CALL DHELS(A, LDA, N, Q, B)                                       SLAP........236000
C                                                                        SLAP........236100
C *Arguments:                                                            SLAP........236200
C A       :IN       Double Precision A(LDA,N)                            SLAP........236300
C          The output from DHEQR which contains the upper                SLAP........236400
C          triangular factor R in the QR decomposition of A.             SLAP........236500
C LDA     :IN       Integer                                              SLAP........236600
C          The leading dimension of the array A.                         SLAP........236700
C N       :IN       Integer                                              SLAP........236800
C          A is originally an (N+1) by N matrix.                         SLAP........236900
C Q       :IN       Double Precision Q(2*N)                              SLAP........237000
C          The coefficients of the N Givens rotations                    SLAP........237100
C          used in the QR factorization of A.                            SLAP........237200
C B       :INOUT    Double Precision B(N+1)                              SLAP........237300
C          On input, B is the right hand side vector.                    SLAP........237400
C          On output, B is the solution vector X.                        SLAP........237500
C                                                                        SLAP........237600
C***SEE ALSO  DGMRES                                                     SLAP........237700
C***ROUTINES CALLED  DAXPY                                               SLAP........237800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........237900
C   890404  DATE WRITTEN                                                 SLAP........238000
C   890404  Previous REVISION DATE                                       SLAP........238100
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........238200
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........238300
C           standard.  (FNF)                                             SLAP........238400
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........238500
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........238600
C   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)            SLAP........238700
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........238800
C   920511  Added complete declaration section.  (WRB)                   SLAP........238900
C***END PROLOGUE  DHELS                                                  SLAP........239000
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........239100
CLLL. OPTIMIZE                                                           SLAP........239200
C     .. Scalar Arguments ..                                             SLAP........239300
      INTEGER LDA, N                                                     SLAP........239400
C     .. Array Arguments ..                                              SLAP........239500
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)                              SLAP........239600
C     .. Local Scalars ..                                                SLAP........239700
      DOUBLE PRECISION C, S, T, T1, T2                                   SLAP........239800
      INTEGER IQ, K, KB, KP1                                             SLAP........239900
C     .. External Subroutines ..                                         SLAP........240000
      EXTERNAL DAXPY                                                     SLAP........240100
C***FIRST EXECUTABLE STATEMENT  DHELS                                    SLAP........240200
C                                                                        SLAP........240300
C         Minimize(B-A*X,B-A*X).  First form Q*B.                        SLAP........240400
C                                                                        SLAP........240500
      DO 20 K = 1, N                                                     SLAP........240600
         KP1 = K + 1                                                     SLAP........240700
         IQ = 2*(K-1) + 1                                                SLAP........240800
         C = Q(IQ)                                                       SLAP........240900
         S = Q(IQ+1)                                                     SLAP........241000
         T1 = B(K)                                                       SLAP........241100
         T2 = B(KP1)                                                     SLAP........241200
         B(K) = C*T1 - S*T2                                              SLAP........241300
         B(KP1) = S*T1 + C*T2                                            SLAP........241400
 20   CONTINUE                                                           SLAP........241500
C                                                                        SLAP........241600
C         Now solve  R*X = Q*B.                                          SLAP........241700
C                                                                        SLAP........241800
      DO 40 KB = 1, N                                                    SLAP........241900
         K = N + 1 - KB                                                  SLAP........242000
         B(K) = B(K)/A(K,K)                                              SLAP........242100
         T = -B(K)                                                       SLAP........242200
         CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)                          SLAP........242300
 40   CONTINUE                                                           SLAP........242400
      RETURN                                                             SLAP........242500
C------------- LAST LINE OF DHELS FOLLOWS ----------------------------   SLAP........242600
      END                                                                SLAP........242700
*DECK DHEQR                                                              SLAP........242800
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)                        SLAP........242900
C***BEGIN PROLOGUE  DHEQR                                                SLAP........243000
C***SUBSIDIARY                                                           SLAP........243100
C***PURPOSE  Internal routine for DGMRES.                                SLAP........243200
C***LIBRARY   SLATEC (SLAP)                                              SLAP........243300
C***CATEGORY  D2A4, D2B4                                                 SLAP........243400
C***TYPE      DOUBLE PRECISION (SHEQR-S, DHEQR-D)                        SLAP........243500
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........243600
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........243700
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........243800
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........243900
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........244000
C             Lawrence Livermore National Laboratory                     SLAP........244100
C             PO Box 808, L-60                                           SLAP........244200
C             Livermore, CA 94550 (510) 423-3141                         SLAP........244300
C***DESCRIPTION                                                          SLAP........244400
C        This   routine  performs  a QR   decomposition  of an  upper    SLAP........244500
C        Hessenberg matrix A using Givens  rotations.  There  are two    SLAP........244600
C        options  available: 1)  Performing  a fresh decomposition 2)    SLAP........244700
C        updating the QR factors by adding a row and  a column to the    SLAP........244800
C        matrix A.                                                       SLAP........244900
C                                                                        SLAP........245000
C *Usage:                                                                SLAP........245100
C      INTEGER LDA, N, INFO, IJOB                                        SLAP........245200
C      DOUBLE PRECISION A(LDA,N), Q(2*N)                                 SLAP........245300
C                                                                        SLAP........245400
C      CALL DHEQR(A, LDA, N, Q, INFO, IJOB)                              SLAP........245500
C                                                                        SLAP........245600
C *Arguments:                                                            SLAP........245700
C A      :INOUT    Double Precision A(LDA,N)                             SLAP........245800
C         On input, the matrix to be decomposed.                         SLAP........245900
C         On output, the upper triangular matrix R.                      SLAP........246000
C         The factorization can be written Q*A = R, where                SLAP........246100
C         Q is a product of Givens rotations and R is upper              SLAP........246200
C         triangular.                                                    SLAP........246300
C LDA    :IN       Integer                                               SLAP........246400
C         The leading dimension of the array A.                          SLAP........246500
C N      :IN       Integer                                               SLAP........246600
C         A is an (N+1) by N Hessenberg matrix.                          SLAP........246700
C Q      :OUT      Double Precision Q(2*N)                               SLAP........246800
C         The factors c and s of each Givens rotation used               SLAP........246900
C         in decomposing A.                                              SLAP........247000
C INFO   :OUT      Integer                                               SLAP........247100
C         = 0  normal value.                                             SLAP........247200
C         = K  if  A(K,K) .eq. 0.0 .  This is not an error               SLAP........247300
C           condition for this subroutine, but it does                   SLAP........247400
C           indicate that DHELS will divide by zero                      SLAP........247500
C           if called.                                                   SLAP........247600
C IJOB   :IN       Integer                                               SLAP........247700
C         = 1     means that a fresh decomposition of the                SLAP........247800
C                 matrix A is desired.                                   SLAP........247900
C         .ge. 2  means that the current decomposition of A              SLAP........248000
C                 will be updated by the addition of a row               SLAP........248100
C                 and a column.                                          SLAP........248200
C                                                                        SLAP........248300
C***SEE ALSO  DGMRES                                                     SLAP........248400
C***ROUTINES CALLED  (NONE)                                              SLAP........248500
C***REVISION HISTORY  (YYMMDD)                                           SLAP........248600
C   890404  DATE WRITTEN                                                 SLAP........248700
C   890404  Previous REVISION DATE                                       SLAP........248800
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........248900
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........249000
C           standard.  (FNF)                                             SLAP........249100
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........249200
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........249300
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........249400
C   920511  Added complete declaration section.  (WRB)                   SLAP........249500
C***END PROLOGUE  DHEQR                                                  SLAP........249600
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........249700
CLLL. OPTIMIZE                                                           SLAP........249800
C     .. Scalar Arguments ..                                             SLAP........249900
      INTEGER IJOB, INFO, LDA, N                                         SLAP........250000
C     .. Array Arguments ..                                              SLAP........250100
      DOUBLE PRECISION A(LDA,*), Q(*)                                    SLAP........250200
C     .. Local Scalars ..                                                SLAP........250300
      DOUBLE PRECISION C, S, T, T1, T2                                   SLAP........250400
      INTEGER I, IQ, J, K, KM1, KP1, NM1                                 SLAP........250500
C     .. Intrinsic Functions ..                                          SLAP........250600
      INTRINSIC ABS, SQRT                                                SLAP........250700
C***FIRST EXECUTABLE STATEMENT  DHEQR                                    SLAP........250800
      IF (IJOB .GT. 1) GO TO 70                                          SLAP........250900
C   -------------------------------------------------------------------  SLAP........251000
C         A new factorization is desired.                                SLAP........251100
C   -------------------------------------------------------------------  SLAP........251200
C         QR decomposition without pivoting.                             SLAP........251300
C                                                                        SLAP........251400
      INFO = 0                                                           SLAP........251500
      DO 60 K = 1, N                                                     SLAP........251600
         KM1 = K - 1                                                     SLAP........251700
         KP1 = K + 1                                                     SLAP........251800
C                                                                        SLAP........251900
C           Compute K-th column of R.                                    SLAP........252000
C           First, multiply the K-th column of A by the previous         SLAP........252100
C           K-1 Givens rotations.                                        SLAP........252200
C                                                                        SLAP........252300
         IF (KM1 .LT. 1) GO TO 20                                        SLAP........252400
         DO 10 J = 1, KM1                                                SLAP........252500
            I = 2*(J-1) + 1                                              SLAP........252600
            T1 = A(J,K)                                                  SLAP........252700
            T2 = A(J+1,K)                                                SLAP........252800
            C = Q(I)                                                     SLAP........252900
            S = Q(I+1)                                                   SLAP........253000
            A(J,K) = C*T1 - S*T2                                         SLAP........253100
            A(J+1,K) = S*T1 + C*T2                                       SLAP........253200
 10      CONTINUE                                                        SLAP........253300
C                                                                        SLAP........253400
C         Compute Givens components C and S.                             SLAP........253500
C                                                                        SLAP........253600
 20      CONTINUE                                                        SLAP........253700
         IQ = 2*KM1 + 1                                                  SLAP........253800
         T1 = A(K,K)                                                     SLAP........253900
         T2 = A(KP1,K)                                                   SLAP........254000
         IF( T2.EQ.0.0D0 ) THEN                                          SLAP........254100
            C = 1                                                        SLAP........254200
            S = 0                                                        SLAP........254300
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN                               SLAP........254400
            T = T1/T2                                                    SLAP........254500
            S = -1.0D0/SQRT(1.0D0+T*T)                                   SLAP........254600
            C = -S*T                                                     SLAP........254700
         ELSE                                                            SLAP........254800
            T = T2/T1                                                    SLAP........254900
            C = 1.0D0/SQRT(1.0D0+T*T)                                    SLAP........255000
            S = -C*T                                                     SLAP........255100
         ENDIF                                                           SLAP........255200
         Q(IQ) = C                                                       SLAP........255300
         Q(IQ+1) = S                                                     SLAP........255400
         A(K,K) = C*T1 - S*T2                                            SLAP........255500
         IF( A(K,K).EQ.0.0D0 ) INFO = K                                  SLAP........255600
 60   CONTINUE                                                           SLAP........255700
      RETURN                                                             SLAP........255800
C   -------------------------------------------------------------------  SLAP........255900
C         The old factorization of a will be updated.  A row and a       SLAP........256000
C         column has been added to the matrix A.  N by N-1 is now        SLAP........256100
C         the old size of the matrix.                                    SLAP........256200
C   -------------------------------------------------------------------  SLAP........256300
 70   CONTINUE                                                           SLAP........256400
      NM1 = N - 1                                                        SLAP........256500
C   -------------------------------------------------------------------  SLAP........256600
C         Multiply the new column by the N previous Givens rotations.    SLAP........256700
C   -------------------------------------------------------------------  SLAP........256800
      DO 100 K = 1,NM1                                                   SLAP........256900
         I = 2*(K-1) + 1                                                 SLAP........257000
         T1 = A(K,N)                                                     SLAP........257100
         T2 = A(K+1,N)                                                   SLAP........257200
         C = Q(I)                                                        SLAP........257300
         S = Q(I+1)                                                      SLAP........257400
         A(K,N) = C*T1 - S*T2                                            SLAP........257500
         A(K+1,N) = S*T1 + C*T2                                          SLAP........257600
 100  CONTINUE                                                           SLAP........257700
C   -------------------------------------------------------------------  SLAP........257800
C         Complete update of decomposition by forming last Givens        SLAP........257900
C         rotation, and multiplying it times the column                  SLAP........258000
C         vector(A(N,N),A(NP1,N)).                                       SLAP........258100
C   -------------------------------------------------------------------  SLAP........258200
      INFO = 0                                                           SLAP........258300
      T1 = A(N,N)                                                        SLAP........258400
      T2 = A(N+1,N)                                                      SLAP........258500
      IF ( T2.EQ.0.0D0 ) THEN                                            SLAP........258600
         C = 1                                                           SLAP........258700
         S = 0                                                           SLAP........258800
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN                                  SLAP........258900
         T = T1/T2                                                       SLAP........259000
         S = -1.0D0/SQRT(1.0D0+T*T)                                      SLAP........259100
         C = -S*T                                                        SLAP........259200
      ELSE                                                               SLAP........259300
         T = T2/T1                                                       SLAP........259400
         C = 1.0D0/SQRT(1.0D0+T*T)                                       SLAP........259500
         S = -C*T                                                        SLAP........259600
      ENDIF                                                              SLAP........259700
      IQ = 2*N - 1                                                       SLAP........259800
      Q(IQ) = C                                                          SLAP........259900
      Q(IQ+1) = S                                                        SLAP........260000
      A(N,N) = C*T1 - S*T2                                               SLAP........260100
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                    SLAP........260200
      RETURN                                                             SLAP........260300
C------------- LAST LINE OF DHEQR FOLLOWS ----------------------------   SLAP........260400
      END                                                                SLAP........260500
*DECK DLLTI2                                                             SLAP........260600
      SUBROUTINE DLLTI2 (N, B, X, NEL, IEL, JEL, EL, DINV)               SLAP........260700
C***BEGIN PROLOGUE  DLLTI2                                               SLAP........260800
C***PURPOSE  SLAP Backsolve routine for LDL' Factorization.              SLAP........260900
C            Routine to solve a system of the form  L*D*L' X = B,        SLAP........261000
C            where L is a unit lower triangular matrix and D is a        SLAP........261100
C            diagonal matrix and ' means transpose.                      SLAP........261200
C***LIBRARY   SLATEC (SLAP)                                              SLAP........261300
C***CATEGORY  D2E                                                        SLAP........261400
C***TYPE      DOUBLE PRECISION (SLLTI2-S, DLLTI2-D)                      SLAP........261500
C***KEYWORDS  INCOMPLETE FACTORIZATION, ITERATIVE PRECONDITION, SLAP,    SLAP........261600
C             SPARSE, SYMMETRIC LINEAR SYSTEM SOLVE                      SLAP........261700
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........261800
C           Seager, Mark K., (LLNL)                                      SLAP........261900
C             Lawrence Livermore National Laboratory                     SLAP........262000
C             PO BOX 808, L-60                                           SLAP........262100
C             Livermore, CA 94550 (510) 423-3141                         SLAP........262200
C             seager@llnl.gov                                            SLAP........262300
C***DESCRIPTION                                                          SLAP........262400
C                                                                        SLAP........262500
C *Usage:                                                                SLAP........262600
C     INTEGER N, NEL, IEL(NEL), JEL(NEL)                                 SLAP........262700
C     DOUBLE PRECISION B(N), X(N), EL(NEL), DINV(N)                      SLAP........262800
C                                                                        SLAP........262900
C     CALL DLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )                    SLAP........263000
C                                                                        SLAP........263100
C *Arguments:                                                            SLAP........263200
C N      :IN       Integer                                               SLAP........263300
C         Order of the Matrix.                                           SLAP........263400
C B      :IN       Double Precision B(N).                                SLAP........263500
C         Right hand side vector.                                        SLAP........263600
C X      :OUT      Double Precision X(N).                                SLAP........263700
C         Solution to L*D*L' x = b.                                      SLAP........263800
C NEL    :IN       Integer.                                              SLAP........263900
C         Number of non-zeros in the EL array.                           SLAP........264000
C IEL    :IN       Integer IEL(NEL).                                     SLAP........264100
C JEL    :IN       Integer JEL(NEL).                                     SLAP........264200
C EL     :IN       Double Precision     EL(NEL).                         SLAP........264300
C         IEL, JEL, EL contain the unit lower triangular factor   of     SLAP........264400
C         the incomplete decomposition   of the A  matrix  stored in     SLAP........264500
C         SLAP Row format.   The diagonal of ones *IS* stored.  This     SLAP........264600
C         structure can be set  up  by  the DS2LT routine.  See  the     SLAP........264700
C         "Description", below for more details about the  SLAP  Row     SLAP........264800
C         format.                                                        SLAP........264900
C DINV   :IN       Double Precision DINV(N).                             SLAP........265000
C         Inverse of the diagonal matrix D.                              SLAP........265100
C                                                                        SLAP........265200
C *Description:                                                          SLAP........265300
C       This routine is supplied with  the SLAP package as a routine     SLAP........265400
C       to perform the MSOLVE operation in the SCG iteration routine     SLAP........265500
C       for  the driver  routine DSICCG.   It must be called via the     SLAP........265600
C       SLAP  MSOLVE calling sequence  convention  interface routine     SLAP........265700
C       DSLLI.                                                           SLAP........265800
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****          SLAP........265900
C               **** SLAP MSOLVE CALLING CONVENTION ****                 SLAP........266000
C                                                                        SLAP........266100
C       IEL, JEL, EL should contain the unit lower triangular factor     SLAP........266200
C       of  the incomplete decomposition of  the A matrix  stored in     SLAP........266300
C       SLAP Row format.   This IC factorization  can be computed by     SLAP........266400
C       the  DSICS routine.  The  diagonal  (which is all one's) is      SLAP........266500
C       stored.                                                          SLAP........266600
C                                                                        SLAP........266700
C       ==================== S L A P Row format ====================     SLAP........266800
C                                                                        SLAP........266900
C       This routine requires  that the matrix A  be  stored  in the     SLAP........267000
C       SLAP  Row format.   In this format  the non-zeros are stored     SLAP........267100
C       counting across  rows (except for the diagonal  entry, which     SLAP........267200
C       must  appear first  in each  "row")  and  are stored  in the     SLAP........267300
C       double precision  array A.  In other words, for each row  in     SLAP........267400
C       the matrix  put the diagonal  entry in A.   Then put in  the     SLAP........267500
C       other  non-zero elements  going across  the row  (except the     SLAP........267600
C       diagonal) in order.  The JA array holds the column index for     SLAP........267700
C       each non-zero.  The IA array holds the offsets  into the JA,     SLAP........267800
C       A  arrays  for  the   beginning  of  each  row.    That  is,     SLAP........267900
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-     SLAP........268000
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)     SLAP........268100
C       are  the last elements  of the  IROW-th row.   Note  that we     SLAP........268200
C       always have  IA(N+1) = NELT+1, where N is the number of rows     SLAP........268300
C       in the matrix  and  NELT is the  number of non-zeros  in the     SLAP........268400
C       matrix.                                                          SLAP........268500
C                                                                        SLAP........268600
C       Here is an example of the SLAP Row storage format for a  5x5     SLAP........268700
C       Matrix (in the A and JA arrays '|' denotes the end of a row):    SLAP........268800
C                                                                        SLAP........268900
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.   SLAP........269000
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........269100
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53   SLAP........269200
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........269300
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12                      SLAP........269400
C       | 0  0  0 44  0|                                                 SLAP........269500
C       |51  0 53  0 55|                                                 SLAP........269600
C                                                                        SLAP........269700
C       With  the SLAP  Row format  the "inner loop" of this routine     SLAP........269800
C       should vectorize   on machines with   hardware  support  for     SLAP........269900
C       vector gather/scatter operations.  Your compiler may require     SLAP........270000
C       a  compiler directive  to  convince   it that there  are  no     SLAP........270100
C       implicit vector  dependencies.  Compiler directives  for the     SLAP........270200
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied     SLAP........270300
C       with the standard SLAP distribution.                             SLAP........270400
C                                                                        SLAP........270500
C***SEE ALSO  DSICCG, DSICS                                              SLAP........270600
C***REFERENCES  (NONE)                                                   SLAP........270700
C***ROUTINES CALLED  (NONE)                                              SLAP........270800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........270900
C   871119  DATE WRITTEN                                                 SLAP........271000
C   881213  Previous REVISION DATE                                       SLAP........271100
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........271200
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........271300
C           standard.  (FNF)                                             SLAP........271400
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........271500
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........271600
C   920511  Added complete declaration section.  (WRB)                   SLAP........271700
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........271800
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........271900
C***END PROLOGUE  DLLTI2                                                 SLAP........272000
C     .. Scalar Arguments ..                                             SLAP........272100
      INTEGER N, NEL                                                     SLAP........272200
C     .. Array Arguments ..                                              SLAP........272300
      DOUBLE PRECISION B(N), DINV(N), EL(NEL), X(N)                      SLAP........272400
      INTEGER IEL(NEL), JEL(NEL)                                         SLAP........272500
C     .. Local Scalars ..                                                SLAP........272600
      INTEGER I, IBGN, IEND, IROW                                        SLAP........272700
C***FIRST EXECUTABLE STATEMENT  DLLTI2                                   SLAP........272800
C                                                                        SLAP........272900
C         Solve  L*y = b,  storing result in x.                          SLAP........273000
C                                                                        SLAP........273100
      DO 10 I=1,N                                                        SLAP........273200
         X(I) = B(I)                                                     SLAP........273300
 10   CONTINUE                                                           SLAP........273400
      DO 30 IROW = 1, N                                                  SLAP........273500
         IBGN = IEL(IROW) + 1                                            SLAP........273600
         IEND = IEL(IROW+1) - 1                                          SLAP........273700
         IF( IBGN.LE.IEND ) THEN                                         SLAP........273800
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........273900
CDIR$ IVDEP                                                              SLAP........274000
CVD$ NOCONCUR                                                            SLAP........274100
CVD$ NODEPCHK                                                            SLAP........274200
            DO 20 I = IBGN, IEND                                         SLAP........274300
               X(IROW) = X(IROW) - EL(I)*X(JEL(I))                       SLAP........274400
 20         CONTINUE                                                     SLAP........274500
         ENDIF                                                           SLAP........274600
 30   CONTINUE                                                           SLAP........274700
C                                                                        SLAP........274800
C         Solve  D*Z = Y,  storing result in X.                          SLAP........274900
C                                                                        SLAP........275000
      DO 40 I=1,N                                                        SLAP........275100
         X(I) = X(I)*DINV(I)                                             SLAP........275200
 40   CONTINUE                                                           SLAP........275300
C                                                                        SLAP........275400
C         Solve  L-trans*X = Z.                                          SLAP........275500
C                                                                        SLAP........275600
      DO 60 IROW = N, 2, -1                                              SLAP........275700
         IBGN = IEL(IROW) + 1                                            SLAP........275800
         IEND = IEL(IROW+1) - 1                                          SLAP........275900
         IF( IBGN.LE.IEND ) THEN                                         SLAP........276000
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........276100
CDIR$ IVDEP                                                              SLAP........276200
CVD$ NOCONCUR                                                            SLAP........276300
CVD$ NODEPCHK                                                            SLAP........276400
            DO 50 I = IBGN, IEND                                         SLAP........276500
               X(JEL(I)) = X(JEL(I)) - EL(I)*X(IROW)                     SLAP........276600
 50         CONTINUE                                                     SLAP........276700
         ENDIF                                                           SLAP........276800
 60   CONTINUE                                                           SLAP........276900
C                                                                        SLAP........277000
      RETURN                                                             SLAP........277100
C------------- LAST LINE OF DLLTI2 FOLLOWS ----------------------------  SLAP........277200
      END                                                                SLAP........277300
*DECK DNRM2                                                              SLAP........277400
      DOUBLE PRECISION FUNCTION DNRM2 (N, DX, INCX)                      SLAP........277500
C***BEGIN PROLOGUE  DNRM2                                                SLAP........277600
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.         SLAP........277700
C***LIBRARY   SLATEC (BLAS)                                              SLAP........277800
C***CATEGORY  D1A3B                                                      SLAP........277900
C***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)              SLAP........278000
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,                SLAP........278100
C             LINEAR ALGEBRA, UNITARY, VECTOR                            SLAP........278200
C***AUTHOR  Lawson, C. L., (JPL)                                         SLAP........278300
C           Hanson, R. J., (SNLA)                                        SLAP........278400
C           Kincaid, D. R., (U. of Texas)                                SLAP........278500
C           Krogh, F. T., (JPL)                                          SLAP........278600
C***DESCRIPTION                                                          SLAP........278700
C                                                                        SLAP........278800
C                B L A S  Subprogram                                     SLAP........278900
C    Description of parameters                                           SLAP........279000
C                                                                        SLAP........279100
C     --Input--                                                          SLAP........279200
C        N  number of elements in input vector(s)                        SLAP........279300
C       DX  double precision vector with N elements                      SLAP........279400
C     INCX  storage spacing between elements of DX                       SLAP........279500
C                                                                        SLAP........279600
C     --Output--                                                         SLAP........279700
C    DNRM2  double precision result (zero if N .LE. 0)                   SLAP........279800
C                                                                        SLAP........279900
C     Euclidean norm of the N-vector stored in DX with storage           SLAP........280000
C     increment INCX.                                                    SLAP........280100
C     If N .LE. 0, return with result = 0.                               SLAP........280200
C     If N .GE. 1, then INCX must be .GE. 1                              SLAP........280300
C                                                                        SLAP........280400
C     Four phase method using two built-in constants that are            SLAP........280500
C     hopefully applicable to all machines.                              SLAP........280600
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.      SLAP........280700
C         CUTHI = minimum of  SQRT(V)      over all known machines.      SLAP........280800
C     where                                                              SLAP........280900
C         EPS = smallest no. such that EPS + 1. .GT. 1.                  SLAP........281000
C         U   = smallest positive no.   (underflow limit)                SLAP........281100
C         V   = largest  no.            (overflow  limit)                SLAP........281200
C                                                                        SLAP........281300
C     Brief outline of algorithm.                                        SLAP........281400
C                                                                        SLAP........281500
C     Phase 1 scans zero components.                                     SLAP........281600
C     move to phase 2 when a component is nonzero and .LE. CUTLO         SLAP........281700
C     move to phase 3 when a component is .GT. CUTLO                     SLAP........281800
C     move to phase 4 when a component is .GE. CUTHI/M                   SLAP........281900
C     where M = N for X() real and M = 2*N for complex.                  SLAP........282000
C                                                                        SLAP........282100
C     Values for CUTLO and CUTHI.                                        SLAP........282200
C     From the environmental parameters listed in the IMSL converter     SLAP........282300
C     document the limiting values are as follows:                       SLAP........282400
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are SLAP........282500
C                   Univac and DEC at 2**(-103)                          SLAP........282600
C                   Thus CUTLO = 2**(-51) = 4.44089E-16                  SLAP........282700
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.           SLAP........282800
C                   Thus CUTHI = 2**(63.5) = 1.30438E19                  SLAP........282900
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.              SLAP........283000
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11                SLAP........283100
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19                     SLAP........283200
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/                           SLAP........283300
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/                           SLAP........283400
C                                                                        SLAP........283500
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.      SLAP........283600
C                 Krogh, Basic linear algebra subprograms for Fortran    SLAP........283700
C                 usage, Algorithm No. 539, Transactions on Mathematical SLAP........283800
C                 Software 5, 3 (September 1979), pp. 308-323.           SLAP........283900
C***ROUTINES CALLED  (NONE)                                              SLAP........284000
C***REVISION HISTORY  (YYMMDD)                                           SLAP........284100
C   791001  DATE WRITTEN                                                 SLAP........284200
C   890531  Changed all specific intrinsics to generic.  (WRB)           SLAP........284300
C   890831  Modified array declarations.  (WRB)                          SLAP........284400
C   890831  REVISION DATE from Version 3.2                               SLAP........284500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........284600
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........284700
C***END PROLOGUE  DNRM2                                                  SLAP........284800
      INTEGER NEXT                                                       SLAP........284900
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,     SLAP........285000
     +                 ONE                                               SLAP........285100
      SAVE CUTLO, CUTHI, ZERO, ONE                                       SLAP........285200
      DATA ZERO, ONE /0.0D0, 1.0D0/                                      SLAP........285300
C                                                                        SLAP........285400
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/                           SLAP........285500
C***FIRST EXECUTABLE STATEMENT  DNRM2                                    SLAP........285600
      IF (N .GT. 0) GO TO 10                                             SLAP........285700
         DNRM2  = ZERO                                                   SLAP........285800
         GO TO 300                                                       SLAP........285900
C                                                                        SLAP........286000
   10 ASSIGN 30 TO NEXT                                                  SLAP........286100
      SUM = ZERO                                                         SLAP........286200
      NN = N * INCX                                                      SLAP........286300
C                                                                        SLAP........286400
C                                                 BEGIN MAIN LOOP        SLAP........286500
C                                                                        SLAP........286600
      I = 1                                                              SLAP........286700
   20    GO TO NEXT,(30, 50, 70, 110)                                    SLAP........286800
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85                                SLAP........286900
      ASSIGN 50 TO NEXT                                                  SLAP........287000
      XMAX = ZERO                                                        SLAP........287100
C                                                                        SLAP........287200
C                        PHASE 1.  SUM IS ZERO                           SLAP........287300
C                                                                        SLAP........287400
   50 IF (DX(I) .EQ. ZERO) GO TO 200                                     SLAP........287500
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85                                SLAP........287600
C                                                                        SLAP........287700
C                                PREPARE FOR PHASE 2.                    SLAP........287800
C                                                                        SLAP........287900
      ASSIGN 70 TO NEXT                                                  SLAP........288000
      GO TO 105                                                          SLAP........288100
C                                                                        SLAP........288200
C                                PREPARE FOR PHASE 4.                    SLAP........288300
C                                                                        SLAP........288400
  100 I = J                                                              SLAP........288500
      ASSIGN 110 TO NEXT                                                 SLAP........288600
      SUM = (SUM / DX(I)) / DX(I)                                        SLAP........288700
  105 XMAX = ABS(DX(I))                                                  SLAP........288800
      GO TO 115                                                          SLAP........288900
C                                                                        SLAP........289000
C                   PHASE 2.  SUM IS SMALL.                              SLAP........289100
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.      SLAP........289200
C                                                                        SLAP........289300
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75                                SLAP........289400
C                                                                        SLAP........289500
C                     COMMON CODE FOR PHASES 2 AND 4.                    SLAP........289600
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. SLAP........289700
C                                                                        SLAP........289800
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115                                SLAP........289900
         SUM = ONE + SUM * (XMAX / DX(I))**2                             SLAP........290000
         XMAX = ABS(DX(I))                                               SLAP........290100
         GO TO 200                                                       SLAP........290200
C                                                                        SLAP........290300
  115 SUM = SUM + (DX(I)/XMAX)**2                                        SLAP........290400
      GO TO 200                                                          SLAP........290500
C                                                                        SLAP........290600
C                  PREPARE FOR PHASE 3.                                  SLAP........290700
C                                                                        SLAP........290800
   75 SUM = (SUM * XMAX) * XMAX                                          SLAP........290900
C                                                                        SLAP........291000
C     FOR REAL OR D.P. SET HITEST = CUTHI/N                              SLAP........291100
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)                          SLAP........291200
C                                                                        SLAP........291300
   85 HITEST = CUTHI / N                                                 SLAP........291400
C                                                                        SLAP........291500
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.             SLAP........291600
C                                                                        SLAP........291700
      DO 95 J = I,NN,INCX                                                SLAP........291800
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100                              SLAP........291900
   95    SUM = SUM + DX(J)**2                                            SLAP........292000
      DNRM2 = SQRT(SUM)                                                  SLAP........292100
      GO TO 300                                                          SLAP........292200
C                                                                        SLAP........292300
  200 CONTINUE                                                           SLAP........292400
      I = I + INCX                                                       SLAP........292500
      IF (I .LE. NN) GO TO 20                                            SLAP........292600
C                                                                        SLAP........292700
C              END OF MAIN LOOP.                                         SLAP........292800
C                                                                        SLAP........292900
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.               SLAP........293000
C                                                                        SLAP........293100
      DNRM2 = XMAX * SQRT(SUM)                                           SLAP........293200
  300 CONTINUE                                                           SLAP........293300
      RETURN                                                             SLAP........293400
      END                                                                SLAP........293500
*DECK DOMN                                                               SLAP........293600
      SUBROUTINE DOMN (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,   SLAP........293700
     +   NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,   SLAP........293800
     +   EMAP, DZ, CSAV, RWORK, IWORK)                                   SLAP........293900
C***BEGIN PROLOGUE  DOMN                                                 SLAP........294000
C***PURPOSE  Preconditioned Orthomin Sparse Iterative Ax=b Solver.       SLAP........294100
C            Routine to solve a general linear system  Ax = b  using     SLAP........294200
C            the Preconditioned Orthomin method.                         SLAP........294300
C***LIBRARY   SLATEC (SLAP)                                              SLAP........294400
C***CATEGORY  D2A4, D2B4                                                 SLAP........294500
C***TYPE      DOUBLE PRECISION (SOMN-S, DOMN-D)                          SLAP........294600
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM,       SLAP........294700
C             ORTHOMIN, SLAP, SPARSE                                     SLAP........294800
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........294900
C           Seager, Mark K., (LLNL)                                      SLAP........295000
C             Lawrence Livermore National Laboratory                     SLAP........295100
C             PO BOX 808, L-60                                           SLAP........295200
C             Livermore, CA 94550 (510) 423-3141                         SLAP........295300
C             seager@llnl.gov                                            SLAP........295400
C***DESCRIPTION                                                          SLAP........295500
C                                                                        SLAP........295600
C *Usage:                                                                SLAP........295700
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX     SLAP........295800
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)                    SLAP........295900
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)         SLAP........296000
C     DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE)      SLAP........296100
C     DOUBLE PRECISION DZ(N), CSAV(NSAVE), RWORK(USER DEFINED)           SLAP........296200
C     EXTERNAL MATVEC, MSOLVE                                            SLAP........296300
C                                                                        SLAP........296400
C     CALL DOMN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,          SLAP........296500
C    $     NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R,           SLAP........296600
C    $     Z, P, AP, EMAP, DZ, CSAV, RWORK, IWORK)                       SLAP........296700
C                                                                        SLAP........296800
C *Arguments:                                                            SLAP........296900
C N      :IN       Integer.                                              SLAP........297000
C         Order of the Matrix.                                           SLAP........297100
C B      :IN       Double Precision B(N).                                SLAP........297200
C         Right-hand side vector.                                        SLAP........297300
C X      :INOUT    Double Precision X(N).                                SLAP........297400
C         On input X is your initial guess for solution vector.          SLAP........297500
C         On output X is the final approximate solution.                 SLAP........297600
C NELT   :IN       Integer.                                              SLAP........297700
C         Number of Non-Zeros stored in A.                               SLAP........297800
C IA     :IN       Integer IA(NELT).                                     SLAP........297900
C JA     :IN       Integer JA(NELT).                                     SLAP........298000
C A      :IN       Double Precision A(NELT).                             SLAP........298100
C         These arrays contain the matrix data structure for A.          SLAP........298200
C         It could take any form.  See "Description", below, for more    SLAP........298300
C         details.                                                       SLAP........298400
C ISYM   :IN       Integer.                                              SLAP........298500
C         Flag to indicate symmetric storage format.                     SLAP........298600
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........298700
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........298800
C         or lower triangle of the matrix is stored.                     SLAP........298900
C MATVEC :EXT      External.                                             SLAP........299000
C         Name of a routine which performs the matrix vector multiply    SLAP........299100
C         Y = A*X given A and X.  The name of the MATVEC routine must    SLAP........299200
C         be declared external in the calling program.  The calling      SLAP........299300
C         sequence to MATVEC is:                                         SLAP........299400
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )              SLAP........299500
C         Where N is the number of unknowns, Y is the product A*X        SLAP........299600
C         upon return X is an input vector, NELT is the number of        SLAP........299700
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.      SLAP........299800
C         ISYM is a flag which, if non-zero, denotest that A is          SLAP........299900
C         symmetric and only the lower or upper triangle is stored.      SLAP........300000
C MSOLVE :EXT      External.                                             SLAP........300100
C         Name of a routine which solves a linear system MZ = R for      SLAP........300200
C         Z given R with the preconditioning matrix M (M is supplied via SLAP........300300
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  SLAP........300400
C         be declared external in the calling program.  The calling      SLAP........300500
C         sequence to MSOLVE is:                                         SLAP........300600
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........300700
C         Where N is the number of unknowns, R is the right-hand side    SLAP........300800
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........300900
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........301000
C         that can be used to pass necessary preconditioning information SLAP........301100
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........301200
C         for the same purpose as RWORK.                                 SLAP........301300
C NSAVE  :IN       Integer.                                              SLAP........301400
C         Number of  direction vectors to save and orthogonalize         SLAP........301500
C         against.  NSAVE >= 0.                                          SLAP........301600
C ITOL   :IN       Integer.                                              SLAP........301700
C         Flag to indicate type of convergence criterion.                SLAP........301800
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........301900
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........302000
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........302100
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........302200
C         side is less than TOL, where M-inv is the inverse of the       SLAP........302300
C         diagonal of A.                                                 SLAP........302400
C         ITOL=11 is often useful for checking and comparing different   SLAP........302500
C         routines.  For this case, the user must supply the "exact"     SLAP........302600
C         solution or a very accurate approximation (one with an error   SLAP........302700
C         much less than TOL) through a common block,                    SLAP........302800
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........302900
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........303000
C         between the iterative approximation and the user-supplied      SLAP........303100
C         solution divided by the 2-norm of the user-supplied solution   SLAP........303200
C         is less than TOL.  Note that this requires the user to set up  SLAP........303300
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........303400
C         The routine with this declaration should be loaded before the  SLAP........303500
C         stop test so that the correct length is used by the loader.    SLAP........303600
C         This procedure is not standard Fortran and may not work        SLAP........303700
C         correctly on your system (although it has worked on every      SLAP........303800
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........303900
C         common block is indeed standard Fortran.                       SLAP........304000
C TOL    :INOUT    Double Precision.                                     SLAP........304100
C         Convergence criterion, as described above.  (Reset if IERR=4.) SLAP........304200
C ITMAX  :IN       Integer.                                              SLAP........304300
C         Maximum number of iterations.                                  SLAP........304400
CC..THE NEXT FOUR COMMENT LINES *NOT* MARKED WITH "CC" WERE MODIFIED     SLAP........304500
CC    DURING INTEGRATION OF SLAP WITH SUTRA.  THEY ORIGINALLY READ:      SLAP........304600
CC  C ITER   :OUT      Integer.                                          SLAP........304700
CC  C         Number of iterations required to reach convergence, or     SLAP........304800
CC  C         ITMAX+1 if convergence criterion could not be achieved in  SLAP........304900
CC  C         ITMAX iterations.                                          SLAP........305000
C ITER   :OUT      Integer.                                              SLAP........305100
C         Number of iterations required to reach convergence, or         SLAP........305200
C         ITMAX if convergence criterion could not be achieved in        SLAP........305300
C         ITMAX iterations.                                              SLAP........305400
C ERR    :OUT      Double Precision.                                     SLAP........305500
C         Error estimate of error in final approximate solution, as      SLAP........305600
C         defined by ITOL.                                               SLAP........305700
C IERR   :OUT      Integer.                                              SLAP........305800
C         Return error flag.                                             SLAP........305900
C           IERR = 0 => All went well.                                   SLAP........306000
C           IERR = 1 => Insufficient space allocated for WORK or IWORK.  SLAP........306100
C           IERR = 2 => Method failed to converge in ITMAX steps.        SLAP........306200
C           IERR = 3 => Error in user input.                             SLAP........306300
C                       Check input values of N, ITOL.                   SLAP........306400
C           IERR = 4 => User error tolerance set too tight.              SLAP........306500
C                       Reset to 500*D1MACH(3).  Iteration proceeded.    SLAP........306600
C           IERR = 5 => Preconditioning matrix, M, is not positive       SLAP........306700
C                       definite.  (r,z) < 0.                            SLAP........306800
C           IERR = 6 => Breakdown of method detected.                    SLAP........306900
C                       (p,Ap) < epsilon**2.                             SLAP........307000
C IUNIT  :IN       Integer.                                              SLAP........307100
C         Unit number on which to write the error at each iteration,     SLAP........307200
C         if this is desired for monitoring convergence.  If unit        SLAP........307300
C         number is 0, no writing will occur.                            SLAP........307400
C R      :WORK     Double Precision R(N).                                SLAP........307500
C Z      :WORK     Double Precision Z(N).                                SLAP........307600
C P      :WORK     Double Precision P(N,0:NSAVE).                        SLAP........307700
C AP     :WORK     Double Precision AP(N,0:NSAVE).                       SLAP........307800
C EMAP   :WORK     Double Precision EMAP(N,0:NSAVE).                     SLAP........307900
C DZ     :WORK     Double Precision DZ(N).                               SLAP........308000
C CSAV   :WORK     Double Precision CSAV(NSAVE)                          SLAP........308100
C         Double Precision arrays used for workspace.                    SLAP........308200
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........308300
C         Double Precision array that can be used for workspace in       SLAP........308400
C         MSOLVE.                                                        SLAP........308500
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........308600
C         Integer array that can be used for workspace in MSOLVE.        SLAP........308700
C                                                                        SLAP........308800
C *Description                                                           SLAP........308900
C       This routine does  not care  what matrix data   structure is     SLAP........309000
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE     SLAP........309100
C       routines, with  the arguments as  described above.  The user     SLAP........309200
C       could write any type of structure and the appropriate MATVEC     SLAP........309300
C       and MSOLVE routines.  It is assumed  that A is stored in the     SLAP........309400
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is     SLAP........309500
C       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP     SLAP........309600
C       routines DSDOMN and DSLUOM are examples of this procedure.       SLAP........309700
C                                                                        SLAP........309800
C       Two  examples  of  matrix  data structures  are the: 1) SLAP     SLAP........309900
C       Triad  format and 2) SLAP Column format.                         SLAP........310000
C                                                                        SLAP........310100
C       =================== S L A P Triad format ===================     SLAP........310200
C       In  this   format only the  non-zeros are  stored.  They may     SLAP........310300
C       appear  in *ANY* order.   The user  supplies three arrays of     SLAP........310400
C       length NELT, where  NELT  is the number  of non-zeros in the     SLAP........310500
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero     SLAP........310600
C       the  user puts   the row  and  column index   of that matrix     SLAP........310700
C       element in the IA and JA arrays.  The  value of the non-zero     SLAP........310800
C       matrix  element is  placed in  the corresponding location of     SLAP........310900
C       the A  array.  This is  an extremely easy data  structure to     SLAP........311000
C       generate.  On  the other hand it  is  not too  efficient  on     SLAP........311100
C       vector  computers   for the  iterative  solution  of  linear     SLAP........311200
C       systems.  Hence, SLAP  changes this input  data structure to     SLAP........311300
C       the SLAP   Column  format for the  iteration (but   does not     SLAP........311400
C       change it back).                                                 SLAP........311500
C                                                                        SLAP........311600
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........311700
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........311800
C                                                                        SLAP........311900
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........312000
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........312100
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........312200
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........312300
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........312400
C       | 0  0  0 44  0|                                                 SLAP........312500
C       |51  0 53  0 55|                                                 SLAP........312600
C                                                                        SLAP........312700
C       =================== S L A P Column format ==================     SLAP........312800
C                                                                        SLAP........312900
C       In  this format   the non-zeros are    stored counting  down     SLAP........313000
C       columns (except  for the diagonal  entry, which must  appear     SLAP........313100
C       first  in each "column") and are  stored in the  double pre-     SLAP........313200
C       cision array  A. In  other  words,  for each  column  in the     SLAP........313300
C       matrix  first put  the diagonal entry in A.  Then put in the     SLAP........313400
C       other non-zero  elements going  down the column  (except the     SLAP........313500
C       diagonal)  in order.  The IA array  holds the  row index for     SLAP........313600
C       each non-zero.  The JA array  holds the offsets into the IA,     SLAP........313700
C       A  arrays  for  the  beginning  of  each  column.  That  is,     SLAP........313800
C       IA(JA(ICOL)),A(JA(ICOL)) are the first elements of the ICOL-     SLAP........313900
C       th column in IA and A, and IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)     SLAP........314000
C       are  the last elements of the ICOL-th column.   Note that we     SLAP........314100
C       always have JA(N+1)=NELT+1, where N is the number of columns     SLAP........314200
C       in the matrix  and NELT  is the number  of non-zeros  in the     SLAP........314300
C       matrix.                                                          SLAP........314400
C                                                                        SLAP........314500
C       Here is an example of the  SLAP Column  storage format for a     SLAP........314600
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........314700
C       column):                                                         SLAP........314800
C                                                                        SLAP........314900
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........315000
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........315100
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........315200
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........315300
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........315400
C       | 0  0  0 44  0|                                                 SLAP........315500
C       |51  0 53  0 55|                                                 SLAP........315600
C                                                                        SLAP........315700
C *Cautions:                                                             SLAP........315800
C     This routine will attempt to write to the Fortran logical output   SLAP........315900
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........316000
C     this logical unit is attached to a file or terminal before calling SLAP........316100
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........316200
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........316300
C                                                                        SLAP........316400
C***SEE ALSO  DSDOMN, DSLUOM, ISDOMN                                     SLAP........316500
C***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in             SLAP........316600
C                  G. F. Carey, Ed., Parallel Supercomputing: Methods,   SLAP........316700
C                  Algorithms and Applications, Wiley, 1989, pp.135-155. SLAP........316800
C***ROUTINES CALLED  D1MACH, DAXPY, DCOPY, DDOT, ISDOMN                  SLAP........316900
C***REVISION HISTORY  (YYMMDD)                                           SLAP........317000
C   890404  DATE WRITTEN                                                 SLAP........317100
C   890404  Previous REVISION DATE                                       SLAP........317200
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........317300
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........317400
C           standard.  (FNF)                                             SLAP........317500
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........317600
C   891004  Added new reference.                                         SLAP........317700
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........317800
C   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)  SLAP........317900
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........318000
C   920511  Added complete declaration section.  (WRB)                   SLAP........318100
C   920929  Corrected format of reference.  (FNF)                        SLAP........318200
C   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)     SLAP........318300
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........318400
C   930326  Removed unused variable.  (FNF)                              SLAP........318500
C***END PROLOGUE  DOMN                                                   SLAP........318600
C     .. Scalar Arguments ..                                             SLAP........318700
      DOUBLE PRECISION ERR, TOL                                          SLAP........318800
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT, NSAVE       SLAP........318900
C     .. Array Arguments ..                                              SLAP........319000
      DOUBLE PRECISION A(NELT), AP(N,0:NSAVE), B(N), CSAV(NSAVE),        SLAP........319100
     +                 DZ(N), EMAP(N,0:NSAVE), P(N,0:NSAVE), R(N),       SLAP........319200
     +                 RWORK(*), X(N), Z(N)                              SLAP........319300
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........319400
C     .. Subroutine Arguments ..                                         SLAP........319500
      EXTERNAL MATVEC, MSOLVE                                            SLAP........319600
C     .. Local Scalars ..                                                SLAP........319700
      DOUBLE PRECISION AK, AKDEN, AKNUM, BKL, BNRM, FUZZ, SOLNRM         SLAP........319800
      INTEGER I, IP, IPO, K, L, LMAX                                     SLAP........319900
C     .. External Functions ..                                           SLAP........320000
      DOUBLE PRECISION D1MACH, DDOT                                      SLAP........320100
      INTEGER ISDOMN                                                     SLAP........320200
      EXTERNAL D1MACH, DDOT, ISDOMN                                      SLAP........320300
C     .. External Subroutines ..                                         SLAP........320400
      EXTERNAL DAXPY, DCOPY                                              SLAP........320500
C     .. Intrinsic Functions ..                                          SLAP........320600
      INTRINSIC ABS, MIN, MOD                                            SLAP........320700
C***FIRST EXECUTABLE STATEMENT  DOMN                                     SLAP........320800
C                                                                        SLAP........320900
C         Check some of the input data.                                  SLAP........321000
C                                                                        SLAP........321100
      ITER = 0                                                           SLAP........321200
      IERR = 0                                                           SLAP........321300
      IF( N.LT.1 ) THEN                                                  SLAP........321400
         IERR = 3                                                        SLAP........321500
         RETURN                                                          SLAP........321600
      ENDIF                                                              SLAP........321700
      FUZZ = D1MACH(3)                                                   SLAP........321800
      IF( TOL.LT.500*FUZZ ) THEN                                         SLAP........321900
         TOL = 500*FUZZ                                                  SLAP........322000
         IERR = 4                                                        SLAP........322100
      ENDIF                                                              SLAP........322200
      FUZZ = FUZZ*FUZZ                                                   SLAP........322300
C                                                                        SLAP........322400
C         Calculate initial residual and pseudo-residual, and check      SLAP........322500
C         stopping criterion.                                            SLAP........322600
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)                        SLAP........322700
      DO 10 I = 1, N                                                     SLAP........322800
         R(I)  = B(I) - R(I)                                             SLAP........322900
 10   CONTINUE                                                           SLAP........323000
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)          SLAP........323100
C                                                                        SLAP........323200
      IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,          SLAP........323300
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,                     SLAP........323400
     $     R, Z, P, AP, EMAP, DZ, CSAV,                                  SLAP........323500
     $     RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200            SLAP........323600
      IF( IERR.NE.0 ) RETURN                                             SLAP........323700
C                                                                        SLAP........323800
C                                                                        SLAP........323900
C         ***** iteration loop *****                                     SLAP........324000
C                                                                        SLAP........324100
CVD$R NOVECTOR                                                           SLAP........324200
CVD$R NOCONCUR                                                           SLAP........324300
      DO 100 K = 1, ITMAX                                                SLAP........324400
         ITER = K                                                        SLAP........324500
         IP = MOD( ITER-1, NSAVE+1 )                                     SLAP........324600
C                                                                        SLAP........324700
C         calculate direction vector p, a*p, and (m-inv)*a*p,            SLAP........324800
C         and save if desired.                                           SLAP........324900
         CALL DCOPY(N, Z, 1, P(1,IP), 1)                                 SLAP........325000
         CALL MATVEC(N, P(1,IP), AP(1,IP), NELT, IA, JA, A, ISYM)        SLAP........325100
         CALL MSOLVE(N, AP(1,IP), EMAP(1,IP), NELT, IA, JA, A, ISYM,     SLAP........325200
     $        RWORK, IWORK)                                              SLAP........325300
         IF( NSAVE.EQ.0 ) THEN                                           SLAP........325400
            AKDEN = DDOT(N, EMAP, 1, EMAP, 1)                            SLAP........325500
         ELSE                                                            SLAP........325600
            IF( ITER.GT.1 ) THEN                                         SLAP........325700
               LMAX = MIN( NSAVE, ITER-1 )                               SLAP........325800
               DO 20 L = 1, LMAX                                         SLAP........325900
                  IPO = MOD(IP+(NSAVE+1-L),NSAVE+1)                      SLAP........326000
                  BKL = DDOT(N, EMAP(1,IP), 1, EMAP(1,IPO), 1)           SLAP........326100
                  BKL = BKL*CSAV(L)                                      SLAP........326200
                  CALL DAXPY(N, -BKL,    P(1,IPO), 1,    P(1,IP), 1)     SLAP........326300
                  CALL DAXPY(N, -BKL,   AP(1,IPO), 1,   AP(1,IP), 1)     SLAP........326400
                  CALL DAXPY(N, -BKL, EMAP(1,IPO), 1, EMAP(1,IP), 1)     SLAP........326500
 20            CONTINUE                                                  SLAP........326600
               IF( NSAVE.GT.1 ) THEN                                     SLAP........326700
                  DO 30 L = NSAVE-1, 1, -1                               SLAP........326800
                     CSAV(L+1) = CSAV(L)                                 SLAP........326900
 30               CONTINUE                                               SLAP........327000
               ENDIF                                                     SLAP........327100
            ENDIF                                                        SLAP........327200
            AKDEN = DDOT(N, EMAP(1,IP), 1, EMAP(1,IP), 1)                SLAP........327300
            IF( ABS(AKDEN).LT.FUZZ ) THEN                                SLAP........327400
               IERR = 6                                                  SLAP........327500
               RETURN                                                    SLAP........327600
            ENDIF                                                        SLAP........327700
            CSAV(1) = 1.0D0/AKDEN                                        SLAP........327800
C                                                                        SLAP........327900
C         calculate coefficient ak, new iterate x, new residual r, and   SLAP........328000
C         new pseudo-residual z.                                         SLAP........328100
         ENDIF                                                           SLAP........328200
         AKNUM = DDOT(N, Z, 1, EMAP(1,IP), 1)                            SLAP........328300
         AK = AKNUM/AKDEN                                                SLAP........328400
         CALL DAXPY(N,  AK,    P(1,IP), 1, X, 1)                         SLAP........328500
         CALL DAXPY(N, -AK,   AP(1,IP), 1, R, 1)                         SLAP........328600
         CALL DAXPY(N, -AK, EMAP(1,IP), 1, Z, 1)                         SLAP........328700
C                                                                        SLAP........328800
C         check stopping criterion.                                      SLAP........328900
         IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,       SLAP........329000
     $        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,                  SLAP........329100
     $        R, Z, P, AP, EMAP, DZ, CSAV,                               SLAP........329200
     $        RWORK, IWORK, AK, BNRM, SOLNRM) .NE. 0 ) GO TO 200         SLAP........329300
C                                                                        SLAP........329400
 100  CONTINUE                                                           SLAP........329500
C                                                                        SLAP........329600
C         *****   end of loop  *****                                     SLAP........329700
C                                                                        SLAP........329800
C         Stopping criterion not satisfied.                              SLAP........329900
C.....THE COMMENT LINE MARKED WITH "CCC" BELOW WAS ORIGINALLY ACTIVE     SLAP........330000
C        CODE AND WAS COMMENTED OUT DURING INTEGRATION OF SLAP           SLAP........330100
C        WITH SUTRA.                                                     SLAP........330200
CCC   ITER = ITMAX + 1                                                   SLAP........330300
      IERR = 2                                                           SLAP........330400
C                                                                        SLAP........330500
 200  RETURN                                                             SLAP........330600
C------------- LAST LINE OF DOMN FOLLOWS ----------------------------    SLAP........330700
      END                                                                SLAP........330800
*DECK DORTH                                                              SLAP........330900
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)         SLAP........331000
C***BEGIN PROLOGUE  DORTH                                                SLAP........331100
C***SUBSIDIARY                                                           SLAP........331200
C***PURPOSE  Internal routine for DGMRES.                                SLAP........331300
C***LIBRARY   SLATEC (SLAP)                                              SLAP........331400
C***CATEGORY  D2A4, D2B4                                                 SLAP........331500
C***TYPE      DOUBLE PRECISION (SORTH-S, DORTH-D)                        SLAP........331600
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........331700
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........331800
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........331900
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........332000
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........332100
C             Lawrence Livermore National Laboratory                     SLAP........332200
C             PO Box 808, L-60                                           SLAP........332300
C             Livermore, CA 94550 (510) 423-3141                         SLAP........332400
C***DESCRIPTION                                                          SLAP........332500
C        This routine  orthogonalizes  the  vector  VNEW  against the    SLAP........332600
C        previous KMP  vectors in the   V array.  It uses  a modified    SLAP........332700
C        Gram-Schmidt   orthogonalization procedure with  conditional    SLAP........332800
C        reorthogonalization.                                            SLAP........332900
C                                                                        SLAP........333000
C *Usage:                                                                SLAP........333100
C      INTEGER N, LL, LDHES, KMP                                         SLAP........333200
C      DOUBLE PRECISION VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW          SLAP........333300
C                                                                        SLAP........333400
C      CALL DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)               SLAP........333500
C                                                                        SLAP........333600
C *Arguments:                                                            SLAP........333700
C VNEW   :INOUT    Double Precision VNEW(N)                              SLAP........333800
C         On input, the vector of length N containing a scaled           SLAP........333900
C         product of the Jacobian and the vector V(*,LL).                SLAP........334000
C         On output, the new vector orthogonal to V(*,i0) to V(*,LL),    SLAP........334100
C         where i0 = max(1, LL-KMP+1).                                   SLAP........334200
C V      :IN       Double Precision V(N,LL)                              SLAP........334300
C         The N x LL array containing the previous LL                    SLAP........334400
C         orthogonal vectors V(*,1) to V(*,LL).                          SLAP........334500
C HES    :INOUT    Double Precision HES(LDHES,LL)                        SLAP........334600
C         On input, an LL x LL upper Hessenberg matrix containing,       SLAP........334700
C         in HES(I,K), K.lt.LL, the scaled inner products of             SLAP........334800
C         A*V(*,K) and V(*,i).                                           SLAP........334900
C         On return, column LL of HES is filled in with                  SLAP........335000
C         the scaled inner products of A*V(*,LL) and V(*,i).             SLAP........335100
C N      :IN       Integer                                               SLAP........335200
C         The order of the matrix A, and the length of VNEW.             SLAP........335300
C LL     :IN       Integer                                               SLAP........335400
C         The current order of the matrix HES.                           SLAP........335500
C LDHES  :IN       Integer                                               SLAP........335600
C         The leading dimension of the HES array.                        SLAP........335700
C KMP    :IN       Integer                                               SLAP........335800
C         The number of previous vectors the new vector VNEW             SLAP........335900
C         must be made orthogonal to (KMP .le. MAXL).                    SLAP........336000
C SNORMW :OUT      DOUBLE PRECISION                                      SLAP........336100
C         Scalar containing the l-2 norm of VNEW.                        SLAP........336200
C                                                                        SLAP........336300
C***SEE ALSO  DGMRES                                                     SLAP........336400
C***ROUTINES CALLED  DAXPY, DDOT, DNRM2                                  SLAP........336500
C***REVISION HISTORY  (YYMMDD)                                           SLAP........336600
C   890404  DATE WRITTEN                                                 SLAP........336700
C   890404  Previous REVISION DATE                                       SLAP........336800
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........336900
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........337000
C           standard.  (FNF)                                             SLAP........337100
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........337200
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........337300
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........337400
C   920511  Added complete declaration section.  (WRB)                   SLAP........337500
C***END PROLOGUE  DORTH                                                  SLAP........337600
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........337700
CLLL. OPTIMIZE                                                           SLAP........337800
C     .. Scalar Arguments ..                                             SLAP........337900
      DOUBLE PRECISION SNORMW                                            SLAP........338000
      INTEGER KMP, LDHES, LL, N                                          SLAP........338100
C     .. Array Arguments ..                                              SLAP........338200
      DOUBLE PRECISION HES(LDHES,*), V(N,*), VNEW(*)                     SLAP........338300
C     .. Local Scalars ..                                                SLAP........338400
      DOUBLE PRECISION ARG, SUMDSQ, TEM, VNRM                            SLAP........338500
      INTEGER I, I0                                                      SLAP........338600
C     .. External Functions ..                                           SLAP........338700
      DOUBLE PRECISION DDOT, DNRM2                                       SLAP........338800
      EXTERNAL DDOT, DNRM2                                               SLAP........338900
C     .. External Subroutines ..                                         SLAP........339000
      EXTERNAL DAXPY                                                     SLAP........339100
C     .. Intrinsic Functions ..                                          SLAP........339200
      INTRINSIC MAX, SQRT                                                SLAP........339300
C***FIRST EXECUTABLE STATEMENT  DORTH                                    SLAP........339400
C                                                                        SLAP........339500
C         Get norm of unaltered VNEW for later use.                      SLAP........339600
C                                                                        SLAP........339700
      VNRM = DNRM2(N, VNEW, 1)                                           SLAP........339800
C   -------------------------------------------------------------------  SLAP........339900
C         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).  SLAP........340000
C         Scaled inner products give new column of HES.                  SLAP........340100
C         Projections of earlier vectors are subtracted from VNEW.       SLAP........340200
C   -------------------------------------------------------------------  SLAP........340300
      I0 = MAX(1,LL-KMP+1)                                               SLAP........340400
      DO 10 I = I0,LL                                                    SLAP........340500
         HES(I,LL) = DDOT(N, V(1,I), 1, VNEW, 1)                         SLAP........340600
         TEM = -HES(I,LL)                                                SLAP........340700
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)                          SLAP........340800
 10   CONTINUE                                                           SLAP........340900
C   -------------------------------------------------------------------  SLAP........341000
C         Compute SNORMW = norm of VNEW.  If VNEW is small compared      SLAP........341100
C         to its input value (in norm), then reorthogonalize VNEW to     SLAP........341200
C         V(*,1) through V(*,LL).  Correct if relative correction        SLAP........341300
C         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using   SLAP........341400
C         the dot products involved.                                     SLAP........341500
C   -------------------------------------------------------------------  SLAP........341600
      SNORMW = DNRM2(N, VNEW, 1)                                         SLAP........341700
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN                        SLAP........341800
      SUMDSQ = 0                                                         SLAP........341900
      DO 30 I = I0,LL                                                    SLAP........342000
         TEM = -DDOT(N, V(1,I), 1, VNEW, 1)                              SLAP........342100
         IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30            SLAP........342200
         HES(I,LL) = HES(I,LL) - TEM                                     SLAP........342300
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)                          SLAP........342400
         SUMDSQ = SUMDSQ + TEM**2                                        SLAP........342500
 30   CONTINUE                                                           SLAP........342600
      IF (SUMDSQ .EQ. 0.0D0) RETURN                                      SLAP........342700
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)                                SLAP........342800
      SNORMW = SQRT(ARG)                                                 SLAP........342900
C                                                                        SLAP........343000
      RETURN                                                             SLAP........343100
C------------- LAST LINE OF DORTH FOLLOWS ----------------------------   SLAP........343200
      END                                                                SLAP........343300
*DECK DPIGMR                                                             SLAP........343400
      SUBROUTINE DPIGMR (N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS, SLAP........343500
     +   JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK, SLAP........343600
     +   DL, RHOL, NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A,    SLAP........343700
C........THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION OF        SLAP........343800
C           SLAP WITH SUTRA.  IT ORIGINALLY READ:                        SLAP........343900
C           +   ISYM, IUNIT, IFLAG, ERR)                                 SLAP........344000
     +   ISYM, IUNIT, IFLAG, ERR, ITMAX)                                 SLAP........344100
C***BEGIN PROLOGUE  DPIGMR                                               SLAP........344200
C***SUBSIDIARY                                                           SLAP........344300
C***PURPOSE  Internal routine for DGMRES.                                SLAP........344400
C***LIBRARY   SLATEC (SLAP)                                              SLAP........344500
C***CATEGORY  D2A4, D2B4                                                 SLAP........344600
C***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D)                      SLAP........344700
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........344800
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........344900
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........345000
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........345100
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........345200
C             Lawrence Livermore National Laboratory                     SLAP........345300
C             PO Box 808, L-60                                           SLAP........345400
C             Livermore, CA 94550 (510) 423-3141                         SLAP........345500
C***DESCRIPTION                                                          SLAP........345600
C         This routine solves the linear system A * Z = R0 using a       SLAP........345700
C         scaled preconditioned version of the generalized minimum       SLAP........345800
C         residual method.  An initial guess of Z = 0 is assumed.        SLAP........345900
C                                                                        SLAP........346000
C *Usage:                                                                SLAP........346100
C      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR      SLAP........346200
C      INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT) SLAP........346300
C      INTEGER ISYM, IUNIT, IFLAG                                        SLAP........346400
C      DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1),          SLAP........346500
C     $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED), SLAP........346600
C     $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N),     SLAP........346700
C     $                 TOL, A(NELT), ERR                                SLAP........346800
C      EXTERNAL MATVEC, MSOLVE                                           SLAP........346900
C                                                                        SLAP........347000
C      CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP,              SLAP........347100
C     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,       SLAP........347200
C     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,             SLAP........347300
CC..........THE NEXT COMMENT LINE *NOT* MARKED WITH "CC" WAS MODIFIED    SLAP........347400
CC             DURING INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY     SLAP........347500
CC             READ:                                                     SLAP........347600
CC    C     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)   SLAP........347700
C     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR, ITMAX)  SLAP........347800
C                                                                        SLAP........347900
C *Arguments:                                                            SLAP........348000
C N      :IN       Integer                                               SLAP........348100
C         The order of the matrix A, and the lengths                     SLAP........348200
C         of the vectors SR, SZ, R0 and Z.                               SLAP........348300
C R0     :IN       Double Precision R0(N)                                SLAP........348400
C         R0 = the right hand side of the system A*Z = R0.               SLAP........348500
C         R0 is also used as workspace when computing                    SLAP........348600
C         the final approximation.                                       SLAP........348700
C         (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.)         SLAP........348800
C SR     :IN       Double Precision SR(N)                                SLAP........348900
C         SR is a vector of length N containing the non-zero             SLAP........349000
C         elements of the diagonal scaling matrix for R0.                SLAP........349100
C SZ     :IN       Double Precision SZ(N)                                SLAP........349200
C         SZ is a vector of length N containing the non-zero             SLAP........349300
C         elements of the diagonal scaling matrix for Z.                 SLAP........349400
C JSCAL  :IN       Integer                                               SLAP........349500
C         A flag indicating whether arrays SR and SZ are used.           SLAP........349600
C         JSCAL=0 means SR and SZ are not used and the                   SLAP........349700
C                 algorithm will perform as if all                       SLAP........349800
C                 SR(i) = 1 and SZ(i) = 1.                               SLAP........349900
C         JSCAL=1 means only SZ is used, and the algorithm               SLAP........350000
C                 performs as if all SR(i) = 1.                          SLAP........350100
C         JSCAL=2 means only SR is used, and the algorithm               SLAP........350200
C                 performs as if all SZ(i) = 1.                          SLAP........350300
C         JSCAL=3 means both SR and SZ are used.                         SLAP........350400
C MAXL   :IN       Integer                                               SLAP........350500
C         The maximum allowable order of the matrix H.                   SLAP........350600
C MAXLP1 :IN       Integer                                               SLAP........350700
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.       SLAP........350800
C KMP    :IN       Integer                                               SLAP........350900
C         The number of previous vectors the new vector VNEW             SLAP........351000
C         must be made orthogonal to.  (KMP .le. MAXL)                   SLAP........351100
C NRSTS  :IN       Integer                                               SLAP........351200
C         Counter for the number of restarts on the current              SLAP........351300
C         call to DGMRES.  If NRSTS .gt. 0, then the residual            SLAP........351400
C         R0 is already scaled, and so scaling of it is                  SLAP........351500
C         not necessary.                                                 SLAP........351600
C JPRE   :IN       Integer                                               SLAP........351700
C         Preconditioner type flag.                                      SLAP........351800
C MATVEC :EXT      External.                                             SLAP........351900
C         Name of a routine which performs the matrix vector multiply    SLAP........352000
C         Y = A*X given A and X.  The name of the MATVEC routine must    SLAP........352100
C         be declared external in the calling program.  The calling      SLAP........352200
C         sequence to MATVEC is:                                         SLAP........352300
C             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)                SLAP........352400
C         where N is the number of unknowns, Y is the product A*X        SLAP........352500
C         upon return, X is an input vector, and NELT is the number of   SLAP........352600
C         non-zeros in the SLAP IA, JA, A storage for the matrix A.      SLAP........352700
C         ISYM is a flag which, if non-zero, denotes that A is           SLAP........352800
C         symmetric and only the lower or upper triangle is stored.      SLAP........352900
C MSOLVE :EXT      External.                                             SLAP........353000
C         Name of the routine which solves a linear system Mz = r for    SLAP........353100
C         z given r with the preconditioning matrix M (M is supplied via SLAP........353200
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must     SLAP........353300
C         be declared external in the calling program.  The calling      SLAP........353400
C         sequence to MSOLVE is:                                         SLAP........353500
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)    SLAP........353600
C         Where N is the number of unknowns, R is the right-hand side    SLAP........353700
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........353800
C         ISYM are defined as below.  RPAR is a double precision array   SLAP........353900
C         that can be used to pass necessary preconditioning information SLAP........354000
C         and/or workspace to MSOLVE.  IPAR is an integer work array     SLAP........354100
C         for the same purpose as RPAR.                                  SLAP........354200
C NMSL   :OUT      Integer                                               SLAP........354300
C         The number of calls to MSOLVE.                                 SLAP........354400
C Z      :OUT      Double Precision Z(N)                                 SLAP........354500
C         The final computed approximation to the solution               SLAP........354600
C         of the system A*Z = R0.                                        SLAP........354700
C V      :OUT      Double Precision V(N,MAXLP1)                          SLAP........354800
C         The N by (LGMR+1) array containing the LGMR                    SLAP........354900
C         orthogonal vectors V(*,1) to V(*,LGMR).                        SLAP........355000
C HES    :OUT      Double Precision HES(MAXLP1,MAXL)                     SLAP........355100
C         The upper triangular factor of the QR decomposition            SLAP........355200
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose          SLAP........355300
C         entries are the scaled inner-products of A*V(*,I)              SLAP........355400
C         and V(*,K).                                                    SLAP........355500
C Q      :OUT      Double Precision Q(2*MAXL)                            SLAP........355600
C         A double precision array of length 2*MAXL containing the       SLAP........355700
C         components of the Givens rotations used in the QR              SLAP........355800
C         decomposition of HES.  It is loaded in DHEQR and used in       SLAP........355900
C         DHELS.                                                         SLAP........356000
C LGMR   :OUT      Integer                                               SLAP........356100
C         The number of iterations performed and                         SLAP........356200
C         the current order of the upper Hessenberg                      SLAP........356300
C         matrix HES.                                                    SLAP........356400
C RPAR   :IN       Double Precision RPAR(USER DEFINED)                   SLAP........356500
C         Double Precision workspace passed directly to the MSOLVE       SLAP........356600
C         routine.                                                       SLAP........356700
C IPAR   :IN       Integer IPAR(USER DEFINED)                            SLAP........356800
C         Integer workspace passed directly to the MSOLVE routine.       SLAP........356900
C WK     :IN       Double Precision WK(N)                                SLAP........357000
C         A double precision work array of length N used by routines     SLAP........357100
C         MATVEC and MSOLVE.                                             SLAP........357200
C DL     :INOUT    Double Precision DL(N)                                SLAP........357300
C         On input, a double precision work array of length N used for   SLAP........357400
C         calculation of the residual norm RHO when the method is        SLAP........357500
C         incomplete (KMP.lt.MAXL), and/or when using restarting.        SLAP........357600
C         On output, the scaled residual vector RL.  It is only loaded   SLAP........357700
C         when performing restarts of the Krylov iteration.              SLAP........357800
C RHOL   :OUT      Double Precision                                      SLAP........357900
C         A double precision scalar containing the norm of the final     SLAP........358000
C         residual.                                                      SLAP........358100
C NRMAX  :IN       Integer                                               SLAP........358200
C         The maximum number of restarts of the Krylov iteration.        SLAP........358300
C         NRMAX .gt. 0 means restarting is active, while                 SLAP........358400
C         NRMAX = 0 means restarting is not being used.                  SLAP........358500
C B      :IN       Double Precision B(N)                                 SLAP........358600
C         The right hand side of the linear system A*X = b.              SLAP........358700
C BNRM   :IN       Double Precision                                      SLAP........358800
C         The scaled norm of b.                                          SLAP........358900
C X      :IN       Double Precision X(N)                                 SLAP........359000
C         The current approximate solution as of the last                SLAP........359100
C         restart.                                                       SLAP........359200
C XL     :IN       Double Precision XL(N)                                SLAP........359300
C         An array of length N used to hold the approximate              SLAP........359400
C         solution X(L) when ITOL=11.                                    SLAP........359500
C ITOL   :IN       Integer                                               SLAP........359600
C         A flag to indicate the type of convergence criterion           SLAP........359700
C         used.  See the driver for its description.                     SLAP........359800
C TOL    :IN       Double Precision                                      SLAP........359900
C         The tolerance on residuals R0-A*Z in scaled norm.              SLAP........360000
C NELT   :IN       Integer                                               SLAP........360100
C         The length of arrays IA, JA and A.                             SLAP........360200
C IA     :IN       Integer IA(NELT)                                      SLAP........360300
C         An integer array of length NELT containing matrix data.        SLAP........360400
C         It is passed directly to the MATVEC and MSOLVE routines.       SLAP........360500
C JA     :IN       Integer JA(NELT)                                      SLAP........360600
C         An integer array of length NELT containing matrix data.        SLAP........360700
C         It is passed directly to the MATVEC and MSOLVE routines.       SLAP........360800
C A      :IN       Double Precision A(NELT)                              SLAP........360900
C         A double precision array of length NELT containing matrix      SLAP........361000
C         data. It is passed directly to the MATVEC and MSOLVE routines. SLAP........361100
C ISYM   :IN       Integer                                               SLAP........361200
C         A flag to indicate symmetric matrix storage.                   SLAP........361300
C         If ISYM=0, all non-zero entries of the matrix are              SLAP........361400
C         stored.  If ISYM=1, the matrix is symmetric and                SLAP........361500
C         only the upper or lower triangular part is stored.             SLAP........361600
C IUNIT  :IN       Integer                                               SLAP........361700
C         The i/o unit number for writing intermediate residual          SLAP........361800
C         norm values.                                                   SLAP........361900
C IFLAG  :OUT      Integer                                               SLAP........362000
C         An integer error flag..                                        SLAP........362100
C         0 means convergence in LGMR iterations, LGMR.le.MAXL.          SLAP........362200
C         1 means the convergence test did not pass in MAXL              SLAP........362300
C           iterations, but the residual norm is .lt. norm(R0),          SLAP........362400
C           and so Z is computed.                                        SLAP........362500
C         2 means the convergence test did not pass in MAXL              SLAP........362600
C           iterations, residual .ge. norm(R0), and Z = 0.               SLAP........362700
C ERR    :OUT      Double Precision.                                     SLAP........362800
C         Error estimate of error in final approximate solution, as      SLAP........362900
C         defined by ITOL.                                               SLAP........363000
C                                                                        SLAP........363100
C *Cautions:                                                             SLAP........363200
C     This routine will attempt to write to the Fortran logical output   SLAP........363300
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........363400
C     this logical unit is attached to a file or terminal before calling SLAP........363500
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........363600
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........363700
C                                                                        SLAP........363800
C***SEE ALSO  DGMRES                                                     SLAP........363900
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL,   SLAP........364000
C                    DSCAL, ISDGMR                                       SLAP........364100
C***REVISION HISTORY  (YYMMDD)                                           SLAP........364200
C   890404  DATE WRITTEN                                                 SLAP........364300
C   890404  Previous REVISION DATE                                       SLAP........364400
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........364500
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........364600
C           standard.  (FNF)                                             SLAP........364700
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........364800
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........364900
C   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)  SLAP........365000
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........365100
C   920511  Added complete declaration section.  (WRB)                   SLAP........365200
C***END PROLOGUE  DPIGMR                                                 SLAP........365300
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........365400
CLLL. OPTIMIZE                                                           SLAP........365500
C     .. Scalar Arguments ..                                             SLAP........365600
      DOUBLE PRECISION BNRM, ERR, RHOL, TOL                              SLAP........365700
      INTEGER IFLAG, ISYM, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, MAXL,    SLAP........365800
     +        MAXLP1, N, NELT, NMSL, NRMAX, NRSTS                        SLAP........365900
C     .. Array Arguments ..                                              SLAP........366000
      DOUBLE PRECISION A(NELT), B(*), DL(*), HES(MAXLP1,*), Q(*), R0(*), SLAP........366100
     +                 RPAR(*), SR(*), SZ(*), V(N,*), WK(*), X(*),       SLAP........366200
     +                 XL(*), Z(*)                                       SLAP........366300
      INTEGER IA(NELT), IPAR(*), JA(NELT)                                SLAP........366400
C     .. Subroutine Arguments ..                                         SLAP........366500
      EXTERNAL MATVEC, MSOLVE                                            SLAP........366600
C     .. Local Scalars ..                                                SLAP........366700
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM        SLAP........366800
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1              SLAP........366900
C     .. External Functions ..                                           SLAP........367000
      DOUBLE PRECISION DNRM2                                             SLAP........367100
      INTEGER ISDGMR                                                     SLAP........367200
      EXTERNAL DNRM2, ISDGMR                                             SLAP........367300
C     .. External Subroutines ..                                         SLAP........367400
      EXTERNAL DAXPY, DCOPY, DHELS, DHEQR, DORTH, DRLCAL, DSCAL          SLAP........367500
C     .. Intrinsic Functions ..                                          SLAP........367600
      INTRINSIC ABS                                                      SLAP........367700
C***FIRST EXECUTABLE STATEMENT  DPIGMR                                   SLAP........367800
C                                                                        SLAP........367900
C         Zero out the Z array.                                          SLAP........368000
C                                                                        SLAP........368100
      DO 5 I = 1,N                                                       SLAP........368200
         Z(I) = 0                                                        SLAP........368300
 5    CONTINUE                                                           SLAP........368400
C                                                                        SLAP........368500
      IFLAG = 0                                                          SLAP........368600
      LGMR = 0                                                           SLAP........368700
      NMSL = 0                                                           SLAP........368800
C         Load ITMAX, the maximum number of iterations.                  SLAP........368900
C.....THE COMMENT LINE MARKED "CCC" BELOW WAS ORIGINALLY ACTIVE CODE     SLAP........369000
C        AND WAS COMMENTED OUT DURING INTEGRATION OF SLAP WITH SUTRA.    SLAP........369100
CCC   ITMAX =(NRMAX+1)*MAXL                                              SLAP........369200
C   -------------------------------------------------------------------  SLAP........369300
C         The initial residual is the vector R0.                         SLAP........369400
C         Apply left precon. if JPRE < 0 and this is not a restart.      SLAP........369500
C         Apply scaling to R0 if JSCAL = 2 or 3.                         SLAP........369600
C   -------------------------------------------------------------------  SLAP........369700
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN                        SLAP........369800
         CALL DCOPY(N, R0, 1, WK, 1)                                     SLAP........369900
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)       SLAP........370000
         NMSL = NMSL + 1                                                 SLAP........370100
      ENDIF                                                              SLAP........370200
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN        SLAP........370300
         DO 10 I = 1,N                                                   SLAP........370400
            V(I,1) = R0(I)*SR(I)                                         SLAP........370500
 10      CONTINUE                                                        SLAP........370600
      ELSE                                                               SLAP........370700
         DO 20 I = 1,N                                                   SLAP........370800
            V(I,1) = R0(I)                                               SLAP........370900
 20      CONTINUE                                                        SLAP........371000
      ENDIF                                                              SLAP........371100
      R0NRM = DNRM2(N, V, 1)                                             SLAP........371200
      ITER = NRSTS*MAXL                                                  SLAP........371300
C                                                                        SLAP........371400
C         Call stopping routine ISDGMR.                                  SLAP........371500
C                                                                        SLAP........371600
      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,             SLAP........371700
     $    NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK,       SLAP........371800
     $    RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL,                        SLAP........371900
     $    KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,            SLAP........372000
     $    HES, JPRE) .NE. 0) RETURN                                      SLAP........372100
      TEM = 1.0D0/R0NRM                                                  SLAP........372200
      CALL DSCAL(N, TEM, V(1,1), 1)                                      SLAP........372300
C                                                                        SLAP........372400
C         Zero out the HES array.                                        SLAP........372500
C                                                                        SLAP........372600
      DO 50 J = 1,MAXL                                                   SLAP........372700
         DO 40 I = 1,MAXLP1                                              SLAP........372800
            HES(I,J) = 0                                                 SLAP........372900
 40      CONTINUE                                                        SLAP........373000
 50   CONTINUE                                                           SLAP........373100
C   -------------------------------------------------------------------  SLAP........373200
C         Main loop to compute the vectors V(*,2) to V(*,MAXL).          SLAP........373300
C         The running product PROD is needed for the convergence test.   SLAP........373400
C   -------------------------------------------------------------------  SLAP........373500
      PROD = 1                                                           SLAP........373600
      DO 90 LL = 1,MAXL                                                  SLAP........373700
         LGMR = LL                                                       SLAP........373800
C   -------------------------------------------------------------------  SLAP........373900
C        Unscale  the  current V(LL)  and store  in WK.  Call routine    SLAP........374000
C        MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the    SLAP........374100
C        preconditioner matrix.  Save the answer in Z.   Call routine    SLAP........374200
C        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system    SLAP........374300
C        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call    SLAP........374400
C        routine DORTH  to  orthogonalize the    new vector VNEW   =     SLAP........374500
C        V(*,LL+1).  Call routine DHEQR to update the factors of HES.    SLAP........374600
C   -------------------------------------------------------------------  SLAP........374700
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN                      SLAP........374800
           DO 60 I = 1,N                                                 SLAP........374900
              WK(I) = V(I,LL)/SZ(I)                                      SLAP........375000
 60        CONTINUE                                                      SLAP........375100
        ELSE                                                             SLAP........375200
           CALL DCOPY(N, V(1,LL), 1, WK, 1)                              SLAP........375300
        ENDIF                                                            SLAP........375400
        IF (JPRE .GT. 0) THEN                                            SLAP........375500
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)      SLAP........375600
           NMSL = NMSL + 1                                               SLAP........375700
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)           SLAP........375800
        ELSE                                                             SLAP........375900
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)          SLAP........376000
        ENDIF                                                            SLAP........376100
        IF (JPRE .LT. 0) THEN                                            SLAP........376200
           CALL DCOPY(N, V(1,LL+1), 1, WK, 1)                            SLAP........376300
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)       SLAP........376400
           NMSL = NMSL + 1                                               SLAP........376500
        ENDIF                                                            SLAP........376600
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN                      SLAP........376700
           DO 65 I = 1,N                                                 SLAP........376800
              V(I,LL+1) = V(I,LL+1)*SR(I)                                SLAP........376900
 65        CONTINUE                                                      SLAP........377000
        ENDIF                                                            SLAP........377100
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)        SLAP........377200
        HES(LL+1,LL) = SNORMW                                            SLAP........377300
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)                         SLAP........377400
        IF (INFO .EQ. LL) GO TO 120                                      SLAP........377500
C   -------------------------------------------------------------------  SLAP........377600
C         Update RHO, the estimate of the norm of the residual R0-A*ZL.  SLAP........377700
C         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not  SLAP........377800
C         necessarily orthogonal for LL > KMP.  The vector DL must then  SLAP........377900
C         be computed, and its norm used in the calculation of RHO.      SLAP........378000
C   -------------------------------------------------------------------  SLAP........378100
        PROD = PROD*Q(2*LL)                                              SLAP........378200
        RHO = ABS(PROD*R0NRM)                                            SLAP........378300
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN                         SLAP........378400
           IF (LL .EQ. KMP+1) THEN                                       SLAP........378500
              CALL DCOPY(N, V(1,1), 1, DL, 1)                            SLAP........378600
              DO 75 I = 1,KMP                                            SLAP........378700
                 IP1 = I + 1                                             SLAP........378800
                 I2 = I*2                                                SLAP........378900
                 S = Q(I2)                                               SLAP........379000
                 C = Q(I2-1)                                             SLAP........379100
                 DO 70 K = 1,N                                           SLAP........379200
                    DL(K) = S*DL(K) + C*V(K,IP1)                         SLAP........379300
 70              CONTINUE                                                SLAP........379400
 75           CONTINUE                                                   SLAP........379500
           ENDIF                                                         SLAP........379600
           S = Q(2*LL)                                                   SLAP........379700
           C = Q(2*LL-1)/SNORMW                                          SLAP........379800
           LLP1 = LL + 1                                                 SLAP........379900
           DO 80 K = 1,N                                                 SLAP........380000
              DL(K) = S*DL(K) + C*V(K,LLP1)                              SLAP........380100
 80        CONTINUE                                                      SLAP........380200
           DLNRM = DNRM2(N, DL, 1)                                       SLAP........380300
           RHO = RHO*DLNRM                                               SLAP........380400
        ENDIF                                                            SLAP........380500
        RHOL = RHO                                                       SLAP........380600
C   -------------------------------------------------------------------  SLAP........380700
C         Test for convergence.  If passed, compute approximation ZL.    SLAP........380800
C         If failed and LL < MAXL, then continue iterating.              SLAP........380900
C   -------------------------------------------------------------------  SLAP........381000
        ITER = NRSTS*MAXL + LGMR                                         SLAP........381100
        IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,           SLAP........381200
     $      NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK,         SLAP........381300
     $      RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL,                       SLAP........381400
     $      KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,          SLAP........381500
     $      HES, JPRE) .NE. 0) GO TO 200                                 SLAP........381600
C.......THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION OF SLAP    SLAP........381700
C          WITH SUTRA.  IT ORGINALLY READ:                               SLAP........381800
C          IF (LL .EQ. MAXL) GO TO 100                                   SLAP........381900
        IF ((ITER.EQ.ITMAX).OR.(LL .EQ. MAXL)) GO TO 100                 SLAP........382000
C   -------------------------------------------------------------------  SLAP........382100
C         Rescale so that the norm of V(1,LL+1) is one.                  SLAP........382200
C   -------------------------------------------------------------------  SLAP........382300
        TEM = 1.0D0/SNORMW                                               SLAP........382400
        CALL DSCAL(N, TEM, V(1,LL+1), 1)                                 SLAP........382500
 90   CONTINUE                                                           SLAP........382600
 100  CONTINUE                                                           SLAP........382700
      IF (RHO .LT. R0NRM) GO TO 150                                      SLAP........382800
 120  CONTINUE                                                           SLAP........382900
      IFLAG = 2                                                          SLAP........383000
C                                                                        SLAP........383100
C         Load approximate solution with zero.                           SLAP........383200
C                                                                        SLAP........383300
      DO 130 I = 1,N                                                     SLAP........383400
         Z(I) = 0                                                        SLAP........383500
 130  CONTINUE                                                           SLAP........383600
      RETURN                                                             SLAP........383700
 150  IFLAG = 1                                                          SLAP........383800
C                                                                        SLAP........383900
C         Tolerance not met, but residual norm reduced.                  SLAP........384000
C                                                                        SLAP........384100
      IF (NRMAX .GT. 0) THEN                                             SLAP........384200
C                                                                        SLAP........384300
C        If performing restarting (NRMAX > 0)  calculate the residual    SLAP........384400
C        vector RL and  store it in the DL  array.  If the incomplete    SLAP........384500
C        version is being used (KMP < MAXL) then DL has  already been    SLAP........384600
C        calculated up to a scaling factor.   Use DRLCAL to calculate    SLAP........384700
C        the scaled residual vector.                                     SLAP........384800
C                                                                        SLAP........384900
         CALL DRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD,         SLAP........385000
     $        R0NRM)                                                     SLAP........385100
      ENDIF                                                              SLAP........385200
C   -------------------------------------------------------------------  SLAP........385300
C         Compute the approximation ZL to the solution.  Since the       SLAP........385400
C         vector Z was used as workspace, and the initial guess          SLAP........385500
C         of the linear iteration is zero, Z must be reset to zero.      SLAP........385600
C   -------------------------------------------------------------------  SLAP........385700
 200  CONTINUE                                                           SLAP........385800
      LL = LGMR                                                          SLAP........385900
      LLP1 = LL + 1                                                      SLAP........386000
      DO 210 K = 1,LLP1                                                  SLAP........386100
         R0(K) = 0                                                       SLAP........386200
 210  CONTINUE                                                           SLAP........386300
      R0(1) = R0NRM                                                      SLAP........386400
      CALL DHELS(HES, MAXLP1, LL, Q, R0)                                 SLAP........386500
      DO 220 K = 1,N                                                     SLAP........386600
         Z(K) = 0                                                        SLAP........386700
 220  CONTINUE                                                           SLAP........386800
      DO 230 I = 1,LL                                                    SLAP........386900
         CALL DAXPY(N, R0(I), V(1,I), 1, Z, 1)                           SLAP........387000
 230  CONTINUE                                                           SLAP........387100
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN                        SLAP........387200
         DO 240 I = 1,N                                                  SLAP........387300
            Z(I) = Z(I)/SZ(I)                                            SLAP........387400
 240     CONTINUE                                                        SLAP........387500
      ENDIF                                                              SLAP........387600
      IF (JPRE .GT. 0) THEN                                              SLAP........387700
         CALL DCOPY(N, Z, 1, WK, 1)                                      SLAP........387800
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)        SLAP........387900
         NMSL = NMSL + 1                                                 SLAP........388000
      ENDIF                                                              SLAP........388100
      RETURN                                                             SLAP........388200
C------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------  SLAP........388300
      END                                                                SLAP........388400
*DECK DRLCAL                                                             SLAP........388500
      SUBROUTINE DRLCAL (N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD,       SLAP........388600
     +   R0NRM)                                                          SLAP........388700
C***BEGIN PROLOGUE  DRLCAL                                               SLAP........388800
C***SUBSIDIARY                                                           SLAP........388900
C***PURPOSE  Internal routine for DGMRES.                                SLAP........389000
C***LIBRARY   SLATEC (SLAP)                                              SLAP........389100
C***CATEGORY  D2A4, D2B4                                                 SLAP........389200
C***TYPE      DOUBLE PRECISION (SRLCAL-S, DRLCAL-D)                      SLAP........389300
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........389400
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........389500
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........389600
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........389700
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........389800
C             Lawrence Livermore National Laboratory                     SLAP........389900
C             PO Box 808, L-60                                           SLAP........390000
C             Livermore, CA 94550 (510) 423-3141                         SLAP........390100
C***DESCRIPTION                                                          SLAP........390200
C         This routine calculates the scaled residual RL from the        SLAP........390300
C         V(I)'s.                                                        SLAP........390400
C *Usage:                                                                SLAP........390500
C      INTEGER N, KMP, LL, MAXL                                          SLAP........390600
C      DOUBLE PRECISION V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM  SLAP........390700
C                                                                        SLAP........390800
C      CALL DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM)      SLAP........390900
C                                                                        SLAP........391000
C *Arguments:                                                            SLAP........391100
C N      :IN       Integer                                               SLAP........391200
C         The order of the matrix A, and the lengths                     SLAP........391300
C         of the vectors SR, SZ, R0 and Z.                               SLAP........391400
C KMP    :IN       Integer                                               SLAP........391500
C         The number of previous V vectors the new vector VNEW           SLAP........391600
C         must be made orthogonal to. (KMP .le. MAXL)                    SLAP........391700
C LL     :IN       Integer                                               SLAP........391800
C         The current dimension of the Krylov subspace.                  SLAP........391900
C MAXL   :IN       Integer                                               SLAP........392000
C         The maximum dimension of the Krylov subspace.                  SLAP........392100
C V      :IN       Double Precision V(N,LL)                              SLAP........392200
C         The N x LL array containing the orthogonal vectors             SLAP........392300
C         V(*,1) to V(*,LL).                                             SLAP........392400
C Q      :IN       Double Precision Q(2*MAXL)                            SLAP........392500
C         A double precision array of length 2*MAXL containing the       SLAP........392600
C         components of the Givens rotations used in the QR              SLAP........392700
C         decomposition of HES.  It is loaded in DHEQR and used in       SLAP........392800
C         DHELS.                                                         SLAP........392900
C RL     :OUT      Double Precision RL(N)                                SLAP........393000
C         The residual vector RL.  This is either SB*(B-A*XL) if         SLAP........393100
C         not preconditioning or preconditioning on the right,           SLAP........393200
C         or SB*(M-inverse)*(B-A*XL) if preconditioning on the           SLAP........393300
C         left.                                                          SLAP........393400
C SNORMW :IN       Double Precision                                      SLAP........393500
C         Scale factor.                                                  SLAP........393600
C PROD   :IN       Double Precision                                      SLAP........393700
C         The product s1*s2*...*sl = the product of the sines of the     SLAP........393800
C         Givens rotations used in the QR factorization of               SLAP........393900
C         the Hessenberg matrix HES.                                     SLAP........394000
C R0NRM  :IN       Double Precision                                      SLAP........394100
C         The scaled norm of initial residual R0.                        SLAP........394200
C                                                                        SLAP........394300
C***SEE ALSO  DGMRES                                                     SLAP........394400
C***ROUTINES CALLED  DCOPY, DSCAL                                        SLAP........394500
C***REVISION HISTORY  (YYMMDD)                                           SLAP........394600
C   890404  DATE WRITTEN                                                 SLAP........394700
C   890404  Previous REVISION DATE                                       SLAP........394800
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........394900
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........395000
C           standard.  (FNF)                                             SLAP........395100
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........395200
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........395300
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........395400
C   920511  Added complete declaration section.  (WRB)                   SLAP........395500
C***END PROLOGUE  DRLCAL                                                 SLAP........395600
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........395700
CLLL. OPTIMIZE                                                           SLAP........395800
C     .. Scalar Arguments ..                                             SLAP........395900
      DOUBLE PRECISION PROD, R0NRM, SNORMW                               SLAP........396000
      INTEGER KMP, LL, MAXL, N                                           SLAP........396100
C     .. Array Arguments ..                                              SLAP........396200
      DOUBLE PRECISION Q(*), RL(N), V(N,*)                               SLAP........396300
C     .. Local Scalars ..                                                SLAP........396400
      DOUBLE PRECISION C, S, TEM                                         SLAP........396500
      INTEGER I, I2, IP1, K, LLM1, LLP1                                  SLAP........396600
C     .. External Subroutines ..                                         SLAP........396700
      EXTERNAL DCOPY, DSCAL                                              SLAP........396800
C***FIRST EXECUTABLE STATEMENT  DRLCAL                                   SLAP........396900
      IF (KMP .EQ. MAXL) THEN                                            SLAP........397000
C                                                                        SLAP........397100
C         calculate RL.  Start by copying V(*,1) into RL.                SLAP........397200
C                                                                        SLAP........397300
         CALL DCOPY(N, V(1,1), 1, RL, 1)                                 SLAP........397400
         LLM1 = LL - 1                                                   SLAP........397500
         DO 20 I = 1,LLM1                                                SLAP........397600
            IP1 = I + 1                                                  SLAP........397700
            I2 = I*2                                                     SLAP........397800
            S = Q(I2)                                                    SLAP........397900
            C = Q(I2-1)                                                  SLAP........398000
            DO 10 K = 1,N                                                SLAP........398100
               RL(K) = S*RL(K) + C*V(K,IP1)                              SLAP........398200
 10         CONTINUE                                                     SLAP........398300
 20      CONTINUE                                                        SLAP........398400
         S = Q(2*LL)                                                     SLAP........398500
         C = Q(2*LL-1)/SNORMW                                            SLAP........398600
         LLP1 = LL + 1                                                   SLAP........398700
         DO 30 K = 1,N                                                   SLAP........398800
            RL(K) = S*RL(K) + C*V(K,LLP1)                                SLAP........398900
 30      CONTINUE                                                        SLAP........399000
      ENDIF                                                              SLAP........399100
C                                                                        SLAP........399200
C         When KMP < MAXL, RL vector already partially calculated.       SLAP........399300
C         Scale RL by R0NRM*PROD to obtain the residual RL.              SLAP........399400
C                                                                        SLAP........399500
      TEM = R0NRM*PROD                                                   SLAP........399600
      CALL DSCAL(N, TEM, RL, 1)                                          SLAP........399700
      RETURN                                                             SLAP........399800
C------------- LAST LINE OF DRLCAL FOLLOWS ----------------------------  SLAP........399900
      END                                                                SLAP........400000
*DECK DS2Y                                                               SLAP........400100
      SUBROUTINE DS2Y (N, NELT, IA, JA, A, ISYM)                         SLAP........400200
C***BEGIN PROLOGUE  DS2Y                                                 SLAP........400300
C***PURPOSE  SLAP Triad to SLAP Column Format Converter.                 SLAP........400400
C            Routine to convert from the SLAP Triad to SLAP Column       SLAP........400500
C            format.                                                     SLAP........400600
C***LIBRARY   SLATEC (SLAP)                                              SLAP........400700
C***CATEGORY  D1B9                                                       SLAP........400800
C***TYPE      DOUBLE PRECISION (SS2Y-S, DS2Y-D)                          SLAP........400900
C***KEYWORDS  LINEAR SYSTEM, SLAP SPARSE                                 SLAP........401000
C***AUTHOR  Seager, Mark K., (LLNL)                                      SLAP........401100
C             Lawrence Livermore National Laboratory                     SLAP........401200
C             PO BOX 808, L-60                                           SLAP........401300
C             Livermore, CA 94550 (510) 423-3141                         SLAP........401400
C             seager@llnl.gov                                            SLAP........401500
C***DESCRIPTION                                                          SLAP........401600
C                                                                        SLAP........401700
C *Usage:                                                                SLAP........401800
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM                          SLAP........401900
C     DOUBLE PRECISION A(NELT)                                           SLAP........402000
C                                                                        SLAP........402100
C     CALL DS2Y( N, NELT, IA, JA, A, ISYM )                              SLAP........402200
C                                                                        SLAP........402300
C *Arguments:                                                            SLAP........402400
C N      :IN       Integer                                               SLAP........402500
C         Order of the Matrix.                                           SLAP........402600
C NELT   :IN       Integer.                                              SLAP........402700
C         Number of non-zeros stored in A.                               SLAP........402800
C IA     :INOUT    Integer IA(NELT).                                     SLAP........402900
C JA     :INOUT    Integer JA(NELT).                                     SLAP........403000
C A      :INOUT    Double Precision A(NELT).                             SLAP........403100
C         These arrays should hold the matrix A in either the SLAP       SLAP........403200
C         Triad format or the SLAP Column format.  See "Description",    SLAP........403300
C         below.  If the SLAP Triad format is used, this format is       SLAP........403400
C         translated to the SLAP Column format by this routine.          SLAP........403500
C ISYM   :IN       Integer.                                              SLAP........403600
C         Flag to indicate symmetric storage format.                     SLAP........403700
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........403800
C         If ISYM=1, the matrix is symmetric, and only the lower         SLAP........403900
C         triangle of the matrix is stored.                              SLAP........404000
C                                                                        SLAP........404100
C *Description:                                                          SLAP........404200
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix     SLAP........404300
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP     SLAP........404400
C       Column format.  The user can hand this routine either of the     SLAP........404500
C       of these data structures.  If the SLAP Triad format is give      SLAP........404600
C       as input then this routine transforms it into SLAP Column        SLAP........404700
C       format.  The way this routine tells which format is given as     SLAP........404800
C       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we        SLAP........404900
C       have the SLAP Column format.  If that equality does not hold     SLAP........405000
C       then it is assumed that the IA, JA, A arrays contain the         SLAP........405100
C       SLAP Triad format.                                               SLAP........405200
C                                                                        SLAP........405300
C       =================== S L A P Triad format ===================     SLAP........405400
C       This routine requires that the  matrix A be   stored in  the     SLAP........405500
C       SLAP  Triad format.  In  this format only the non-zeros  are     SLAP........405600
C       stored.  They may appear in  *ANY* order.  The user supplies     SLAP........405700
C       three arrays of  length NELT, where  NELT is  the number  of     SLAP........405800
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For     SLAP........405900
C       each non-zero the user puts the row and column index of that     SLAP........406000
C       matrix element  in the IA and  JA arrays.  The  value of the     SLAP........406100
C       non-zero   matrix  element is  placed  in  the corresponding     SLAP........406200
C       location of the A array.   This is  an  extremely  easy data     SLAP........406300
C       structure to generate.  On  the  other hand it   is  not too     SLAP........406400
C       efficient on vector computers for  the iterative solution of     SLAP........406500
C       linear systems.  Hence,   SLAP changes   this  input    data     SLAP........406600
C       structure to the SLAP Column format  for  the iteration (but     SLAP........406700
C       does not change it back).                                        SLAP........406800
C                                                                        SLAP........406900
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........407000
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........407100
C                                                                        SLAP........407200
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........407300
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........407400
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........407500
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........407600
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........407700
C       | 0  0  0 44  0|                                                 SLAP........407800
C       |51  0 53  0 55|                                                 SLAP........407900
C                                                                        SLAP........408000
C       =================== S L A P Column format ==================     SLAP........408100
C                                                                        SLAP........408200
C       This routine  requires that  the matrix A  be stored in  the     SLAP........408300
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........408400
C       counting down columns (except for  the diagonal entry, which     SLAP........408500
C       must appear first in each  "column")  and are stored  in the     SLAP........408600
C       double precision array A.   In other words,  for each column     SLAP........408700
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........408800
C       other non-zero  elements going down  the column (except  the     SLAP........408900
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........409000
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........409100
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........409200
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........409300
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........409400
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........409500
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........409600
C       number of columns in  the matrix and NELT  is the number  of     SLAP........409700
C       non-zeros in the matrix.                                         SLAP........409800
C                                                                        SLAP........409900
C       Here is an example of the  SLAP Column  storage format for a     SLAP........410000
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........410100
C       column):                                                         SLAP........410200
C                                                                        SLAP........410300
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........410400
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........410500
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........410600
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........410700
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........410800
C       | 0  0  0 44  0|                                                 SLAP........410900
C       |51  0 53  0 55|                                                 SLAP........411000
C                                                                        SLAP........411100
C***REFERENCES  (NONE)                                                   SLAP........411200
C***ROUTINES CALLED  QS2I1D                                              SLAP........411300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........411400
C   871119  DATE WRITTEN                                                 SLAP........411500
C   881213  Previous REVISION DATE                                       SLAP........411600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........411700
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........411800
C           standard.  (FNF)                                             SLAP........411900
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........412000
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........412100
C   910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF)        SLAP........412200
C   920511  Added complete declaration section.  (WRB)                   SLAP........412300
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........412400
C***END PROLOGUE  DS2Y                                                   SLAP........412500
C     .. Scalar Arguments ..                                             SLAP........412600
      INTEGER ISYM, N, NELT                                              SLAP........412700
C     .. Array Arguments ..                                              SLAP........412800
      DOUBLE PRECISION A(NELT)                                           SLAP........412900
      INTEGER IA(NELT), JA(NELT)                                         SLAP........413000
C     .. Local Scalars ..                                                SLAP........413100
      DOUBLE PRECISION TEMP                                              SLAP........413200
      INTEGER I, IBGN, ICOL, IEND, ITEMP, J                              SLAP........413300
C     .. External Subroutines ..                                         SLAP........413400
      EXTERNAL QS2I1D                                                    SLAP........413500
C***FIRST EXECUTABLE STATEMENT  DS2Y                                     SLAP........413600
C                                                                        SLAP........413700
C         Check to see if the (IA,JA,A) arrays are in SLAP Column        SLAP........413800
C         format.  If it's not then transform from SLAP Triad.           SLAP........413900
C                                                                        SLAP........414000
      IF( JA(N+1).EQ.NELT+1 ) RETURN                                     SLAP........414100
C                                                                        SLAP........414200
C         Sort into ascending order by COLUMN (on the ja array).         SLAP........414300
C         This will line up the columns.                                 SLAP........414400
C                                                                        SLAP........414500
      CALL QS2I1D( JA, IA, A, NELT, 1 )                                  SLAP........414600
C                                                                        SLAP........414700
C         Loop over each column to see where the column indices change   SLAP........414800
C         in the column index array ja.  This marks the beginning of the SLAP........414900
C         next column.                                                   SLAP........415000
C                                                                        SLAP........415100
CVD$R NOVECTOR                                                           SLAP........415200
      JA(1) = 1                                                          SLAP........415300
      DO 20 ICOL = 1, N-1                                                SLAP........415400
         DO 10 J = JA(ICOL)+1, NELT                                      SLAP........415500
            IF( JA(J).NE.ICOL ) THEN                                     SLAP........415600
               JA(ICOL+1) = J                                            SLAP........415700
               GOTO 20                                                   SLAP........415800
            ENDIF                                                        SLAP........415900
 10      CONTINUE                                                        SLAP........416000
 20   CONTINUE                                                           SLAP........416100
      JA(N+1) = NELT+1                                                   SLAP........416200
C                                                                        SLAP........416300
C         Mark the n+2 element so that future calls to a SLAP routine    SLAP........416400
C         utilizing the YSMP-Column storage format will be able to tell. SLAP........416500
C                                                                        SLAP........416600
      JA(N+2) = 0                                                        SLAP........416700
C                                                                        SLAP........416800
C         Now loop through the IA array making sure that the diagonal    SLAP........416900
C         matrix element appears first in the column.  Then sort the     SLAP........417000
C         rest of the column in ascending order.                         SLAP........417100
C                                                                        SLAP........417200
      DO 70 ICOL = 1, N                                                  SLAP........417300
         IBGN = JA(ICOL)                                                 SLAP........417400
         IEND = JA(ICOL+1)-1                                             SLAP........417500
         DO 30 I = IBGN, IEND                                            SLAP........417600
            IF( IA(I).EQ.ICOL ) THEN                                     SLAP........417700
C                                                                        SLAP........417800
C              Swap the diagonal element with the first element in the   SLAP........417900
C              column.                                                   SLAP........418000
C                                                                        SLAP........418100
               ITEMP = IA(I)                                             SLAP........418200
               IA(I) = IA(IBGN)                                          SLAP........418300
               IA(IBGN) = ITEMP                                          SLAP........418400
               TEMP = A(I)                                               SLAP........418500
               A(I) = A(IBGN)                                            SLAP........418600
               A(IBGN) = TEMP                                            SLAP........418700
               GOTO 40                                                   SLAP........418800
            ENDIF                                                        SLAP........418900
 30      CONTINUE                                                        SLAP........419000
 40      IBGN = IBGN + 1                                                 SLAP........419100
         IF( IBGN.LT.IEND ) THEN                                         SLAP........419200
            DO 60 I = IBGN, IEND                                         SLAP........419300
               DO 50 J = I+1, IEND                                       SLAP........419400
                  IF( IA(I).GT.IA(J) ) THEN                              SLAP........419500
                     ITEMP = IA(I)                                       SLAP........419600
                     IA(I) = IA(J)                                       SLAP........419700
                     IA(J) = ITEMP                                       SLAP........419800
                     TEMP = A(I)                                         SLAP........419900
                     A(I) = A(J)                                         SLAP........420000
                     A(J) = TEMP                                         SLAP........420100
                  ENDIF                                                  SLAP........420200
 50            CONTINUE                                                  SLAP........420300
 60         CONTINUE                                                     SLAP........420400
         ENDIF                                                           SLAP........420500
 70   CONTINUE                                                           SLAP........420600
      RETURN                                                             SLAP........420700
C------------- LAST LINE OF DS2Y FOLLOWS ----------------------------    SLAP........420800
      END                                                                SLAP........420900
*DECK DSCAL                                                              SLAP........421000
      SUBROUTINE DSCAL (N, DA, DX, INCX)                                 SLAP........421100
C***BEGIN PROLOGUE  DSCAL                                                SLAP........421200
C***PURPOSE  Multiply a vector by a constant.                            SLAP........421300
C***LIBRARY   SLATEC (BLAS)                                              SLAP........421400
C***CATEGORY  D1A6                                                       SLAP........421500
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)               SLAP........421600
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR                        SLAP........421700
C***AUTHOR  Lawson, C. L., (JPL)                                         SLAP........421800
C           Hanson, R. J., (SNLA)                                        SLAP........421900
C           Kincaid, D. R., (U. of Texas)                                SLAP........422000
C           Krogh, F. T., (JPL)                                          SLAP........422100
C***DESCRIPTION                                                          SLAP........422200
C                                                                        SLAP........422300
C                B L A S  Subprogram                                     SLAP........422400
C    Description of Parameters                                           SLAP........422500
C                                                                        SLAP........422600
C     --Input--                                                          SLAP........422700
C        N  number of elements in input vector(s)                        SLAP........422800
C       DA  double precision scale factor                                SLAP........422900
C       DX  double precision vector with N elements                      SLAP........423000
C     INCX  storage spacing between elements of DX                       SLAP........423100
C                                                                        SLAP........423200
C     --Output--                                                         SLAP........423300
C       DX  double precision result (unchanged if N.LE.0)                SLAP........423400
C                                                                        SLAP........423500
C     Replace double precision DX by double precision DA*DX.             SLAP........423600
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),  SLAP........423700
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.               SLAP........423800
C                                                                        SLAP........423900
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.      SLAP........424000
C                 Krogh, Basic linear algebra subprograms for Fortran    SLAP........424100
C                 usage, Algorithm No. 539, Transactions on Mathematical SLAP........424200
C                 Software 5, 3 (September 1979), pp. 308-323.           SLAP........424300
C***ROUTINES CALLED  (NONE)                                              SLAP........424400
C***REVISION HISTORY  (YYMMDD)                                           SLAP........424500
C   791001  DATE WRITTEN                                                 SLAP........424600
C   890831  Modified array declarations.  (WRB)                          SLAP........424700
C   890831  REVISION DATE from Version 3.2                               SLAP........424800
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........424900
C   900821  Modified to correct problem with a negative increment.       SLAP........425000
C           (WRB)                                                        SLAP........425100
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........425200
C***END PROLOGUE  DSCAL                                                  SLAP........425300
      DOUBLE PRECISION DA, DX(*)                                         SLAP........425400
      INTEGER I, INCX, IX, M, MP1, N                                     SLAP........425500
C***FIRST EXECUTABLE STATEMENT  DSCAL                                    SLAP........425600
      IF (N .LE. 0) RETURN                                               SLAP........425700
      IF (INCX .EQ. 1) GOTO 20                                           SLAP........425800
C                                                                        SLAP........425900
C     Code for increment not equal to 1.                                 SLAP........426000
C                                                                        SLAP........426100
      IX = 1                                                             SLAP........426200
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1                              SLAP........426300
      DO 10 I = 1,N                                                      SLAP........426400
        DX(IX) = DA*DX(IX)                                               SLAP........426500
        IX = IX + INCX                                                   SLAP........426600
   10 CONTINUE                                                           SLAP........426700
      RETURN                                                             SLAP........426800
C                                                                        SLAP........426900
C     Code for increment equal to 1.                                     SLAP........427000
C                                                                        SLAP........427100
C     Clean-up loop so remaining vector length is a multiple of 5.       SLAP........427200
C                                                                        SLAP........427300
   20 M = MOD(N,5)                                                       SLAP........427400
      IF (M .EQ. 0) GOTO 40                                              SLAP........427500
      DO 30 I = 1,M                                                      SLAP........427600
        DX(I) = DA*DX(I)                                                 SLAP........427700
   30 CONTINUE                                                           SLAP........427800
      IF (N .LT. 5) RETURN                                               SLAP........427900
   40 MP1 = M + 1                                                        SLAP........428000
      DO 50 I = MP1,N,5                                                  SLAP........428100
        DX(I) = DA*DX(I)                                                 SLAP........428200
        DX(I+1) = DA*DX(I+1)                                             SLAP........428300
        DX(I+2) = DA*DX(I+2)                                             SLAP........428400
        DX(I+3) = DA*DX(I+3)                                             SLAP........428500
        DX(I+4) = DA*DX(I+4)                                             SLAP........428600
   50 CONTINUE                                                           SLAP........428700
      RETURN                                                             SLAP........428800
      END                                                                SLAP........428900
*DECK DSICCG                                                             SLAP........429000
      SUBROUTINE DSICCG (N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,      SLAP........429100
     +   ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)       SLAP........429200
C***BEGIN PROLOGUE  DSICCG                                               SLAP........429300
C***PURPOSE  Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver.  SLAP........429400
C            Routine to solve a symmetric positive definite linear       SLAP........429500
C            system  Ax = b  using the incomplete Cholesky               SLAP........429600
C            Preconditioned Conjugate Gradient method.                   SLAP........429700
C***LIBRARY   SLATEC (SLAP)                                              SLAP........429800
C***CATEGORY  D2B4                                                       SLAP........429900
C***TYPE      DOUBLE PRECISION (SSICCG-S, DSICCG-D)                      SLAP........430000
C***KEYWORDS  INCOMPLETE CHOLESKY, ITERATIVE PRECONDITION, SLAP, SPARSE, SLAP........430100
C             SYMMETRIC LINEAR SYSTEM                                    SLAP........430200
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........430300
C           Seager, Mark K., (LLNL)                                      SLAP........430400
C             Lawrence Livermore National Laboratory                     SLAP........430500
C             PO BOX 808, L-60                                           SLAP........430600
C             Livermore, CA 94550 (510) 423-3141                         SLAP........430700
C             seager@llnl.gov                                            SLAP........430800
C***DESCRIPTION                                                          SLAP........430900
C                                                                        SLAP........431000
C *Usage:                                                                SLAP........431100
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX             SLAP........431200
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NL+2*N+1), LENIW            SLAP........431300
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(NL+5*N)      SLAP........431400
C                                                                        SLAP........431500
C     CALL DSICCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,             SLAP........431600
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )    SLAP........431700
C                                                                        SLAP........431800
C *Arguments:                                                            SLAP........431900
C N      :IN       Integer.                                              SLAP........432000
C         Order of the Matrix.                                           SLAP........432100
C B      :IN       Double Precision B(N).                                SLAP........432200
C         Right-hand side vector.                                        SLAP........432300
C X      :INOUT    Double Precision X(N).                                SLAP........432400
C         On input X is your initial guess for solution vector.          SLAP........432500
C         On output X is the final approximate solution.                 SLAP........432600
C NELT   :IN       Integer.                                              SLAP........432700
C         Number of Non-Zeros stored in A.                               SLAP........432800
C IA     :INOUT    Integer IA(NELT).                                     SLAP........432900
C JA     :INOUT    Integer JA(NELT).                                     SLAP........433000
C A      :INOUT    Double Precision A(NELT).                             SLAP........433100
C         These arrays should hold the matrix A in either the SLAP       SLAP........433200
C         Triad format or the SLAP Column format.  See "Description",    SLAP........433300
C         below.  If the SLAP Triad format is chosen it is changed       SLAP........433400
C         internally to the SLAP Column format.                          SLAP........433500
C ISYM   :IN       Integer.                                              SLAP........433600
C         Flag to indicate symmetric storage format.                     SLAP........433700
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........433800
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........433900
C         or lower triangle of the matrix is stored.                     SLAP........434000
C ITOL   :IN       Integer.                                              SLAP........434100
C         Flag to indicate type of convergence criterion.                SLAP........434200
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........434300
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........434400
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........434500
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........434600
C         side is less than TOL, where M-inv is the inverse of the       SLAP........434700
C         diagonal of A.                                                 SLAP........434800
C         ITOL=11 is often useful for checking and comparing different   SLAP........434900
C         routines.  For this case, the user must supply the "exact"     SLAP........435000
C         solution or a very accurate approximation (one with an error   SLAP........435100
C         much less than TOL) through a common block,                    SLAP........435200
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........435300
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........435400
C         between the iterative approximation and the user-supplied      SLAP........435500
C         solution divided by the 2-norm of the user-supplied solution   SLAP........435600
C         is less than TOL.  Note that this requires the user to set up  SLAP........435700
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........435800
C         The routine with this declaration should be loaded before the  SLAP........435900
C         stop test so that the correct length is used by the loader.    SLAP........436000
C         This procedure is not standard Fortran and may not work        SLAP........436100
C         correctly on your system (although it has worked on every      SLAP........436200
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........436300
C         common block is indeed standard Fortran.                       SLAP........436400
C TOL    :INOUT    Double Precision.                                     SLAP........436500
C         Convergence criterion, as described above.  (Reset if IERR=4.) SLAP........436600
C ITMAX  :IN       Integer.                                              SLAP........436700
C         Maximum number of iterations.                                  SLAP........436800
C ITER   :OUT      Integer.                                              SLAP........436900
C         Number of iterations required to reach convergence, or         SLAP........437000
C         ITMAX+1 if convergence criterion could not be achieved in      SLAP........437100
C         ITMAX iterations.                                              SLAP........437200
C ERR    :OUT      Double Precision.                                     SLAP........437300
C         Error estimate of error in final approximate solution, as      SLAP........437400
C         defined by ITOL.                                               SLAP........437500
C IERR   :OUT      Integer.                                              SLAP........437600
C         Return error flag.                                             SLAP........437700
C           IERR = 0 => All went well.                                   SLAP........437800
C           IERR = 1 => Insufficient space allocated for WORK or IWORK.  SLAP........437900
C           IERR = 2 => Method failed to converge in ITMAX steps.        SLAP........438000
C           IERR = 3 => Error in user input.                             SLAP........438100
C                       Check input values of N, ITOL.                   SLAP........438200
C           IERR = 4 => User error tolerance set too tight.              SLAP........438300
C                       Reset to 500*D1MACH(3).  Iteration proceeded.    SLAP........438400
C           IERR = 5 => Preconditioning matrix, M, is not positive       SLAP........438500
C                       definite.  (r,z) < 0.                            SLAP........438600
C           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.  SLAP........438700
C           IERR = 7 => Incomplete factorization broke down and was      SLAP........438800
C                       fudged.  Resulting preconditioning may be less   SLAP........438900
C                       than the best.                                   SLAP........439000
C IUNIT  :IN       Integer.                                              SLAP........439100
C         Unit number on which to write the error at each iteration,     SLAP........439200
C         if this is desired for monitoring convergence.  If unit        SLAP........439300
C         number is 0, no writing will occur.                            SLAP........439400
C RWORK  :WORK     Double Precision RWORK(LENW).                         SLAP........439500
C         Double Precision array used for workspace.                     SLAP........439600
C LENW   :IN       Integer.                                              SLAP........439700
C         Length of the double precision workspace, RWORK.               SLAP........439800
C         LENW >= NL+5*N.                                                SLAP........439900
C         NL is the number of non-zeros in the lower triangle of the     SLAP........440000
C         matrix (including the diagonal).                               SLAP........440100
C IWORK  :WORK     Integer IWORK(LENIW).                                 SLAP........440200
C         Integer array used for workspace.                              SLAP........440300
C         Upon return the following locations of IWORK hold information  SLAP........440400
C         which may be of use to the user:                               SLAP........440500
C         IWORK(9)  Amount of Integer workspace actually used.           SLAP........440600
C         IWORK(10) Amount of Double Precision workspace actually used.  SLAP........440700
C LENIW  :IN       Integer.                                              SLAP........440800
C         Length of the integer workspace, IWORK.  LENIW >= NL+N+11.     SLAP........440900
C         NL is the number of non-zeros in the lower triangle of the     SLAP........441000
C         matrix (including the diagonal).                               SLAP........441100
C                                                                        SLAP........441200
C *Description:                                                          SLAP........441300
C       This routine  performs  preconditioned  conjugate   gradient     SLAP........441400
C       method on the   symmetric positive  definite  linear  system     SLAP........441500
C       Ax=b.   The preconditioner  is  the incomplete Cholesky (IC)     SLAP........441600
C       factorization of the matrix A.  See  DSICS for details about     SLAP........441700
C       the incomplete   factorization algorithm.  One   should note     SLAP........441800
C       here however, that the  IC factorization is a  slow  process     SLAP........441900
C       and  that  one should   save  factorizations  for  reuse, if     SLAP........442000
C       possible.  The   MSOLVE operation (handled  in  DSLLTI) does     SLAP........442100
C       vectorize on machines  with  hardware  gather/scatter and is     SLAP........442200
C       quite fast.                                                      SLAP........442300
C                                                                        SLAP........442400
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix     SLAP........442500
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP     SLAP........442600
C       Column format.  The user can hand this routine either of the     SLAP........442700
C       of these data structures and SLAP  will figure out  which on     SLAP........442800
C       is being used and act accordingly.                               SLAP........442900
C                                                                        SLAP........443000
C       =================== S L A P Triad format ===================     SLAP........443100
C                                                                        SLAP........443200
C       This routine requires that the  matrix A be   stored in  the     SLAP........443300
C       SLAP  Triad format.  In  this format only the non-zeros  are     SLAP........443400
C       stored.  They may appear in  *ANY* order.  The user supplies     SLAP........443500
C       three arrays of  length NELT, where  NELT is  the number  of     SLAP........443600
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For     SLAP........443700
C       each non-zero the user puts the row and column index of that     SLAP........443800
C       matrix element  in the IA and  JA arrays.  The  value of the     SLAP........443900
C       non-zero   matrix  element is  placed  in  the corresponding     SLAP........444000
C       location of the A array.   This is  an  extremely  easy data     SLAP........444100
C       structure to generate.  On  the  other hand it   is  not too     SLAP........444200
C       efficient on vector computers for  the iterative solution of     SLAP........444300
C       linear systems.  Hence,   SLAP changes   this  input    data     SLAP........444400
C       structure to the SLAP Column format  for  the iteration (but     SLAP........444500
C       does not change it back).                                        SLAP........444600
C                                                                        SLAP........444700
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........444800
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........444900
C                                                                        SLAP........445000
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........445100
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........445200
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........445300
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........445400
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........445500
C       | 0  0  0 44  0|                                                 SLAP........445600
C       |51  0 53  0 55|                                                 SLAP........445700
C                                                                        SLAP........445800
C       =================== S L A P Column format ==================     SLAP........445900
C                                                                        SLAP........446000
C       This routine  requires that  the matrix A  be stored in  the     SLAP........446100
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........446200
C       counting down columns (except for  the diagonal entry, which     SLAP........446300
C       must appear first in each  "column")  and are stored  in the     SLAP........446400
C       double precision array A.   In other words,  for each column     SLAP........446500
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........446600
C       other non-zero  elements going down  the column (except  the     SLAP........446700
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........446800
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........446900
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........447000
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........447100
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........447200
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........447300
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........447400
C       number of columns in  the matrix and NELT  is the number  of     SLAP........447500
C       non-zeros in the matrix.                                         SLAP........447600
C                                                                        SLAP........447700
C       Here is an example of the  SLAP Column  storage format for a     SLAP........447800
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........447900
C       column):                                                         SLAP........448000
C                                                                        SLAP........448100
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........448200
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........448300
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........448400
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........448500
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........448600
C       | 0  0  0 44  0|                                                 SLAP........448700
C       |51  0 53  0 55|                                                 SLAP........448800
C                                                                        SLAP........448900
C *Side Effects:                                                         SLAP........449000
C       The SLAP Triad format (IA, JA, A) is modified internally to be   SLAP........449100
C       the SLAP Column format.  See above.                              SLAP........449200
C                                                                        SLAP........449300
C *Cautions:                                                             SLAP........449400
C     This routine will attempt to write to the Fortran logical output   SLAP........449500
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........449600
C     this logical unit is attached to a file or terminal before calling SLAP........449700
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........449800
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........449900
C                                                                        SLAP........450000
C***SEE ALSO  DCG, DSLLTI                                                SLAP........450100
C***REFERENCES  1. Louis Hageman and David Young, Applied Iterative      SLAP........450200
C                  Methods, Academic Press, New York, 1981.              SLAP........450300
C               2. Concus, Golub and O'Leary, A Generalized Conjugate    SLAP........450400
C                  Gradient Method for the Numerical Solution of         SLAP........450500
C                  Elliptic Partial Differential Equations, in Sparse    SLAP........450600
C                  Matrix Computations, Bunch and Rose, Eds., Academic   SLAP........450700
C                  Press, New York, 1979.                                SLAP........450800
C***ROUTINES CALLED  DCG, DCHKW, DS2Y, DSICS, DSLLTI, DSMV, XERMSG       SLAP........450900
C***REVISION HISTORY  (YYMMDD)                                           SLAP........451000
C   890404  DATE WRITTEN                                                 SLAP........451100
C   890404  Previous REVISION DATE                                       SLAP........451200
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........451300
C   890921  Removed TeX from comments.  (FNF)                            SLAP........451400
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........451500
C           standard.  (FNF)                                             SLAP........451600
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........451700
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)              SLAP........451800
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........451900
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........452000
C   920511  Added complete declaration section.  (WRB)                   SLAP........452100
C   920929  Corrected format of references.  (FNF)                       SLAP........452200
C   921019  Corrected NEL to NL.  (FNF)                                  SLAP........452300
C***END PROLOGUE  DSICCG                                                 SLAP........452400
C     .. Parameters ..                                                   SLAP........452500
      INTEGER LOCRB, LOCIB                                               SLAP........452600
      PARAMETER (LOCRB=1, LOCIB=11)                                      SLAP........452700
C     .. Scalar Arguments ..                                             SLAP........452800
      DOUBLE PRECISION ERR, TOL                                          SLAP........452900
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N, NELT SLAP........453000
C     .. Array Arguments ..                                              SLAP........453100
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N)                  SLAP........453200
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)                           SLAP........453300
C     .. Local Scalars ..                                                SLAP........453400
      INTEGER LOCDIN, LOCDZ, LOCEL, LOCIEL, LOCIW, LOCJEL, LOCP, LOCR,   SLAP........453500
     +        LOCW, LOCZ, NL                                             SLAP........453600
      CHARACTER XERN1*8                                                  SLAP........453700
C     .. External Subroutines ..                                         SLAP........453800
      EXTERNAL DCG, DCHKW, DS2Y, DSICS, DSLLTI, DSMV, XERMSG             SLAP........453900
C***FIRST EXECUTABLE STATEMENT  DSICCG                                   SLAP........454000
C                                                                        SLAP........454100
      IERR = 0                                                           SLAP........454200
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN                                   SLAP........454300
         IERR = 3                                                        SLAP........454400
         RETURN                                                          SLAP........454500
      ENDIF                                                              SLAP........454600
C                                                                        SLAP........454700
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.  SLAP........454800
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )                              SLAP........454900
C                                                                        SLAP........455000
C         Count number of elements in lower triangle of the matrix.      SLAP........455100
C         Then set up the work arrays.                                   SLAP........455200
      IF( ISYM.EQ.0 ) THEN                                               SLAP........455300
         NL = (NELT + N)/2                                               SLAP........455400
      ELSE                                                               SLAP........455500
         NL = NELT                                                       SLAP........455600
      ENDIF                                                              SLAP........455700
C                                                                        SLAP........455800
      LOCJEL = LOCIB                                                     SLAP........455900
      LOCIEL = LOCJEL + NL                                               SLAP........456000
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........456100
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........456200
C        LOCIW = LOCIEL + N + 1                                          SLAP........456300
      LOCIW = LOCIEL + NL                                                SLAP........456400
C                                                                        SLAP........456500
      LOCEL = LOCRB                                                      SLAP........456600
      LOCDIN = LOCEL + NL                                                SLAP........456700
      LOCR = LOCDIN + N                                                  SLAP........456800
      LOCZ = LOCR + N                                                    SLAP........456900
      LOCP = LOCZ + N                                                    SLAP........457000
      LOCDZ = LOCP + N                                                   SLAP........457100
      LOCW = LOCDZ + N                                                   SLAP........457200
C                                                                        SLAP........457300
C         Check the workspace allocations.                               SLAP........457400
      CALL DCHKW( 'DSICCG', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )  SLAP........457500
      IF( IERR.NE.0 ) RETURN                                             SLAP........457600
C                                                                        SLAP........457700
      IWORK(1) = NL                                                      SLAP........457800
      IWORK(2) = LOCJEL                                                  SLAP........457900
      IWORK(3) = LOCIEL                                                  SLAP........458000
      IWORK(4) = LOCEL                                                   SLAP........458100
      IWORK(5) = LOCDIN                                                  SLAP........458200
      IWORK(9) = LOCIW                                                   SLAP........458300
      IWORK(10) = LOCW                                                   SLAP........458400
C                                                                        SLAP........458500
C         Compute the Incomplete Cholesky decomposition.                 SLAP........458600
C                                                                        SLAP........458700
      CALL DSICS(N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIEL),            SLAP........458800
     $     IWORK(LOCJEL), RWORK(LOCEL), RWORK(LOCDIN),                   SLAP........458900
     $     RWORK(LOCR), IERR )                                           SLAP........459000
      IF( IERR.NE.0 ) THEN                                               SLAP........459100
         WRITE (XERN1, '(I8)') IERR                                      SLAP........459200
         CALL XERMSG ('SLATEC', 'DSICCG',                                SLAP........459300
     $      'IC factorization broke down on step ' // XERN1 //           SLAP........459400
     $      '.  Diagonal was set to unity and factorization proceeded.', SLAP........459500
     $      1, 1)                                                        SLAP........459600
         IERR = 7                                                        SLAP........459700
      ENDIF                                                              SLAP........459800
C                                                                        SLAP........459900
C         Do the Preconditioned Conjugate Gradient.                      SLAP........460000
      CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLLTI,             SLAP........460100
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR),        SLAP........460200
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCDZ), RWORK(1),             SLAP........460300
     $     IWORK(1))                                                     SLAP........460400
      RETURN                                                             SLAP........460500
C------------- LAST LINE OF DSICCG FOLLOWS ----------------------------  SLAP........460600
      END                                                                SLAP........460700
*DECK DSICS                                                              SLAP........460800
      SUBROUTINE DSICS (N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D,  SLAP........460900
     +   R, IWARN)                                                       SLAP........461000
C***BEGIN PROLOGUE  DSICS                                                SLAP........461100
C***PURPOSE  Incompl. Cholesky Decomposition Preconditioner SLAP Set Up. SLAP........461200
C            Routine to generate the Incomplete Cholesky decomposition,  SLAP........461300
C            L*D*L-trans, of a symmetric positive definite matrix, A,    SLAP........461400
C            which is stored in SLAP Column format.  The unit lower      SLAP........461500
C            triangular matrix L is stored by rows, and the inverse of   SLAP........461600
C            the diagonal matrix D is stored.                            SLAP........461700
C***LIBRARY   SLATEC (SLAP)                                              SLAP........461800
C***CATEGORY  D2E                                                        SLAP........461900
C***TYPE      DOUBLE PRECISION (SSICS-S, DSICS-D)                        SLAP........462000
C***KEYWORDS  INCOMPLETE CHOLESKY FACTORIZATION,                         SLAP........462100
C             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE         SLAP........462200
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........462300
C           Seager, Mark K., (LLNL)                                      SLAP........462400
C             Lawrence Livermore National Laboratory                     SLAP........462500
C             PO BOX 808, L-60                                           SLAP........462600
C             Livermore, CA 94550 (510) 423-3141                         SLAP........462700
C             seager@llnl.gov                                            SLAP........462800
C***DESCRIPTION                                                          SLAP........462900
C                                                                        SLAP........463000
C *Usage:                                                                SLAP........463100
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM                          SLAP........463200
C     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN                             SLAP........463300
C     DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)                      SLAP........463400
C                                                                        SLAP........463500
C     CALL DSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,     SLAP........463600
C    $    IWARN )                                                        SLAP........463700
C                                                                        SLAP........463800
C *Arguments:                                                            SLAP........463900
C N      :IN       Integer.                                              SLAP........464000
C         Order of the Matrix.                                           SLAP........464100
C NELT   :IN       Integer.                                              SLAP........464200
C         Number of elements in arrays IA, JA, and A.                    SLAP........464300
C IA     :INOUT    Integer IA(NELT).                                     SLAP........464400
C JA     :INOUT    Integer JA(NELT).                                     SLAP........464500
C A      :INOUT    Double Precision A(NELT).                             SLAP........464600
C         These arrays should hold the matrix A in the SLAP Column       SLAP........464700
C         format.  See "Description", below.                             SLAP........464800
C ISYM   :IN       Integer.                                              SLAP........464900
C         Flag to indicate symmetric storage format.                     SLAP........465000
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........465100
C         If ISYM=1, the matrix is symmetric, and only the lower         SLAP........465200
C         triangle of the matrix is stored.                              SLAP........465300
C NEL    :OUT      Integer.                                              SLAP........465400
C         Number of non-zeros in the lower triangle of A.   Also         SLAP........465500
C         corresponds to the length of the IEL, JEL, EL arrays.          SLAP........465600
C IEL    :OUT      Integer IEL(NEL).                                     SLAP........465700
C JEL    :OUT      Integer JEL(NEL).                                     SLAP........465800
C EL     :OUT      Double Precision EL(NEL).                             SLAP........465900
C         IEL, JEL, EL contain the unit lower triangular factor  of the  SLAP........466000
C         incomplete decomposition   of the A  matrix  stored  in  SLAP  SLAP........466100
C         Row format.   The Diagonal of   ones   *IS*   stored.     See  SLAP........466200
C         "Description", below for more details about the SLAP Row fmt.  SLAP........466300
C D      :OUT      Double Precision D(N)                                 SLAP........466400
C         Upon return this array holds D(I) = 1./DIAG(A).                SLAP........466500
C R      :WORK     Double Precision R(N).                                SLAP........466600
C         Temporary double precision workspace needed for the            SLAP........466700
C         factorization.                                                 SLAP........466800
C IWARN  :OUT      Integer.                                              SLAP........466900
C         This is a warning variable and is zero if the IC factoriza-    SLAP........467000
C         tion goes well.  It is set to the row index corresponding to   SLAP........467100
C         the last zero pivot found.  See "Description", below.          SLAP........467200
C                                                                        SLAP........467300
C *Description                                                           SLAP........467400
C       =================== S L A P Column format ==================     SLAP........467500
C       This routine  requires that  the matrix A  be stored in  the     SLAP........467600
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........467700
C       counting down columns (except for  the diagonal entry, which     SLAP........467800
C       must appear first in each  "column")  and are stored  in the     SLAP........467900
C       double precision array A.   In other words,  for each column     SLAP........468000
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........468100
C       other non-zero  elements going down  the column (except  the     SLAP........468200
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........468300
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........468400
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........468500
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........468600
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........468700
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........468800
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........468900
C       number of columns in  the matrix and NELT  is the number  of     SLAP........469000
C       non-zeros in the matrix.                                         SLAP........469100
C                                                                        SLAP........469200
C       Here is an example of the  SLAP Column  storage format for a     SLAP........469300
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........469400
C       column):                                                         SLAP........469500
C                                                                        SLAP........469600
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........469700
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........469800
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........469900
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........470000
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........470100
C       | 0  0  0 44  0|                                                 SLAP........470200
C       |51  0 53  0 55|                                                 SLAP........470300
C                                                                        SLAP........470400
C       ==================== S L A P Row format ====================     SLAP........470500
C                                                                        SLAP........470600
C       This routine requires  that the matrix A  be  stored  in the     SLAP........470700
C       SLAP  Row format.   In this format  the non-zeros are stored     SLAP........470800
C       counting across  rows (except for the diagonal  entry, which     SLAP........470900
C       must  appear first  in each  "row")  and  are stored  in the     SLAP........471000
C       double precision  array A.  In other words, for each row  in     SLAP........471100
C       the matrix  put the diagonal  entry in A.   Then put in  the     SLAP........471200
C       other  non-zero elements  going across  the row  (except the     SLAP........471300
C       diagonal) in order.  The JA array holds the column index for     SLAP........471400
C       each non-zero.  The IA array holds the offsets  into the JA,     SLAP........471500
C       A  arrays  for  the   beginning  of  each  row.    That  is,     SLAP........471600
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-     SLAP........471700
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)     SLAP........471800
C       are  the last elements  of the  IROW-th row.   Note  that we     SLAP........471900
C       always have  IA(N+1) = NELT+1, where N is the number of rows     SLAP........472000
C       in the matrix  and  NELT is the  number of non-zeros  in the     SLAP........472100
C       matrix.                                                          SLAP........472200
C                                                                        SLAP........472300
C       Here is an example of the SLAP Row storage format for a  5x5     SLAP........472400
C       Matrix (in the A and JA arrays '|' denotes the end of a row):    SLAP........472500
C                                                                        SLAP........472600
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.   SLAP........472700
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........472800
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53   SLAP........472900
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........473000
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12                      SLAP........473100
C       | 0  0  0 44  0|                                                 SLAP........473200
C       |51  0 53  0 55|                                                 SLAP........473300
C                                                                        SLAP........473400
C       With the SLAP  format some  of  the   "inner  loops" of this     SLAP........473500
C       routine should vectorize  on  machines with hardware support     SLAP........473600
C       for vector   gather/scatter  operations.  Your compiler  may     SLAP........473700
C       require a compiler directive to  convince it that  there are     SLAP........473800
C       no  implicit  vector  dependencies.  Compiler directives for     SLAP........473900
C       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are     SLAP........474000
C       supplied with the standard SLAP distribution.                    SLAP........474100
C                                                                        SLAP........474200
C       The IC factorization does not always exist for SPD matrices.     SLAP........474300
C       In the event that a zero pivot is found it is set  to be 1.0     SLAP........474400
C       and the factorization proceeds.   The integer variable IWARN     SLAP........474500
C       is set to the last row where the Diagonal was fudged.  This      SLAP........474600
C       eventuality hardly ever occurs in practice.                      SLAP........474700
C                                                                        SLAP........474800
C***SEE ALSO  DCG, DSICCG                                                SLAP........474900
C***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations, SLAP........475000
C                  Johns Hopkins University Press, Baltimore, Maryland,  SLAP........475100
C                  1983.                                                 SLAP........475200
C***ROUTINES CALLED  XERMSG                                              SLAP........475300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........475400
C   890404  DATE WRITTEN                                                 SLAP........475500
C   890404  Previous REVISION DATE                                       SLAP........475600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........475700
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........475800
C           standard.  (FNF)                                             SLAP........475900
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........476000
C   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)              SLAP........476100
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........476200
C   920511  Added complete declaration section.  (WRB)                   SLAP........476300
C   920929  Corrected format of reference.  (FNF)                        SLAP........476400
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........476500
C***END PROLOGUE  DSICS                                                  SLAP........476600
C     .. Scalar Arguments ..                                             SLAP........476700
      INTEGER ISYM, IWARN, N, NEL, NELT                                  SLAP........476800
C     .. Array Arguments ..                                              SLAP........476900
      DOUBLE PRECISION A(NELT), D(N), EL(NEL), R(N)                      SLAP........477000
      INTEGER IA(NELT), IEL(NEL), JA(NELT), JEL(NEL)                     SLAP........477100
C     .. Local Scalars ..                                                SLAP........477200
      DOUBLE PRECISION ELTMP                                             SLAP........477300
      INTEGER I, IBGN, IC, ICBGN, ICEND, ICOL, IEND, IR, IRBGN, IREND,   SLAP........477400
     +        IROW, IRR, J, JBGN, JELTMP, JEND                           SLAP........477500
      CHARACTER XERN1*8                                                  SLAP........477600
C     .. External Subroutines ..                                         SLAP........477700
      EXTERNAL XERMSG                                                    SLAP........477800
C***FIRST EXECUTABLE STATEMENT  DSICS                                    SLAP........477900
C                                                                        SLAP........478000
C         Set the lower triangle in IEL, JEL, EL                         SLAP........478100
C                                                                        SLAP........478200
      IWARN = 0                                                          SLAP........478300
C                                                                        SLAP........478400
C         All matrix elements stored in IA, JA, A.  Pick out the lower   SLAP........478500
C         triangle (making sure that the Diagonal of EL is one) and      SLAP........478600
C         store by rows.                                                 SLAP........478700
C                                                                        SLAP........478800
      NEL = 1                                                            SLAP........478900
      IEL(1) = 1                                                         SLAP........479000
      JEL(1) = 1                                                         SLAP........479100
      EL(1) = 1                                                          SLAP........479200
      D(1) = A(1)                                                        SLAP........479300
CVD$R NOCONCUR                                                           SLAP........479400
      DO 30 IROW = 2, N                                                  SLAP........479500
C         Put in the Diagonal.                                           SLAP........479600
         NEL = NEL + 1                                                   SLAP........479700
         IEL(IROW) = NEL                                                 SLAP........479800
         JEL(NEL) = IROW                                                 SLAP........479900
         EL(NEL) = 1                                                     SLAP........480000
         D(IROW) = A(JA(IROW))                                           SLAP........480100
C                                                                        SLAP........480200
C         Look in all the lower triangle columns for a matching row.     SLAP........480300
C         Since the matrix is symmetric, we can look across the          SLAP........480400
C         IROW-th row by looking down the IROW-th column (if it is       SLAP........480500
C         stored ISYM=0)...                                              SLAP........480600
         IF( ISYM.EQ.0 ) THEN                                            SLAP........480700
            ICBGN = JA(IROW)                                             SLAP........480800
            ICEND = JA(IROW+1)-1                                         SLAP........480900
         ELSE                                                            SLAP........481000
            ICBGN = 1                                                    SLAP........481100
            ICEND = IROW-1                                               SLAP........481200
         ENDIF                                                           SLAP........481300
         DO 20 IC = ICBGN, ICEND                                         SLAP........481400
            IF( ISYM.EQ.0 ) THEN                                         SLAP........481500
               ICOL = IA(IC)                                             SLAP........481600
               IF( ICOL.GE.IROW ) GOTO 20                                SLAP........481700
            ELSE                                                         SLAP........481800
               ICOL = IC                                                 SLAP........481900
            ENDIF                                                        SLAP........482000
            JBGN = JA(ICOL)+1                                            SLAP........482100
            JEND = JA(ICOL+1)-1                                          SLAP........482200
            IF( JBGN.LE.JEND .AND. IA(JEND).GE.IROW ) THEN               SLAP........482300
CVD$ NOVECTOR                                                            SLAP........482400
               DO 10 J = JBGN, JEND                                      SLAP........482500
                  IF( IA(J).EQ.IROW ) THEN                               SLAP........482600
                     NEL = NEL + 1                                       SLAP........482700
                     JEL(NEL) = ICOL                                     SLAP........482800
                     EL(NEL)  = A(J)                                     SLAP........482900
                     GOTO 20                                             SLAP........483000
                  ENDIF                                                  SLAP........483100
 10            CONTINUE                                                  SLAP........483200
            ENDIF                                                        SLAP........483300
 20      CONTINUE                                                        SLAP........483400
 30   CONTINUE                                                           SLAP........483500
      IEL(N+1) = NEL+1                                                   SLAP........483600
C                                                                        SLAP........483700
C         Sort ROWS of lower triangle into descending order (count out   SLAP........483800
C         along rows out from Diagonal).                                 SLAP........483900
C                                                                        SLAP........484000
      DO 60 IROW = 2, N                                                  SLAP........484100
         IBGN = IEL(IROW)+1                                              SLAP........484200
         IEND = IEL(IROW+1)-1                                            SLAP........484300
         IF( IBGN.LT.IEND ) THEN                                         SLAP........484400
            DO 50 I = IBGN, IEND-1                                       SLAP........484500
CVD$ NOVECTOR                                                            SLAP........484600
               DO 40 J = I+1, IEND                                       SLAP........484700
                  IF( JEL(I).GT.JEL(J) ) THEN                            SLAP........484800
                     JELTMP = JEL(J)                                     SLAP........484900
                     JEL(J) = JEL(I)                                     SLAP........485000
                     JEL(I) = JELTMP                                     SLAP........485100
                     ELTMP = EL(J)                                       SLAP........485200
                     EL(J) = EL(I)                                       SLAP........485300
                     EL(I) = ELTMP                                       SLAP........485400
                  ENDIF                                                  SLAP........485500
 40            CONTINUE                                                  SLAP........485600
 50         CONTINUE                                                     SLAP........485700
         ENDIF                                                           SLAP........485800
 60   CONTINUE                                                           SLAP........485900
C                                                                        SLAP........486000
C         Perform the Incomplete Cholesky decomposition by looping       SLAP........486100
C         over the rows.                                                 SLAP........486200
C         Scale the first column.  Use the structure of A to pick out    SLAP........486300
C         the rows with something in column 1.                           SLAP........486400
C                                                                        SLAP........486500
      IRBGN = JA(1)+1                                                    SLAP........486600
      IREND = JA(2)-1                                                    SLAP........486700
      DO 65 IRR = IRBGN, IREND                                           SLAP........486800
         IR = IA(IRR)                                                    SLAP........486900
C         Find the index into EL for EL(1,IR).                           SLAP........487000
C         Hint: it's the second entry.                                   SLAP........487100
         I = IEL(IR)+1                                                   SLAP........487200
         EL(I) = EL(I)/D(1)                                              SLAP........487300
 65   CONTINUE                                                           SLAP........487400
C                                                                        SLAP........487500
      DO 110 IROW = 2, N                                                 SLAP........487600
C                                                                        SLAP........487700
C         Update the IROW-th diagonal.                                   SLAP........487800
C                                                                        SLAP........487900
         DO 66 I = 1, IROW-1                                             SLAP........488000
            R(I) = 0                                                     SLAP........488100
 66      CONTINUE                                                        SLAP........488200
         IBGN = IEL(IROW)+1                                              SLAP........488300
         IEND = IEL(IROW+1)-1                                            SLAP........488400
         IF( IBGN.LE.IEND ) THEN                                         SLAP........488500
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........488600
CDIR$ IVDEP                                                              SLAP........488700
CVD$ NODEPCHK                                                            SLAP........488800
            DO 70 I = IBGN, IEND                                         SLAP........488900
               R(JEL(I)) = EL(I)*D(JEL(I))                               SLAP........489000
               D(IROW) = D(IROW) - EL(I)*R(JEL(I))                       SLAP........489100
 70         CONTINUE                                                     SLAP........489200
C                                                                        SLAP........489300
C         Check to see if we have a problem with the diagonal.           SLAP........489400
C                                                                        SLAP........489500
            IF( D(IROW).LE.0.0D0 ) THEN                                  SLAP........489600
               IF( IWARN.EQ.0 ) IWARN = IROW                             SLAP........489700
               D(IROW) = 1                                               SLAP........489800
            ENDIF                                                        SLAP........489900
         ENDIF                                                           SLAP........490000
C                                                                        SLAP........490100
C         Update each EL(IROW+1:N,IROW), if there are any.               SLAP........490200
C         Use the structure of A to determine the Non-zero elements      SLAP........490300
C         of the IROW-th column of EL.                                   SLAP........490400
C                                                                        SLAP........490500
         IRBGN = JA(IROW)                                                SLAP........490600
         IREND = JA(IROW+1)-1                                            SLAP........490700
         DO 100 IRR = IRBGN, IREND                                       SLAP........490800
            IR = IA(IRR)                                                 SLAP........490900
            IF( IR.LE.IROW ) GOTO 100                                    SLAP........491000
C         Find the index into EL for EL(IR,IROW)                         SLAP........491100
            IBGN = IEL(IR)+1                                             SLAP........491200
            IEND = IEL(IR+1)-1                                           SLAP........491300
            IF( JEL(IBGN).GT.IROW ) GOTO 100                             SLAP........491400
            DO 90 I = IBGN, IEND                                         SLAP........491500
               IF( JEL(I).EQ.IROW ) THEN                                 SLAP........491600
                  ICEND = IEND                                           SLAP........491700
 91               IF( JEL(ICEND).GE.IROW ) THEN                          SLAP........491800
                     ICEND = ICEND - 1                                   SLAP........491900
                     GOTO 91                                             SLAP........492000
                  ENDIF                                                  SLAP........492100
C         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.          SLAP........492200
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........492300
CDIR$ IVDEP                                                              SLAP........492400
CVD$ NODEPCHK                                                            SLAP........492500
                  DO 80 IC = IBGN, ICEND                                 SLAP........492600
                     EL(I) = EL(I) - EL(IC)*R(JEL(IC))                   SLAP........492700
 80               CONTINUE                                               SLAP........492800
                  EL(I) = EL(I)/D(IROW)                                  SLAP........492900
                  GOTO 100                                               SLAP........493000
               ENDIF                                                     SLAP........493100
 90         CONTINUE                                                     SLAP........493200
C                                                                        SLAP........493300
C         If we get here, we have real problems...                       SLAP........493400
            WRITE (XERN1, '(I8)') IROW                                   SLAP........493500
            CALL XERMSG ('SLATEC', 'DSICS',                              SLAP........493600
     $         'A and EL data structure mismatch in row '// XERN1, 1, 2) SLAP........493700
 100     CONTINUE                                                        SLAP........493800
 110  CONTINUE                                                           SLAP........493900
C                                                                        SLAP........494000
C         Replace diagonals by their inverses.                           SLAP........494100
C                                                                        SLAP........494200
CVD$ CONCUR                                                              SLAP........494300
      DO 120 I =1, N                                                     SLAP........494400
         D(I) = 1.0D0/D(I)                                               SLAP........494500
 120  CONTINUE                                                           SLAP........494600
      RETURN                                                             SLAP........494700
C------------- LAST LINE OF DSICS FOLLOWS ----------------------------   SLAP........494800
      END                                                                SLAP........494900
*DECK DSILUS                                                             SLAP........495000
      SUBROUTINE DSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV,  SLAP........495100
     +   NU, IU, JU, U, NROW, NCOL)                                      SLAP........495200
C***BEGIN PROLOGUE  DSILUS                                               SLAP........495300
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.     SLAP........495400
C            Routine to generate the incomplete LDU decomposition of a   SLAP........495500
C            matrix.  The unit lower triangular factor L is stored by    SLAP........495600
C            rows and the unit upper triangular factor U is stored by    SLAP........495700
C            columns.  The inverse of the diagonal matrix D is stored.   SLAP........495800
C            No fill in is allowed.                                      SLAP........495900
C***LIBRARY   SLATEC (SLAP)                                              SLAP........496000
C***CATEGORY  D2E                                                        SLAP........496100
C***TYPE      DOUBLE PRECISION (SSILUS-S, DSILUS-D)                      SLAP........496200
C***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION,       SLAP........496300
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........496400
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........496500
C           Seager, Mark K., (LLNL)                                      SLAP........496600
C             Lawrence Livermore National Laboratory                     SLAP........496700
C             PO BOX 808, L-60                                           SLAP........496800
C             Livermore, CA 94550 (510) 423-3141                         SLAP........496900
C             seager@llnl.gov                                            SLAP........497000
C***DESCRIPTION                                                          SLAP........497100
C                                                                        SLAP........497200
C *Usage:                                                                SLAP........497300
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM                          SLAP........497400
C     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)                     SLAP........497500
C     INTEGER NROW(N), NCOL(N)                                           SLAP........497600
C     DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)                    SLAP........497700
C                                                                        SLAP........497800
C     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,              SLAP........497900
C    $    DINV, NU, IU, JU, U, NROW, NCOL )                              SLAP........498000
C                                                                        SLAP........498100
C *Arguments:                                                            SLAP........498200
C N      :IN       Integer                                               SLAP........498300
C         Order of the Matrix.                                           SLAP........498400
C NELT   :IN       Integer.                                              SLAP........498500
C         Number of elements in arrays IA, JA, and A.                    SLAP........498600
C IA     :IN       Integer IA(NELT).                                     SLAP........498700
C JA     :IN       Integer JA(NELT).                                     SLAP........498800
C A      :IN       Double Precision A(NELT).                             SLAP........498900
C         These arrays should hold the matrix A in the SLAP Column       SLAP........499000
C         format.  See "Description", below.                             SLAP........499100
C ISYM   :IN       Integer.                                              SLAP........499200
C         Flag to indicate symmetric storage format.                     SLAP........499300
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........499400
C         If ISYM=1, the matrix is symmetric, and only the lower         SLAP........499500
C         triangle of the matrix is stored.                              SLAP........499600
C NL     :OUT      Integer.                                              SLAP........499700
C         Number of non-zeros in the L array.                            SLAP........499800
C IL     :OUT      Integer IL(NL).                                       SLAP........499900
C JL     :OUT      Integer JL(NL).                                       SLAP........500000
C L      :OUT      Double Precision L(NL).                               SLAP........500100
C         IL, JL, L  contain the unit lower triangular factor of  the    SLAP........500200
C         incomplete decomposition  of some  matrix stored  in   SLAP    SLAP........500300
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See    SLAP........500400
C         "DESCRIPTION", below for more details about the SLAP format.   SLAP........500500
C NU     :OUT      Integer.                                              SLAP........500600
C         Number of non-zeros in the U array.                            SLAP........500700
C IU     :OUT      Integer IU(NU).                                       SLAP........500800
C JU     :OUT      Integer JU(NU).                                       SLAP........500900
C U      :OUT      Double Precision     U(NU).                           SLAP........501000
C         IU, JU, U contain   the unit upper triangular factor of the    SLAP........501100
C         incomplete  decomposition    of some matrix  stored in SLAP    SLAP........501200
C         Column  format.   The Diagonal of ones   *IS*  stored.  See    SLAP........501300
C         "Description", below  for  more  details  about  the   SLAP    SLAP........501400
C         format.                                                        SLAP........501500
C NROW   :WORK     Integer NROW(N).                                      SLAP........501600
C         NROW(I) is the number of non-zero elements in the I-th row     SLAP........501700
C         of L.                                                          SLAP........501800
C NCOL   :WORK     Integer NCOL(N).                                      SLAP........501900
C         NCOL(I) is the number of non-zero elements in the I-th         SLAP........502000
C         column of U.                                                   SLAP........502100
C                                                                        SLAP........502200
C *Description                                                           SLAP........502300
C       IL, JL, L should contain the unit  lower triangular factor of    SLAP........502400
C       the incomplete decomposition of the A matrix  stored in SLAP     SLAP........502500
C       Row format.  IU, JU, U should contain  the unit upper factor     SLAP........502600
C       of the  incomplete decomposition of  the A matrix  stored in     SLAP........502700
C       SLAP Column format This ILU factorization can be computed by     SLAP........502800
C       the DSILUS routine. The diagonals (which are all one's) are      SLAP........502900
C       stored.                                                          SLAP........503000
C                                                                        SLAP........503100
C       =================== S L A P Column format ==================     SLAP........503200
C                                                                        SLAP........503300
C       This routine  requires that  the matrix A  be stored in  the     SLAP........503400
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........503500
C       counting down columns (except for  the diagonal entry, which     SLAP........503600
C       must appear first in each  "column")  and are stored  in the     SLAP........503700
C       double precision array A.   In other words,  for each column     SLAP........503800
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........503900
C       other non-zero  elements going down  the column (except  the     SLAP........504000
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........504100
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........504200
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........504300
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........504400
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........504500
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........504600
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........504700
C       number of columns in  the matrix and NELT  is the number  of     SLAP........504800
C       non-zeros in the matrix.                                         SLAP........504900
C                                                                        SLAP........505000
C       Here is an example of the  SLAP Column  storage format for a     SLAP........505100
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........505200
C       column):                                                         SLAP........505300
C                                                                        SLAP........505400
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........505500
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........505600
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........505700
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........505800
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........505900
C       | 0  0  0 44  0|                                                 SLAP........506000
C       |51  0 53  0 55|                                                 SLAP........506100
C                                                                        SLAP........506200
C       ==================== S L A P Row format ====================     SLAP........506300
C                                                                        SLAP........506400
C       This routine requires  that the matrix A  be  stored  in the     SLAP........506500
C       SLAP  Row format.   In this format  the non-zeros are stored     SLAP........506600
C       counting across  rows (except for the diagonal  entry, which     SLAP........506700
C       must  appear first  in each  "row")  and  are stored  in the     SLAP........506800
C       double precision  array A.  In other words, for each row  in     SLAP........506900
C       the matrix  put the diagonal  entry in A.   Then put in  the     SLAP........507000
C       other  non-zero elements  going across  the row  (except the     SLAP........507100
C       diagonal) in order.  The JA array holds the column index for     SLAP........507200
C       each non-zero.  The IA array holds the offsets  into the JA,     SLAP........507300
C       A  arrays  for  the   beginning  of  each  row.    That  is,     SLAP........507400
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-     SLAP........507500
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)     SLAP........507600
C       are  the last elements  of the  IROW-th row.   Note  that we     SLAP........507700
C       always have  IA(N+1) = NELT+1, where N is the number of rows     SLAP........507800
C       in the matrix  and  NELT is the  number of non-zeros  in the     SLAP........507900
C       matrix.                                                          SLAP........508000
C                                                                        SLAP........508100
C       Here is an example of the SLAP Row storage format for a  5x5     SLAP........508200
C       Matrix (in the A and JA arrays '|' denotes the end of a row):    SLAP........508300
C                                                                        SLAP........508400
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.   SLAP........508500
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........508600
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53   SLAP........508700
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........508800
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12                      SLAP........508900
C       | 0  0  0 44  0|                                                 SLAP........509000
C       |51  0 53  0 55|                                                 SLAP........509100
C                                                                        SLAP........509200
C***SEE ALSO  SILUR                                                      SLAP........509300
C***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations, SLAP........509400
C                  Johns Hopkins University Press, Baltimore, Maryland,  SLAP........509500
C                  1983.                                                 SLAP........509600
C***ROUTINES CALLED  (NONE)                                              SLAP........509700
C***REVISION HISTORY  (YYMMDD)                                           SLAP........509800
C   890404  DATE WRITTEN                                                 SLAP........509900
C   890404  Previous REVISION DATE                                       SLAP........510000
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........510100
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........510200
C           standard.  (FNF)                                             SLAP........510300
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........510400
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........510500
C   920511  Added complete declaration section.  (WRB)                   SLAP........510600
C   920929  Corrected format of reference.  (FNF)                        SLAP........510700
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........510800
C***END PROLOGUE  DSILUS                                                 SLAP........510900
C     .. Scalar Arguments ..                                             SLAP........511000
      INTEGER ISYM, N, NELT, NL, NU                                      SLAP........511100
C     .. Array Arguments ..                                              SLAP........511200
      DOUBLE PRECISION A(NELT), DINV(N), L(NL), U(NU)                    SLAP........511300
      INTEGER IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU),        SLAP........511400
     +        NCOL(N), NROW(N)                                           SLAP........511500
C     .. Local Scalars ..                                                SLAP........511600
      DOUBLE PRECISION TEMP                                              SLAP........511700
      INTEGER I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2,   SLAP........511800
     +        INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC,  SLAP........511900
     +        KR                                                         SLAP........512000
C***FIRST EXECUTABLE STATEMENT  DSILUS                                   SLAP........512100
C                                                                        SLAP........512200
C         Count number of elements in each row of the lower triangle.    SLAP........512300
C                                                                        SLAP........512400
      DO 10 I=1,N                                                        SLAP........512500
         NROW(I) = 0                                                     SLAP........512600
         NCOL(I) = 0                                                     SLAP........512700
 10   CONTINUE                                                           SLAP........512800
CVD$R NOCONCUR                                                           SLAP........512900
CVD$R NOVECTOR                                                           SLAP........513000
      DO 30 ICOL = 1, N                                                  SLAP........513100
         JBGN = JA(ICOL)+1                                               SLAP........513200
         JEND = JA(ICOL+1)-1                                             SLAP........513300
         IF( JBGN.LE.JEND ) THEN                                         SLAP........513400
            DO 20 J = JBGN, JEND                                         SLAP........513500
               IF( IA(J).LT.ICOL ) THEN                                  SLAP........513600
                  NCOL(ICOL) = NCOL(ICOL) + 1                            SLAP........513700
               ELSE                                                      SLAP........513800
                  NROW(IA(J)) = NROW(IA(J)) + 1                          SLAP........513900
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1          SLAP........514000
               ENDIF                                                     SLAP........514100
 20         CONTINUE                                                     SLAP........514200
         ENDIF                                                           SLAP........514300
 30   CONTINUE                                                           SLAP........514400
      JU(1) = 1                                                          SLAP........514500
      IL(1) = 1                                                          SLAP........514600
      DO 40 ICOL = 1, N                                                  SLAP........514700
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)                              SLAP........514800
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)                              SLAP........514900
         NROW(ICOL) = IL(ICOL)                                           SLAP........515000
         NCOL(ICOL) = JU(ICOL)                                           SLAP........515100
 40   CONTINUE                                                           SLAP........515200
C                                                                        SLAP........515300
C         Copy the matrix A into the L and U structures.                 SLAP........515400
      DO 60 ICOL = 1, N                                                  SLAP........515500
         DINV(ICOL) = A(JA(ICOL))                                        SLAP........515600
         JBGN = JA(ICOL)+1                                               SLAP........515700
         JEND = JA(ICOL+1)-1                                             SLAP........515800
         IF( JBGN.LE.JEND ) THEN                                         SLAP........515900
            DO 50 J = JBGN, JEND                                         SLAP........516000
               IROW = IA(J)                                              SLAP........516100
               IF( IROW.LT.ICOL ) THEN                                   SLAP........516200
C         Part of the upper triangle.                                    SLAP........516300
                  IU(NCOL(ICOL)) = IROW                                  SLAP........516400
                  U(NCOL(ICOL)) = A(J)                                   SLAP........516500
                  NCOL(ICOL) = NCOL(ICOL) + 1                            SLAP........516600
               ELSE                                                      SLAP........516700
C         Part of the lower triangle (stored by row).                    SLAP........516800
                  JL(NROW(IROW)) = ICOL                                  SLAP........516900
                  L(NROW(IROW)) = A(J)                                   SLAP........517000
                  NROW(IROW) = NROW(IROW) + 1                            SLAP........517100
                  IF( ISYM.NE.0 ) THEN                                   SLAP........517200
C         Symmetric...Copy lower triangle into upper triangle as well.   SLAP........517300
                     IU(NCOL(IROW)) = ICOL                               SLAP........517400
                     U(NCOL(IROW)) = A(J)                                SLAP........517500
                     NCOL(IROW) = NCOL(IROW) + 1                         SLAP........517600
                  ENDIF                                                  SLAP........517700
               ENDIF                                                     SLAP........517800
 50         CONTINUE                                                     SLAP........517900
         ENDIF                                                           SLAP........518000
 60   CONTINUE                                                           SLAP........518100
C                                                                        SLAP........518200
C         Sort the rows of L and the columns of U.                       SLAP........518300
      DO 110 K = 2, N                                                    SLAP........518400
         JBGN = JU(K)                                                    SLAP........518500
         JEND = JU(K+1)-1                                                SLAP........518600
         IF( JBGN.LT.JEND ) THEN                                         SLAP........518700
            DO 80 J = JBGN, JEND-1                                       SLAP........518800
               DO 70 I = J+1, JEND                                       SLAP........518900
                  IF( IU(J).GT.IU(I) ) THEN                              SLAP........519000
                     ITEMP = IU(J)                                       SLAP........519100
                     IU(J) = IU(I)                                       SLAP........519200
                     IU(I) = ITEMP                                       SLAP........519300
                     TEMP = U(J)                                         SLAP........519400
                     U(J) = U(I)                                         SLAP........519500
                     U(I) = TEMP                                         SLAP........519600
                  ENDIF                                                  SLAP........519700
 70            CONTINUE                                                  SLAP........519800
 80         CONTINUE                                                     SLAP........519900
         ENDIF                                                           SLAP........520000
         IBGN = IL(K)                                                    SLAP........520100
         IEND = IL(K+1)-1                                                SLAP........520200
         IF( IBGN.LT.IEND ) THEN                                         SLAP........520300
            DO 100 I = IBGN, IEND-1                                      SLAP........520400
               DO 90 J = I+1, IEND                                       SLAP........520500
                  IF( JL(I).GT.JL(J) ) THEN                              SLAP........520600
                     JTEMP = JU(I)                                       SLAP........520700
                     JU(I) = JU(J)                                       SLAP........520800
                     JU(J) = JTEMP                                       SLAP........520900
                     TEMP = L(I)                                         SLAP........521000
                     L(I) = L(J)                                         SLAP........521100
                     L(J) = TEMP                                         SLAP........521200
                  ENDIF                                                  SLAP........521300
 90            CONTINUE                                                  SLAP........521400
 100        CONTINUE                                                     SLAP........521500
         ENDIF                                                           SLAP........521600
 110  CONTINUE                                                           SLAP........521700
C                                                                        SLAP........521800
C         Perform the incomplete LDU decomposition.                      SLAP........521900
      DO 300 I=2,N                                                       SLAP........522000
C                                                                        SLAP........522100
C           I-th row of L                                                SLAP........522200
         INDX1 = IL(I)                                                   SLAP........522300
         INDX2 = IL(I+1) - 1                                             SLAP........522400
         IF(INDX1 .GT. INDX2) GO TO 200                                  SLAP........522500
         DO 190 INDX=INDX1,INDX2                                         SLAP........522600
            IF(INDX .EQ. INDX1) GO TO 180                                SLAP........522700
            INDXR1 = INDX1                                               SLAP........522800
            INDXR2 = INDX - 1                                            SLAP........522900
            INDXC1 = JU(JL(INDX))                                        SLAP........523000
            INDXC2 = JU(JL(INDX)+1) - 1                                  SLAP........523100
            IF(INDXC1 .GT. INDXC2) GO TO 180                             SLAP........523200
 160        KR = JL(INDXR1)                                              SLAP........523300
 170        KC = IU(INDXC1)                                              SLAP........523400
            IF(KR .GT. KC) THEN                                          SLAP........523500
               INDXC1 = INDXC1 + 1                                       SLAP........523600
               IF(INDXC1 .LE. INDXC2) GO TO 170                          SLAP........523700
            ELSEIF(KR .LT. KC) THEN                                      SLAP........523800
               INDXR1 = INDXR1 + 1                                       SLAP........523900
               IF(INDXR1 .LE. INDXR2) GO TO 160                          SLAP........524000
            ELSEIF(KR .EQ. KC) THEN                                      SLAP........524100
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)          SLAP........524200
               INDXR1 = INDXR1 + 1                                       SLAP........524300
               INDXC1 = INDXC1 + 1                                       SLAP........524400
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160 SLAP........524500
            ENDIF                                                        SLAP........524600
 180        L(INDX) = L(INDX)/DINV(JL(INDX))                             SLAP........524700
 190     CONTINUE                                                        SLAP........524800
C                                                                        SLAP........524900
C         I-th column of U                                               SLAP........525000
 200     INDX1 = JU(I)                                                   SLAP........525100
         INDX2 = JU(I+1) - 1                                             SLAP........525200
         IF(INDX1 .GT. INDX2) GO TO 260                                  SLAP........525300
         DO 250 INDX=INDX1,INDX2                                         SLAP........525400
            IF(INDX .EQ. INDX1) GO TO 240                                SLAP........525500
            INDXC1 = INDX1                                               SLAP........525600
            INDXC2 = INDX - 1                                            SLAP........525700
            INDXR1 = IL(IU(INDX))                                        SLAP........525800
            INDXR2 = IL(IU(INDX)+1) - 1                                  SLAP........525900
            IF(INDXR1 .GT. INDXR2) GO TO 240                             SLAP........526000
 210        KR = JL(INDXR1)                                              SLAP........526100
 220        KC = IU(INDXC1)                                              SLAP........526200
            IF(KR .GT. KC) THEN                                          SLAP........526300
               INDXC1 = INDXC1 + 1                                       SLAP........526400
               IF(INDXC1 .LE. INDXC2) GO TO 220                          SLAP........526500
            ELSEIF(KR .LT. KC) THEN                                      SLAP........526600
               INDXR1 = INDXR1 + 1                                       SLAP........526700
               IF(INDXR1 .LE. INDXR2) GO TO 210                          SLAP........526800
            ELSEIF(KR .EQ. KC) THEN                                      SLAP........526900
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)          SLAP........527000
               INDXR1 = INDXR1 + 1                                       SLAP........527100
               INDXC1 = INDXC1 + 1                                       SLAP........527200
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210 SLAP........527300
            ENDIF                                                        SLAP........527400
 240        U(INDX) = U(INDX)/DINV(IU(INDX))                             SLAP........527500
 250     CONTINUE                                                        SLAP........527600
C                                                                        SLAP........527700
C         I-th diagonal element                                          SLAP........527800
 260     INDXR1 = IL(I)                                                  SLAP........527900
         INDXR2 = IL(I+1) - 1                                            SLAP........528000
         IF(INDXR1 .GT. INDXR2) GO TO 300                                SLAP........528100
         INDXC1 = JU(I)                                                  SLAP........528200
         INDXC2 = JU(I+1) - 1                                            SLAP........528300
         IF(INDXC1 .GT. INDXC2) GO TO 300                                SLAP........528400
 270     KR = JL(INDXR1)                                                 SLAP........528500
 280     KC = IU(INDXC1)                                                 SLAP........528600
         IF(KR .GT. KC) THEN                                             SLAP........528700
            INDXC1 = INDXC1 + 1                                          SLAP........528800
            IF(INDXC1 .LE. INDXC2) GO TO 280                             SLAP........528900
         ELSEIF(KR .LT. KC) THEN                                         SLAP........529000
            INDXR1 = INDXR1 + 1                                          SLAP........529100
            IF(INDXR1 .LE. INDXR2) GO TO 270                             SLAP........529200
         ELSEIF(KR .EQ. KC) THEN                                         SLAP........529300
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)             SLAP........529400
            INDXR1 = INDXR1 + 1                                          SLAP........529500
            INDXC1 = INDXC1 + 1                                          SLAP........529600
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270    SLAP........529700
         ENDIF                                                           SLAP........529800
C                                                                        SLAP........529900
 300  CONTINUE                                                           SLAP........530000
C                                                                        SLAP........530100
C         Replace diagonal elements by their inverses.                   SLAP........530200
CVD$ VECTOR                                                              SLAP........530300
      DO 430 I=1,N                                                       SLAP........530400
         DINV(I) = 1.0D0/DINV(I)                                         SLAP........530500
 430  CONTINUE                                                           SLAP........530600
C                                                                        SLAP........530700
      RETURN                                                             SLAP........530800
C------------- LAST LINE OF DSILUS FOLLOWS ----------------------------  SLAP........530900
      END                                                                SLAP........531000
*DECK DSLLTI                                                             SLAP........531100
      SUBROUTINE DSLLTI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)   SLAP........531200
C***BEGIN PROLOGUE  DSLLTI                                               SLAP........531300
C***PURPOSE  SLAP MSOLVE for LDL' (IC) Factorization.                    SLAP........531400
C            This routine acts as an interface between the SLAP generic  SLAP........531500
C            MSOLVE calling convention and the routine that actually     SLAP........531600
C                           -1                                           SLAP........531700
C            computes (LDL')  B = X.                                     SLAP........531800
C***LIBRARY   SLATEC (SLAP)                                              SLAP........531900
C***CATEGORY  D2E                                                        SLAP........532000
C***TYPE      DOUBLE PRECISION (SSLLTI-S, DSLLTI-D)                      SLAP........532100
C***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE  SLAP........532200
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........532300
C           Seager, Mark K., (LLNL)                                      SLAP........532400
C             Lawrence Livermore National Laboratory                     SLAP........532500
C             PO BOX 808, L-60                                           SLAP........532600
C             Livermore, CA 94550 (510) 423-3141                         SLAP........532700
C             seager@llnl.gov                                            SLAP........532800
C***DESCRIPTION                                                          SLAP........532900
C       It is assumed that RWORK and IWORK have initialized with         SLAP........533000
C       the information required for DLLTI2:                             SLAP........533100
C          IWORK(1) = NEL                                                SLAP........533200
C          IWORK(2) = Starting location of IEL in IWORK.                 SLAP........533300
C          IWORK(3) = Starting location of JEL in IWORK.                 SLAP........533400
C          IWORK(4) = Starting location of EL in RWORK.                  SLAP........533500
C          IWORK(5) = Starting location of DINV in RWORK.                SLAP........533600
C       See the DESCRIPTION of DLLTI2 for details.                       SLAP........533700
C***REFERENCES  (NONE)                                                   SLAP........533800
C***ROUTINES CALLED  DLLTI2                                              SLAP........533900
C***REVISION HISTORY  (YYMMDD)                                           SLAP........534000
C   871119  DATE WRITTEN                                                 SLAP........534100
C   881213  Previous REVISION DATE                                       SLAP........534200
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........534300
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........534400
C           standard.  (FNF)                                             SLAP........534500
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........534600
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........534700
C   910502  Corrected conversion error.  (FNF)                           SLAP........534800
C   920511  Added complete declaration section.  (WRB)                   SLAP........534900
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........535000
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........535100
C***END PROLOGUE  DSLLTI                                                 SLAP........535200
C     .. Scalar Arguments ..                                             SLAP........535300
      INTEGER ISYM, N, NELT                                              SLAP........535400
C     .. Array Arguments ..                                              SLAP........535500
      DOUBLE PRECISION A(NELT), B(*), RWORK(*), X(*)                     SLAP........535600
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........535700
C     .. Local Scalars ..                                                SLAP........535800
      INTEGER LOCDIN, LOCEL, LOCIEL, LOCJEL, NEL                         SLAP........535900
C     .. External Subroutines ..                                         SLAP........536000
      EXTERNAL DLLTI2                                                    SLAP........536100
C***FIRST EXECUTABLE STATEMENT  DSLLTI                                   SLAP........536200
      NEL = IWORK(1)                                                     SLAP........536300
      LOCIEL = IWORK(3)                                                  SLAP........536400
      LOCJEL = IWORK(2)                                                  SLAP........536500
      LOCEL  = IWORK(4)                                                  SLAP........536600
      LOCDIN = IWORK(5)                                                  SLAP........536700
      CALL DLLTI2(N, B, X, NEL, IWORK(LOCIEL), IWORK(LOCJEL),            SLAP........536800
     $     RWORK(LOCEL), RWORK(LOCDIN))                                  SLAP........536900
C                                                                        SLAP........537000
      RETURN                                                             SLAP........537100
C------------- LAST LINE OF DSLLTI FOLLOWS ----------------------------  SLAP........537200
      END                                                                SLAP........537300
*DECK DSLUGM                                                             SLAP........537400
      SUBROUTINE DSLUGM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,    SLAP........537500
     +   TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)  SLAP........537600
C***BEGIN PROLOGUE  DSLUGM                                               SLAP........537700
C***PURPOSE  Incomplete LU GMRES iterative sparse Ax=b solver.           SLAP........537800
C            This routine uses the generalized minimum residual          SLAP........537900
C            (GMRES) method with incomplete LU factorization for         SLAP........538000
C            preconditioning to solve possibly non-symmetric linear      SLAP........538100
C            systems of the form: Ax = b.                                SLAP........538200
C***LIBRARY   SLATEC (SLAP)                                              SLAP........538300
C***CATEGORY  D2A4, D2B4                                                 SLAP........538400
C***TYPE      DOUBLE PRECISION (SSLUGM-S, DSLUGM-D)                      SLAP........538500
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........538600
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........538700
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........538800
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........538900
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........539000
C             Lawrence Livermore National Laboratory                     SLAP........539100
C             PO Box 808, L-60                                           SLAP........539200
C             Livermore, CA 94550 (510) 423-3141                         SLAP........539300
C***DESCRIPTION                                                          SLAP........539400
C                                                                        SLAP........539500
C *Usage:                                                                SLAP........539600
C      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL          SLAP........539700
C      INTEGER   ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW     SLAP........539800
C      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)       SLAP........539900
C                                                                        SLAP........540000
C      CALL DSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE,                SLAP........540100
C     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,                    SLAP........540200
C     $     RWORK, LENW, IWORK, LENIW)                                   SLAP........540300
C                                                                        SLAP........540400
C *Arguments:                                                            SLAP........540500
C N      :IN       Integer.                                              SLAP........540600
C         Order of the Matrix.                                           SLAP........540700
C B      :IN       Double Precision B(N).                                SLAP........540800
C         Right-hand side vector.                                        SLAP........540900
C X      :INOUT    Double Precision X(N).                                SLAP........541000
C         On input X is your initial guess for solution vector.          SLAP........541100
C         On output X is the final approximate solution.                 SLAP........541200
C NELT   :IN       Integer.                                              SLAP........541300
C         Number of Non-Zeros stored in A.                               SLAP........541400
C IA     :IN       Integer IA(NELT).                                     SLAP........541500
C JA     :IN       Integer JA(NELT).                                     SLAP........541600
C A      :IN       Double Precision A(NELT).                             SLAP........541700
C         These arrays should hold the matrix A in either the SLAP       SLAP........541800
C         Triad format or the SLAP Column format.  See "Description",    SLAP........541900
C         below.  If the SLAP Triad format is chosen it is changed       SLAP........542000
C         internally to the SLAP Column format.                          SLAP........542100
C ISYM   :IN       Integer.                                              SLAP........542200
C         Flag to indicate symmetric storage format.                     SLAP........542300
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........542400
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........542500
C         or lower triangle of the matrix is stored.                     SLAP........542600
C NSAVE  :IN       Integer.                                              SLAP........542700
C         Number of direction vectors to save and orthogonalize against. SLAP........542800
C         Must be greater than 1.                                        SLAP........542900
C ITOL   :IN       Integer.                                              SLAP........543000
C         Flag to indicate the type of convergence criterion used.       SLAP........543100
C         ITOL=0  Means the  iteration stops when the test described     SLAP........543200
C                 below on  the  residual RL  is satisfied.  This is     SLAP........543300
C                 the  "Natural Stopping Criteria" for this routine.     SLAP........543400
C                 Other values  of   ITOL  cause  extra,   otherwise     SLAP........543500
C                 unnecessary, computation per iteration and     are     SLAP........543600
C                 therefore  much less  efficient.  See  ISDGMR (the     SLAP........543700
C                 stop test routine) for more information.               SLAP........543800
C         ITOL=1  Means   the  iteration stops   when the first test     SLAP........543900
C                 described below on  the residual RL  is satisfied,     SLAP........544000
C                 and there  is either right  or  no preconditioning     SLAP........544100
C                 being used.                                            SLAP........544200
C         ITOL=2  Implies     that   the  user    is   using    left     SLAP........544300
C                 preconditioning, and the second stopping criterion     SLAP........544400
C                 below is used.                                         SLAP........544500
C         ITOL=3  Means the  iteration stops   when  the  third test     SLAP........544600
C                 described below on Minv*Residual is satisfied, and     SLAP........544700
C                 there is either left  or no  preconditioning begin     SLAP........544800
C                 used.                                                  SLAP........544900
C         ITOL=11 is    often  useful  for   checking  and comparing     SLAP........545000
C                 different routines.  For this case, the  user must     SLAP........545100
C                 supply  the  "exact" solution or  a  very accurate     SLAP........545200
C                 approximation (one with  an  error much less  than     SLAP........545300
C                 TOL) through a common block,                           SLAP........545400
C                     COMMON /DSLBLK/ SOLN( )                            SLAP........545500
C                 If ITOL=11, iteration stops when the 2-norm of the     SLAP........545600
C                 difference between the iterative approximation and     SLAP........545700
C                 the user-supplied solution  divided by the  2-norm     SLAP........545800
C                 of the  user-supplied solution  is  less than TOL.     SLAP........545900
C                 Note that this requires  the  user to  set up  the     SLAP........546000
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling     SLAP........546100
C                 routine.  The routine with this declaration should     SLAP........546200
C                 be loaded before the stop test so that the correct     SLAP........546300
C                 length is used by  the loader.  This procedure  is     SLAP........546400
C                 not standard Fortran and may not work correctly on     SLAP........546500
C                 your   system (although  it  has  worked  on every     SLAP........546600
C                 system the authors have tried).  If ITOL is not 11     SLAP........546700
C                 then this common block is indeed standard Fortran.     SLAP........546800
C TOL    :INOUT    Double Precision.                                     SLAP........546900
C         Convergence criterion, as described below.  If TOL is set      SLAP........547000
C         to zero on input, then a default value of 500*(the smallest    SLAP........547100
C         positive magnitude, machine epsilon) is used.                  SLAP........547200
C ITMAX  :IN       Integer.                                              SLAP........547300
C         Maximum number of iterations.  This routine uses the default   SLAP........547400
C         of NRMAX = ITMAX/NSAVE to determine the when each restart      SLAP........547500
C         should occur.  See the description of NRMAX and MAXL in        SLAP........547600
C         DGMRES for a full and frightfully interesting discussion of    SLAP........547700
C         this topic.                                                    SLAP........547800
C ITER   :OUT      Integer.                                              SLAP........547900
C         Number of iterations required to reach convergence, or         SLAP........548000
C         ITMAX+1 if convergence criterion could not be achieved in      SLAP........548100
C         ITMAX iterations.                                              SLAP........548200
C ERR    :OUT      Double Precision.                                     SLAP........548300
C         Error estimate of error in final approximate solution, as      SLAP........548400
C         defined by ITOL.  Letting norm() denote the Euclidean          SLAP........548500
C         norm, ERR is defined as follows...                             SLAP........548600
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........548700
C                               for right or no preconditioning, and     SLAP........548800
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........548900
C                                norm(SB*(M-inverse)*B),                 SLAP........549000
C                               for left preconditioning.                SLAP........549100
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........549200
C                               since right or no preconditioning        SLAP........549300
C                               being used.                              SLAP........549400
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........549500
C                                norm(SB*(M-inverse)*B),                 SLAP........549600
C                               since left preconditioning is being      SLAP........549700
C                               used.                                    SLAP........549800
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|        SLAP........549900
C                               i=1,n                                    SLAP........550000
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).     SLAP........550100
C IERR   :OUT      Integer.                                              SLAP........550200
C         Return error flag.                                             SLAP........550300
C               IERR = 0 => All went well.                               SLAP........550400
C               IERR = 1 => Insufficient storage allocated for           SLAP........550500
C                           RGWK or IGWK.                                SLAP........550600
C               IERR = 2 => Routine DPIGMR failed to reduce the norm     SLAP........550700
C                           of the current residual on its last call,    SLAP........550800
C                           and so the iteration has stalled.  In        SLAP........550900
C                           this case, X equals the last computed        SLAP........551000
C                           approximation.  The user must either         SLAP........551100
C                           increase MAXL, or choose a different         SLAP........551200
C                           initial guess.                               SLAP........551300
C               IERR =-1 => Insufficient length for RGWK array.          SLAP........551400
C                           IGWK(6) contains the required minimum        SLAP........551500
C                           length of the RGWK array.                    SLAP........551600
C               IERR =-2 => Inconsistent ITOL and JPRE values.           SLAP........551700
C         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the        SLAP........551800
C         left-hand-side of the relevant stopping test defined           SLAP........551900
C         below associated with the residual for the current             SLAP........552000
C         approximation X(L).                                            SLAP........552100
C IUNIT  :IN       Integer.                                              SLAP........552200
C         Unit number on which to write the error at each iteration,     SLAP........552300
C         if this is desired for monitoring convergence.  If unit        SLAP........552400
C         number is 0, no writing will occur.                            SLAP........552500
C RWORK  :WORK    Double Precision RWORK(LENW).                          SLAP........552600
C         Double Precision array of size LENW.                           SLAP........552700
C LENW   :IN       Integer.                                              SLAP........552800
C         Length of the double precision workspace, RWORK.               SLAP........552900
C         LENW >= 1 + N*(NSAVE+7) +  NSAVE*(NSAVE+3)+NL+NU.              SLAP........553000
C         Here NL is the number of non-zeros in the lower triangle of    SLAP........553100
C         the matrix (including the diagonal) and NU is the number of    SLAP........553200
C         non-zeros in the upper triangle of the matrix (including the   SLAP........553300
C         diagonal).                                                     SLAP........553400
C         For the recommended values,  RWORK  has size at least          SLAP........553500
C         131 + 17*N + NL + NU.                                          SLAP........553600
C IWORK  :INOUT    Integer IWORK(LENIW).                                 SLAP........553700
C         Used to hold pointers into the RWORK array.                    SLAP........553800
C         Upon return the following locations of IWORK hold information  SLAP........553900
C         which may be of use to the user:                               SLAP........554000
C         IWORK(9)  Amount of Integer workspace actually used.           SLAP........554100
C         IWORK(10) Amount of Double Precision workspace actually used.  SLAP........554200
C LENIW  :IN       Integer.                                              SLAP........554300
C         Length of the integer workspace, IWORK.                        SLAP........554400
C         LENIW >= NL+NU+4*N+32.                                         SLAP........554500
C                                                                        SLAP........554600
C *Description:                                                          SLAP........554700
C       DSLUGM solves a linear system A*X = B rewritten in the form:     SLAP........554800
C                                                                        SLAP........554900
C        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B,                SLAP........555000
C                                                                        SLAP........555100
C       with right preconditioning, or                                   SLAP........555200
C                                                                        SLAP........555300
C        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B,      SLAP........555400
C                                                                        SLAP........555500
C       with left preconditioning, where A is an n-by-n double precision SLAP........555600
C       matrix, X and B are N-vectors, SB and SX are  diagonal scaling   SLAP........555700
C       matrices, and M is the Incomplete LU factorization of A.  It     SLAP........555800
C       uses preconditioned  Krylov subpace   methods  based on  the     SLAP........555900
C       generalized minimum residual  method (GMRES).   This routine     SLAP........556000
C       is a  driver  routine  which  assumes a SLAP   matrix   data     SLAP........556100
C       structure   and  sets  up  the  necessary  information to do     SLAP........556200
C       diagonal  preconditioning  and calls the main GMRES  routine     SLAP........556300
C       DGMRES for the   solution   of the linear   system.   DGMRES     SLAP........556400
C       optionally   performs  either  the full    orthogonalization     SLAP........556500
C       version of the  GMRES algorithm or  an incomplete variant of     SLAP........556600
C       it.  Both versions use restarting of the linear iteration by     SLAP........556700
C       default, although the user can disable this feature.             SLAP........556800
C                                                                        SLAP........556900
C       The GMRES  algorithm generates a sequence  of approximations     SLAP........557000
C       X(L) to the  true solution of the above  linear system.  The     SLAP........557100
C       convergence criteria for stopping the  iteration is based on     SLAP........557200
C       the size  of the  scaled norm of  the residual  R(L)  =  B -     SLAP........557300
C       A*X(L).  The actual stopping test is either:                     SLAP........557400
C                                                                        SLAP........557500
C               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B),                 SLAP........557600
C                                                                        SLAP........557700
C       for right preconditioning, or                                    SLAP........557800
C                                                                        SLAP........557900
C               norm(SB*(M-inverse)*(B-A*X(L))) .le.                     SLAP........558000
C                       TOL*norm(SB*(M-inverse)*B),                      SLAP........558100
C                                                                        SLAP........558200
C       for left preconditioning, where norm() denotes the Euclidean     SLAP........558300
C       norm, and TOL is  a positive scalar less  than one  input by     SLAP........558400
C       the user.  If TOL equals zero  when DSLUGM is called, then a     SLAP........558500
C       default  value  of 500*(the   smallest  positive  magnitude,     SLAP........558600
C       machine epsilon) is used.  If the  scaling arrays SB  and SX     SLAP........558700
C       are used, then  ideally they  should be chosen  so  that the     SLAP........558800
C       vectors SX*X(or SX*M*X) and  SB*B have all their  components     SLAP........558900
C       approximately equal  to  one in  magnitude.  If one wants to     SLAP........559000
C       use the same scaling in X  and B, then  SB and SX can be the     SLAP........559100
C       same array in the calling program.                               SLAP........559200
C                                                                        SLAP........559300
C       The following is a list of the other routines and their          SLAP........559400
C       functions used by GMRES:                                         SLAP........559500
C       DGMRES  Contains the matrix structure independent driver         SLAP........559600
C               routine for GMRES.                                       SLAP........559700
C       DPIGMR  Contains the main iteration loop for GMRES.              SLAP........559800
C       DORTH   Orthogonalizes a new vector against older basis vectors. SLAP........559900
C       DHEQR   Computes a QR decomposition of a Hessenberg matrix.      SLAP........560000
C       DHELS   Solves a Hessenberg least-squares system, using QR       SLAP........560100
C               factors.                                                 SLAP........560200
C       RLCALC  Computes the scaled residual RL.                         SLAP........560300
C       XLCALC  Computes the solution XL.                                SLAP........560400
C       ISDGMR  User-replaceable stopping routine.                       SLAP........560500
C                                                                        SLAP........560600
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix     SLAP........560700
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP     SLAP........560800
C       Column format.  The user can hand this routine either of the     SLAP........560900
C       of these data structures and SLAP  will figure out  which on     SLAP........561000
C       is being used and act accordingly.                               SLAP........561100
C                                                                        SLAP........561200
C       =================== S L A P Triad format ===================     SLAP........561300
C       This routine requires that the  matrix A be   stored in  the     SLAP........561400
C       SLAP  Triad format.  In  this format only the non-zeros  are     SLAP........561500
C       stored.  They may appear in  *ANY* order.  The user supplies     SLAP........561600
C       three arrays of  length NELT, where  NELT is  the number  of     SLAP........561700
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For     SLAP........561800
C       each non-zero the user puts the row and column index of that     SLAP........561900
C       matrix element  in the IA and  JA arrays.  The  value of the     SLAP........562000
C       non-zero   matrix  element is  placed  in  the corresponding     SLAP........562100
C       location of the A array.   This is  an  extremely  easy data     SLAP........562200
C       structure to generate.  On  the  other hand it   is  not too     SLAP........562300
C       efficient on vector computers for  the iterative solution of     SLAP........562400
C       linear systems.  Hence,   SLAP changes   this  input    data     SLAP........562500
C       structure to the SLAP Column format  for  the iteration (but     SLAP........562600
C       does not change it back).                                        SLAP........562700
C                                                                        SLAP........562800
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........562900
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........563000
C                                                                        SLAP........563100
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........563200
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........563300
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........563400
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........563500
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........563600
C       | 0  0  0 44  0|                                                 SLAP........563700
C       |51  0 53  0 55|                                                 SLAP........563800
C                                                                        SLAP........563900
C       =================== S L A P Column format ==================     SLAP........564000
C                                                                        SLAP........564100
C       This routine  requires that  the matrix A  be stored in  the     SLAP........564200
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........564300
C       counting down columns (except for  the diagonal entry, which     SLAP........564400
C       must appear first in each  "column")  and are stored  in the     SLAP........564500
C       double precision array A.   In other words,  for each column     SLAP........564600
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........564700
C       other non-zero  elements going down  the column (except  the     SLAP........564800
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........564900
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........565000
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........565100
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........565200
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........565300
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........565400
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........565500
C       number of columns in  the matrix and NELT  is the number  of     SLAP........565600
C       non-zeros in the matrix.                                         SLAP........565700
C                                                                        SLAP........565800
C       Here is an example of the  SLAP Column  storage format for a     SLAP........565900
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........566000
C       column):                                                         SLAP........566100
C                                                                        SLAP........566200
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........566300
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........566400
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........566500
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........566600
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........566700
C       | 0  0  0 44  0|                                                 SLAP........566800
C       |51  0 53  0 55|                                                 SLAP........566900
C                                                                        SLAP........567000
C *Side Effects:                                                         SLAP........567100
C       The SLAP Triad format (IA, JA, A) is modified internally to be   SLAP........567200
C       the SLAP Column format.  See above.                              SLAP........567300
C                                                                        SLAP........567400
C *Cautions:                                                             SLAP........567500
C     This routine will attempt to write to the Fortran logical output   SLAP........567600
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........567700
C     this logical unit is attached to a file or terminal before calling SLAP........567800
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........567900
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........568000
C                                                                        SLAP........568100
C***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, Reduced Storage   SLAP........568200
C                  Matrix Methods in Stiff ODE Systems, Lawrence Liver-  SLAP........568300
C                  more National Laboratory Report UCRL-95088, Rev. 1,   SLAP........568400
C                  Livermore, California, June 1987.                     SLAP........568500
C***ROUTINES CALLED  DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV            SLAP........568600
C***REVISION HISTORY  (YYMMDD)                                           SLAP........568700
C   890404  DATE WRITTEN                                                 SLAP........568800
C   890404  Previous REVISION DATE                                       SLAP........568900
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........569000
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........569100
C           standard.  (FNF)                                             SLAP........569200
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........569300
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........569400
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........569500
C   920511  Added complete declaration section.  (WRB)                   SLAP........569600
C   920929  Corrected format of references.  (FNF)                       SLAP........569700
C   921019  Corrected NEL to NL.  (FNF)                                  SLAP........569800
C***END PROLOGUE  DSLUGM                                                 SLAP........569900
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........570000
CLLL. OPTIMIZE                                                           SLAP........570100
C     .. Parameters ..                                                   SLAP........570200
      INTEGER LOCRB, LOCIB                                               SLAP........570300
      PARAMETER (LOCRB=1, LOCIB=11)                                      SLAP........570400
C     .. Scalar Arguments ..                                             SLAP........570500
      DOUBLE PRECISION ERR, TOL                                          SLAP........570600
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,      SLAP........570700
     +        NELT, NSAVE                                                SLAP........570800
C     .. Array Arguments ..                                              SLAP........570900
      DOUBLE PRECISION A(NELT), B(N), RWORK(LENW), X(N)                  SLAP........571000
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)                           SLAP........571100
C     .. Local Scalars ..                                                SLAP........571200
      INTEGER ICOL, J, JBGN, JEND, LOCDIN, LOCIGW, LOCIL, LOCIU, LOCIW,  SLAP........571300
     +        LOCJL, LOCJU, LOCL, LOCNC, LOCNR, LOCRGW, LOCU, LOCW,      SLAP........571400
     +        MYITOL, NL, NU                                             SLAP........571500
C     .. External Subroutines ..                                         SLAP........571600
      EXTERNAL DCHKW, DGMRES, DS2Y, DSILUS, DSLUI, DSMV                  SLAP........571700
C***FIRST EXECUTABLE STATEMENT  DSLUGM                                   SLAP........571800
C                                                                        SLAP........571900
      IERR = 0                                                           SLAP........572000
      ERR  = 0                                                           SLAP........572100
      IF( NSAVE.LE.1 ) THEN                                              SLAP........572200
         IERR = 3                                                        SLAP........572300
         RETURN                                                          SLAP........572400
      ENDIF                                                              SLAP........572500
C                                                                        SLAP........572600
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.  SLAP........572700
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )                              SLAP........572800
C                                                                        SLAP........572900
C         Count number of Non-Zero elements preconditioner ILU matrix.   SLAP........573000
C         Then set up the work arrays.  We assume MAXL=KMP=NSAVE.        SLAP........573100
      NL = 0                                                             SLAP........573200
      NU = 0                                                             SLAP........573300
      DO 20 ICOL = 1, N                                                  SLAP........573400
C         Don't count diagonal.                                          SLAP........573500
         JBGN = JA(ICOL)+1                                               SLAP........573600
         JEND = JA(ICOL+1)-1                                             SLAP........573700
         IF( JBGN.LE.JEND ) THEN                                         SLAP........573800
CVD$ NOVECTOR                                                            SLAP........573900
            DO 10 J = JBGN, JEND                                         SLAP........574000
               IF( IA(J).GT.ICOL ) THEN                                  SLAP........574100
                  NL = NL + 1                                            SLAP........574200
                  IF( ISYM.NE.0 ) NU = NU + 1                            SLAP........574300
               ELSE                                                      SLAP........574400
                  NU = NU + 1                                            SLAP........574500
               ENDIF                                                     SLAP........574600
 10         CONTINUE                                                     SLAP........574700
         ENDIF                                                           SLAP........574800
 20   CONTINUE                                                           SLAP........574900
C                                                                        SLAP........575000
      LOCIGW = LOCIB                                                     SLAP........575100
      LOCIL = LOCIGW + 20                                                SLAP........575200
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........575300
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........575400
C        LOCJL = LOCIL + N+1                                             SLAP........575500
      LOCJL = LOCIL + NL                                                 SLAP........575600
      LOCIU = LOCJL + NL                                                 SLAP........575700
      LOCJU = LOCIU + NU                                                 SLAP........575800
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........575900
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........576000
C        LOCNR = LOCJU + N+1                                             SLAP........576100
      LOCNR = LOCJU + NU                                                 SLAP........576200
      LOCNC = LOCNR + N                                                  SLAP........576300
      LOCIW = LOCNC + N                                                  SLAP........576400
C                                                                        SLAP........576500
      LOCL = LOCRB                                                       SLAP........576600
      LOCDIN = LOCL + NL                                                 SLAP........576700
      LOCU = LOCDIN + N                                                  SLAP........576800
      LOCRGW = LOCU + NU                                                 SLAP........576900
      LOCW = LOCRGW + 1+N*(NSAVE+6)+NSAVE*(NSAVE+3)                      SLAP........577000
C                                                                        SLAP........577100
C         Check the workspace allocations.                               SLAP........577200
      CALL DCHKW( 'DSLUGM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )  SLAP........577300
      IF( IERR.NE.0 ) RETURN                                             SLAP........577400
C                                                                        SLAP........577500
      IWORK(1) = LOCIL                                                   SLAP........577600
      IWORK(2) = LOCJL                                                   SLAP........577700
      IWORK(3) = LOCIU                                                   SLAP........577800
      IWORK(4) = LOCJU                                                   SLAP........577900
      IWORK(5) = LOCL                                                    SLAP........578000
      IWORK(6) = LOCDIN                                                  SLAP........578100
      IWORK(7) = LOCU                                                    SLAP........578200
      IWORK(9) = LOCIW                                                   SLAP........578300
      IWORK(10) = LOCW                                                   SLAP........578400
C                                                                        SLAP........578500
C         Compute the Incomplete LU decomposition.                       SLAP........578600
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),           SLAP........578700
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),   SLAP........578800
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )       SLAP........578900
C                                                                        SLAP........579000
C         Perform the Incomplete LU Preconditioned Generalized Minimum   SLAP........579100
C         Residual iteration algorithm.  The following DGMRES            SLAP........579200
C         defaults are used MAXL = KMP = NSAVE, JSCAL = 0,               SLAP........579300
C         JPRE = -1, NRMAX = ITMAX/NSAVE                                 SLAP........579400
      IWORK(LOCIGW  ) = NSAVE                                            SLAP........579500
      IWORK(LOCIGW+1) = NSAVE                                            SLAP........579600
      IWORK(LOCIGW+2) = 0                                                SLAP........579700
      IWORK(LOCIGW+3) = -1                                               SLAP........579800
      IWORK(LOCIGW+4) = ITMAX/NSAVE                                      SLAP........579900
      MYITOL = 0                                                         SLAP........580000
C                                                                        SLAP........580100
      CALL DGMRES( N, B, X, NELT, IA, JA, A, ISYM, DSMV, DSLUI,          SLAP........580200
     $     MYITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, RWORK,     SLAP........580300
     $     RWORK(LOCRGW), LENW-LOCRGW, IWORK(LOCIGW), 20,                SLAP........580400
     $     RWORK, IWORK )                                                SLAP........580500
C                                                                        SLAP........580600
C.....THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION OF SLAP      SLAP........580700
C        WITH SUTRA.  IT ORIGINALLY READ:                                SLAP........580800
C        IF( ITER.GT.ITMAX ) IERR = 2                                    SLAP........580900
      IF (IERR.EQ.1) IERR = 2                                            SLAP........581000
      RETURN                                                             SLAP........581100
C------------- LAST LINE OF DSLUGM FOLLOWS ----------------------------  SLAP........581200
      END                                                                SLAP........581300
*DECK DSLUI                                                              SLAP........581400
      SUBROUTINE DSLUI (N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)    SLAP........581500
C***BEGIN PROLOGUE  DSLUI                                                SLAP........581600
C***PURPOSE  SLAP MSOLVE for LDU Factorization.                          SLAP........581700
C            This routine acts as an interface between the SLAP generic  SLAP........581800
C            MSOLVE calling convention and the routine that actually     SLAP........581900
C                           -1                                           SLAP........582000
C            computes  (LDU)  B = X.                                     SLAP........582100
C***LIBRARY   SLATEC (SLAP)                                              SLAP........582200
C***CATEGORY  D2E                                                        SLAP........582300
C***TYPE      DOUBLE PRECISION (SSLUI-S, DSLUI-D)                        SLAP........582400
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE, SLAP........582500
C             SLAP, SPARSE                                               SLAP........582600
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........582700
C           Seager, Mark K., (LLNL)                                      SLAP........582800
C             Lawrence Livermore National Laboratory                     SLAP........582900
C             PO BOX 808, L-60                                           SLAP........583000
C             Livermore, CA 94550 (510) 423-3141                         SLAP........583100
C             seager@llnl.gov                                            SLAP........583200
C***DESCRIPTION                                                          SLAP........583300
C       It is assumed that RWORK and IWORK have initialized with         SLAP........583400
C       the information required for DSLUI2:                             SLAP........583500
C          IWORK(1) = Starting location of IL in IWORK.                  SLAP........583600
C          IWORK(2) = Starting location of JL in IWORK.                  SLAP........583700
C          IWORK(3) = Starting location of IU in IWORK.                  SLAP........583800
C          IWORK(4) = Starting location of JU in IWORK.                  SLAP........583900
C          IWORK(5) = Starting location of L in RWORK.                   SLAP........584000
C          IWORK(6) = Starting location of DINV in RWORK.                SLAP........584100
C          IWORK(7) = Starting location of U in RWORK.                   SLAP........584200
C       See the DESCRIPTION of DSLUI2 for details.                       SLAP........584300
C***REFERENCES  (NONE)                                                   SLAP........584400
C***ROUTINES CALLED  DSLUI2                                              SLAP........584500
C***REVISION HISTORY  (YYMMDD)                                           SLAP........584600
C   871119  DATE WRITTEN                                                 SLAP........584700
C   881213  Previous REVISION DATE                                       SLAP........584800
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........584900
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........585000
C           standard.  (FNF)                                             SLAP........585100
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........585200
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........585300
C   920511  Added complete declaration section.  (WRB)                   SLAP........585400
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........585500
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........585600
C***END PROLOGUE  DSLUI                                                  SLAP........585700
C     .. Scalar Arguments ..                                             SLAP........585800
      INTEGER ISYM, N, NELT                                              SLAP........585900
C     .. Array Arguments ..                                              SLAP........586000
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)                     SLAP........586100
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........586200
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........586300
C        INTEGER IA(NELT), IWORK(10), JA(NELT)                           SLAP........586400
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........586500
C     .. Local Scalars ..                                                SLAP........586600
      INTEGER LOCDIN, LOCIL, LOCIU, LOCJL, LOCJU, LOCL, LOCU             SLAP........586700
C     .. External Subroutines ..                                         SLAP........586800
      EXTERNAL DSLUI2                                                    SLAP........586900
C***FIRST EXECUTABLE STATEMENT  DSLUI                                    SLAP........587000
C                                                                        SLAP........587100
C         Pull out the locations of the arrays holding the ILU           SLAP........587200
C         factorization.                                                 SLAP........587300
C                                                                        SLAP........587400
      LOCIL = IWORK(1)                                                   SLAP........587500
      LOCJL = IWORK(2)                                                   SLAP........587600
      LOCIU = IWORK(3)                                                   SLAP........587700
      LOCJU = IWORK(4)                                                   SLAP........587800
      LOCL = IWORK(5)                                                    SLAP........587900
      LOCDIN = IWORK(6)                                                  SLAP........588000
      LOCU = IWORK(7)                                                    SLAP........588100
C                                                                        SLAP........588200
C         Solve the system LUx = b                                       SLAP........588300
      CALL DSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL),      SLAP........588400
     $     RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) )      SLAP........588500
C                                                                        SLAP........588600
      RETURN                                                             SLAP........588700
C------------- LAST LINE OF DSLUI FOLLOWS ----------------------------   SLAP........588800
      END                                                                SLAP........588900
*DECK DSLUI2                                                             SLAP........589000
      SUBROUTINE DSLUI2 (N, B, X, IL, JL, L, DINV, IU, JU, U)            SLAP........589100
C***BEGIN PROLOGUE  DSLUI2                                               SLAP........589200
C***PURPOSE  SLAP Backsolve for LDU Factorization.                       SLAP........589300
C            Routine to solve a system of the form  L*D*U X = B,         SLAP........589400
C            where L is a unit lower triangular matrix, D is a diagonal  SLAP........589500
C            matrix, and U is a unit upper triangular matrix.            SLAP........589600
C***LIBRARY   SLATEC (SLAP)                                              SLAP........589700
C***CATEGORY  D2E                                                        SLAP........589800
C***TYPE      DOUBLE PRECISION (SSLUI2-S, DSLUI2-D)                      SLAP........589900
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM SOLVE, SLAP........590000
C             SLAP, SPARSE                                               SLAP........590100
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........590200
C           Seager, Mark K., (LLNL)                                      SLAP........590300
C             Lawrence Livermore National Laboratory                     SLAP........590400
C             PO BOX 808, L-60                                           SLAP........590500
C             Livermore, CA 94550 (510) 423-3141                         SLAP........590600
C             seager@llnl.gov                                            SLAP........590700
C***DESCRIPTION                                                          SLAP........590800
C                                                                        SLAP........590900
C *Usage:                                                                SLAP........591000
C     INTEGER N, IL(NL), JL(NL), IU(NU), JU(NU)                          SLAP........591100
C     DOUBLE PRECISION B(N), X(N), L(NL), DINV(N), U(NU)                 SLAP........591200
C                                                                        SLAP........591300
C     CALL DSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U )                 SLAP........591400
C                                                                        SLAP........591500
C *Arguments:                                                            SLAP........591600
C N      :IN       Integer                                               SLAP........591700
C         Order of the Matrix.                                           SLAP........591800
C B      :IN       Double Precision B(N).                                SLAP........591900
C         Right hand side.                                               SLAP........592000
C X      :OUT      Double Precision X(N).                                SLAP........592100
C         Solution of L*D*U x = b.                                       SLAP........592200
C IL     :IN       Integer IL(NL).                                       SLAP........592300
C JL     :IN       Integer JL(NL).                                       SLAP........592400
C L      :IN       Double Precision L(NL).                               SLAP........592500
C         IL, JL, L contain the unit  lower triangular factor of the     SLAP........592600
C         incomplete decomposition of some matrix stored in SLAP Row     SLAP........592700
C         format.  The diagonal of ones *IS* stored.  This structure     SLAP........592800
C         can   be   set  up  by   the  DSILUS  routine.   See   the     SLAP........592900
C         "Description", below  for more   details about   the  SLAP     SLAP........593000
C         format.  (NL is the number of non-zeros in the L array.)       SLAP........593100
C DINV   :IN       Double Precision DINV(N).                             SLAP........593200
C         Inverse of the diagonal matrix D.                              SLAP........593300
C IU     :IN       Integer IU(NU).                                       SLAP........593400
C JU     :IN       Integer JU(NU).                                       SLAP........593500
C U      :IN       Double Precision U(NU).                               SLAP........593600
C         IU, JU, U contain the unit upper triangular factor  of the     SLAP........593700
C         incomplete decomposition  of  some  matrix stored in  SLAP     SLAP........593800
C         Column format.   The diagonal of ones  *IS* stored.   This     SLAP........593900
C         structure can be set up  by the DSILUS routine.  See   the     SLAP........594000
C         "Description", below   for  more   details about  the SLAP     SLAP........594100
C         format.  (NU is the number of non-zeros in the U array.)       SLAP........594200
C                                                                        SLAP........594300
C *Description:                                                          SLAP........594400
C       This routine is supplied with  the SLAP package as a routine     SLAP........594500
C       to  perform  the  MSOLVE operation  in   the  SIR and   SBCG     SLAP........594600
C       iteration routines for  the  drivers DSILUR and DSLUBC.   It     SLAP........594700
C       must  be called  via   the  SLAP  MSOLVE  calling   sequence     SLAP........594800
C       convention interface routine DSLUI.                              SLAP........594900
C         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****          SLAP........595000
C               **** SLAP MSOLVE CALLING CONVENTION ****                 SLAP........595100
C                                                                        SLAP........595200
C       IL, JL, L should contain the unit lower triangular factor of     SLAP........595300
C       the incomplete decomposition of the A matrix  stored in SLAP     SLAP........595400
C       Row format.  IU, JU, U should contain  the unit upper factor     SLAP........595500
C       of the  incomplete decomposition of  the A matrix  stored in     SLAP........595600
C       SLAP Column format This ILU factorization can be computed by     SLAP........595700
C       the DSILUS routine. The diagonals (which are all one's) are      SLAP........595800
C       stored.                                                          SLAP........595900
C                                                                        SLAP........596000
C       =================== S L A P Column format ==================     SLAP........596100
C                                                                        SLAP........596200
C       This routine  requires that  the matrix A  be stored in  the     SLAP........596300
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........596400
C       counting down columns (except for  the diagonal entry, which     SLAP........596500
C       must appear first in each  "column")  and are stored  in the     SLAP........596600
C       double precision array A.   In other words,  for each column     SLAP........596700
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........596800
C       other non-zero  elements going down  the column (except  the     SLAP........596900
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........597000
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........597100
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........597200
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........597300
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........597400
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........597500
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........597600
C       number of columns in  the matrix and NELT  is the number  of     SLAP........597700
C       non-zeros in the matrix.                                         SLAP........597800
C                                                                        SLAP........597900
C       Here is an example of the  SLAP Column  storage format for a     SLAP........598000
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........598100
C       column):                                                         SLAP........598200
C                                                                        SLAP........598300
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........598400
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........598500
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........598600
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........598700
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........598800
C       | 0  0  0 44  0|                                                 SLAP........598900
C       |51  0 53  0 55|                                                 SLAP........599000
C                                                                        SLAP........599100
C       ==================== S L A P Row format ====================     SLAP........599200
C                                                                        SLAP........599300
C       This routine requires  that the matrix A  be  stored  in the     SLAP........599400
C       SLAP  Row format.   In this format  the non-zeros are stored     SLAP........599500
C       counting across  rows (except for the diagonal  entry, which     SLAP........599600
C       must  appear first  in each  "row")  and  are stored  in the     SLAP........599700
C       double precision  array A.  In other words, for each row  in     SLAP........599800
C       the matrix  put the diagonal  entry in A.   Then put in  the     SLAP........599900
C       other  non-zero elements  going across  the row  (except the     SLAP........600000
C       diagonal) in order.  The JA array holds the column index for     SLAP........600100
C       each non-zero.  The IA array holds the offsets  into the JA,     SLAP........600200
C       A  arrays  for  the   beginning  of  each  row.    That  is,     SLAP........600300
C       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-     SLAP........600400
C       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)     SLAP........600500
C       are  the last elements  of the  IROW-th row.   Note  that we     SLAP........600600
C       always have  IA(N+1) = NELT+1, where N is the number of rows     SLAP........600700
C       in the matrix  and  NELT is the  number of non-zeros  in the     SLAP........600800
C       matrix.                                                          SLAP........600900
C                                                                        SLAP........601000
C       Here is an example of the SLAP Row storage format for a  5x5     SLAP........601100
C       Matrix (in the A and JA arrays '|' denotes the end of a row):    SLAP........601200
C                                                                        SLAP........601300
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.   SLAP........601400
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........601500
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53   SLAP........601600
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........601700
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12                      SLAP........601800
C       | 0  0  0 44  0|                                                 SLAP........601900
C       |51  0 53  0 55|                                                 SLAP........602000
C                                                                        SLAP........602100
C       With  the SLAP  format  the "inner  loops" of  this  routine     SLAP........602200
C       should vectorize   on machines with   hardware  support  for     SLAP........602300
C       vector gather/scatter operations.  Your compiler may require     SLAP........602400
C       a  compiler directive  to  convince   it that there  are  no     SLAP........602500
C       implicit vector  dependencies.  Compiler directives  for the     SLAP........602600
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied     SLAP........602700
C       with the standard SLAP distribution.                             SLAP........602800
C                                                                        SLAP........602900
C***SEE ALSO  DSILUS                                                     SLAP........603000
C***REFERENCES  (NONE)                                                   SLAP........603100
C***ROUTINES CALLED  (NONE)                                              SLAP........603200
C***REVISION HISTORY  (YYMMDD)                                           SLAP........603300
C   871119  DATE WRITTEN                                                 SLAP........603400
C   881213  Previous REVISION DATE                                       SLAP........603500
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........603600
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........603700
C           standard.  (FNF)                                             SLAP........603800
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........603900
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........604000
C   920511  Added complete declaration section.  (WRB)                   SLAP........604100
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........604200
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........604300
C***END PROLOGUE  DSLUI2                                                 SLAP........604400
C     .. Scalar Arguments ..                                             SLAP........604500
      INTEGER N                                                          SLAP........604600
C     .. Array Arguments ..                                              SLAP........604700
      DOUBLE PRECISION B(N), DINV(N), L(*), U(*), X(N)                   SLAP........604800
      INTEGER IL(*), IU(*), JL(*), JU(*)                                 SLAP........604900
C     .. Local Scalars ..                                                SLAP........605000
      INTEGER I, ICOL, IROW, J, JBGN, JEND                               SLAP........605100
C***FIRST EXECUTABLE STATEMENT  DSLUI2                                   SLAP........605200
C                                                                        SLAP........605300
C         Solve  L*Y = B,  storing result in X, L stored by rows.        SLAP........605400
C                                                                        SLAP........605500
      DO 10 I = 1, N                                                     SLAP........605600
         X(I) = B(I)                                                     SLAP........605700
 10   CONTINUE                                                           SLAP........605800
      DO 30 IROW = 2, N                                                  SLAP........605900
         JBGN = IL(IROW)                                                 SLAP........606000
         JEND = IL(IROW+1)-1                                             SLAP........606100
         IF( JBGN.LE.JEND ) THEN                                         SLAP........606200
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........606300
CDIR$ IVDEP                                                              SLAP........606400
CVD$ ASSOC                                                               SLAP........606500
CVD$ NODEPCHK                                                            SLAP........606600
            DO 20 J = JBGN, JEND                                         SLAP........606700
               X(IROW) = X(IROW) - L(J)*X(JL(J))                         SLAP........606800
 20         CONTINUE                                                     SLAP........606900
         ENDIF                                                           SLAP........607000
 30   CONTINUE                                                           SLAP........607100
C                                                                        SLAP........607200
C         Solve  D*Z = Y,  storing result in X.                          SLAP........607300
      DO 40 I=1,N                                                        SLAP........607400
         X(I) = X(I)*DINV(I)                                             SLAP........607500
 40   CONTINUE                                                           SLAP........607600
C                                                                        SLAP........607700
C         Solve  U*X = Z, U stored by columns.                           SLAP........607800
      DO 60 ICOL = N, 2, -1                                              SLAP........607900
         JBGN = JU(ICOL)                                                 SLAP........608000
         JEND = JU(ICOL+1)-1                                             SLAP........608100
         IF( JBGN.LE.JEND ) THEN                                         SLAP........608200
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........608300
CDIR$ IVDEP                                                              SLAP........608400
CVD$ NODEPCHK                                                            SLAP........608500
            DO 50 J = JBGN, JEND                                         SLAP........608600
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)                        SLAP........608700
 50         CONTINUE                                                     SLAP........608800
         ENDIF                                                           SLAP........608900
 60   CONTINUE                                                           SLAP........609000
C                                                                        SLAP........609100
      RETURN                                                             SLAP........609200
C------------- LAST LINE OF DSLUI2 FOLLOWS ----------------------------  SLAP........609300
      END                                                                SLAP........609400
*DECK DSLUOM                                                             SLAP........609500
      SUBROUTINE DSLUOM (N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL,    SLAP........609600
     +   TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)  SLAP........609700
C***BEGIN PROLOGUE  DSLUOM                                               SLAP........609800
C***PURPOSE  Incomplete LU Orthomin Sparse Iterative Ax=b Solver.        SLAP........609900
C            Routine to solve a general linear system  Ax = b  using     SLAP........610000
C            the Orthomin method with Incomplete LU decomposition.       SLAP........610100
C***LIBRARY   SLATEC (SLAP)                                              SLAP........610200
C***CATEGORY  D2A4, D2B4                                                 SLAP........610300
C***TYPE      DOUBLE PRECISION (SSLUOM-S, DSLUOM-D)                      SLAP........610400
C***KEYWORDS  ITERATIVE INCOMPLETE LU PRECONDITION,                      SLAP........610500
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........610600
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........610700
C           Seager, Mark K., (LLNL)                                      SLAP........610800
C             Lawrence Livermore National Laboratory                     SLAP........610900
C             PO BOX 808, L-60                                           SLAP........611000
C             Livermore, CA 94550 (510) 423-3141                         SLAP........611100
C             seager@llnl.gov                                            SLAP........611200
C***DESCRIPTION                                                          SLAP........611300
C                                                                        SLAP........611400
C *Usage:                                                                SLAP........611500
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX      SLAP........611600
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NL+NU+4*N+2), LENIW         SLAP........611700
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR                     SLAP........611800
C     DOUBLE PRECISION RWORK(NL+NU+7*N+3*N*NSAVE+NSAVE)                  SLAP........611900
C                                                                        SLAP........612000
C     CALL DSLUOM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, ITOL, TOL,      SLAP........612100
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )    SLAP........612200
C                                                                        SLAP........612300
C *Arguments:                                                            SLAP........612400
C N      :IN       Integer.                                              SLAP........612500
C         Order of the matrix.                                           SLAP........612600
C B      :IN       Double Precision B(N).                                SLAP........612700
C         Right-hand side vector.                                        SLAP........612800
C X      :INOUT    Double Precision X(N).                                SLAP........612900
C         On input X is your initial guess for solution vector.          SLAP........613000
C         On output X is the final approximate solution.                 SLAP........613100
C NELT   :IN       Integer.                                              SLAP........613200
C         Number of Non-Zeros stored in A.                               SLAP........613300
C IA     :INOUT    Integer IA(NELT).                                     SLAP........613400
C JA     :INOUT    Integer JA(NELT).                                     SLAP........613500
C A      :INOUT    Double Precision A(NELT).                             SLAP........613600
C         These arrays should hold the matrix A in either the SLAP       SLAP........613700
C         Triad format or the SLAP Column format.  See "Description",    SLAP........613800
C         below.  If the SLAP Triad format is chosen, it is changed      SLAP........613900
C         internally to the SLAP Column format.                          SLAP........614000
C ISYM   :IN       Integer.                                              SLAP........614100
C         Flag to indicate symmetric storage format.                     SLAP........614200
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........614300
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........614400
C         or lower triangle of the matrix is stored.                     SLAP........614500
C NSAVE  :IN       Integer.                                              SLAP........614600
C         Number of direction vectors to save and orthogonalize against. SLAP........614700
C ITOL   :IN       Integer.                                              SLAP........614800
C         Flag to indicate type of convergence criterion.                SLAP........614900
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........615000
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........615100
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........615200
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........615300
C         side is less than TOL, where M-inv is the inverse of the       SLAP........615400
C         diagonal of A.                                                 SLAP........615500
C         ITOL=11 is often useful for checking and comparing different   SLAP........615600
C         routines.  For this case, the user must supply the "exact"     SLAP........615700
C         solution or a very accurate approximation (one with an error   SLAP........615800
C         much less than TOL) through a common block,                    SLAP........615900
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........616000
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........616100
C         between the iterative approximation and the user-supplied      SLAP........616200
C         solution divided by the 2-norm of the user-supplied solution   SLAP........616300
C         is less than TOL.  Note that this requires the user to set up  SLAP........616400
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........616500
C         The routine with this declaration should be loaded before the  SLAP........616600
C         stop test so that the correct length is used by the loader.    SLAP........616700
C         This procedure is not standard Fortran and may not work        SLAP........616800
C         correctly on your system (although it has worked on every      SLAP........616900
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........617000
C         common block is indeed standard Fortran.                       SLAP........617100
C TOL    :INOUT    Double Precision.                                     SLAP........617200
C         Convergence criterion, as described above.  (Reset if IERR=4.) SLAP........617300
C ITMAX  :IN       Integer.                                              SLAP........617400
C         Maximum number of iterations.                                  SLAP........617500
C ITER   :OUT      Integer.                                              SLAP........617600
C         Number of iterations required to reach convergence, or         SLAP........617700
C         ITMAX+1 if convergence criterion could not be achieved in      SLAP........617800
C         ITMAX iterations.                                              SLAP........617900
C ERR    :OUT      Double Precision.                                     SLAP........618000
C         Error estimate of error in final approximate solution, as      SLAP........618100
C         defined by ITOL.                                               SLAP........618200
C IERR   :OUT      Integer.                                              SLAP........618300
C         Return error flag.                                             SLAP........618400
C           IERR = 0 => All went well.                                   SLAP........618500
C           IERR = 1 => Insufficient space allocated for WORK or IWORK.  SLAP........618600
C           IERR = 2 => Method failed to converge in ITMAX steps.        SLAP........618700
C           IERR = 3 => Error in user input.                             SLAP........618800
C                       Check input values of N, ITOL.                   SLAP........618900
C           IERR = 4 => User error tolerance set too tight.              SLAP........619000
C                       Reset to 500*D1MACH(3).  Iteration proceeded.    SLAP........619100
C           IERR = 5 => Preconditioning matrix, M, is not positive       SLAP........619200
C                       definite.  (r,z) < 0.                            SLAP........619300
C           IERR = 6 => Breakdown of the method detected.                SLAP........619400
C                       (p,Ap) < epsilon**2.                             SLAP........619500
C           IERR = 7 => Incomplete factorization broke down and was      SLAP........619600
C                       fudged.  Resulting preconditioning may be less   SLAP........619700
C                       than the best.                                   SLAP........619800
C IUNIT  :IN       Integer.                                              SLAP........619900
C         Unit number on which to write the error at each iteration,     SLAP........620000
C         if this is desired for monitoring convergence.  If unit        SLAP........620100
C         number is 0, no writing will occur.                            SLAP........620200
C RWORK  :WORK     Double Precision RWORK(LENW).                         SLAP........620300
C         Double Precision array used for workspace.  NL is the number   SLAP........620400
C         of non-zeros in the lower triangle of the matrix (including    SLAP........620500
C         the diagonal).  NU is the number of non-zeros in the upper     SLAP........620600
C         triangle of the matrix (including the diagonal).               SLAP........620700
C LENW   :IN       Integer.                                              SLAP........620800
C         Length of the double precision workspace, RWORK.               SLAP........620900
C         LENW >= NL+NU+4*N+NSAVE*(3*N+1)                                SLAP........621000
C IWORK  :WORK     Integer IWORK(LENIW)                                  SLAP........621100
C         Integer array used for workspace.  NL is the number of non-    SLAP........621200
C         zeros in the lower triangle of the matrix (including the       SLAP........621300
C         diagonal).  NU is the number of non-zeros in the upper         SLAP........621400
C         triangle of the matrix (including the diagonal).               SLAP........621500
C         Upon return the following locations of IWORK hold information  SLAP........621600
C         which may be of use to the user:                               SLAP........621700
C         IWORK(9)  Amount of Integer workspace actually used.           SLAP........621800
C         IWORK(10) Amount of Double Precision workspace actually used.  SLAP........621900
C LENIW  :IN       Integer.                                              SLAP........622000
C         Length of the integer workspace, IWORK.                        SLAP........622100
C         LENIW >= NL+NU+4*N+12.                                         SLAP........622200
C                                                                        SLAP........622300
C *Description:                                                          SLAP........622400
C       This routine is  simply a driver  for  the DOMN routine.  It     SLAP........622500
C       calls the DSILUS routine  to set  up the preconditioning and     SLAP........622600
C       then  calls   DOMN  with the appropriate  MATVEC  and MSOLVE     SLAP........622700
C       routines.                                                        SLAP........622800
C                                                                        SLAP........622900
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix     SLAP........623000
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP     SLAP........623100
C       Column format.  The user can hand this routine either of the     SLAP........623200
C       of these data structures and SLAP  will figure out  which on     SLAP........623300
C       is being used and act accordingly.                               SLAP........623400
C                                                                        SLAP........623500
C       =================== S L A P Triad format ===================     SLAP........623600
C                                                                        SLAP........623700
C       This routine requires that the  matrix A be   stored in  the     SLAP........623800
C       SLAP  Triad format.  In  this format only the non-zeros  are     SLAP........623900
C       stored.  They may appear in  *ANY* order.  The user supplies     SLAP........624000
C       three arrays of  length NELT, where  NELT is  the number  of     SLAP........624100
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For     SLAP........624200
C       each non-zero the user puts the row and column index of that     SLAP........624300
C       matrix element  in the IA and  JA arrays.  The  value of the     SLAP........624400
C       non-zero   matrix  element is  placed  in  the corresponding     SLAP........624500
C       location of the A array.   This is  an  extremely  easy data     SLAP........624600
C       structure to generate.  On  the  other hand it   is  not too     SLAP........624700
C       efficient on vector computers for  the iterative solution of     SLAP........624800
C       linear systems.  Hence,   SLAP changes   this  input    data     SLAP........624900
C       structure to the SLAP Column format  for  the iteration (but     SLAP........625000
C       does not change it back).                                        SLAP........625100
C                                                                        SLAP........625200
C       Here is an example of the  SLAP Triad   storage format for a     SLAP........625300
C       5x5 Matrix.  Recall that the entries may appear in any order.    SLAP........625400
C                                                                        SLAP........625500
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.    SLAP........625600
C                              1  2  3  4  5  6  7  8  9 10 11           SLAP........625700
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21           SLAP........625800
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2           SLAP........625900
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1           SLAP........626000
C       | 0  0  0 44  0|                                                 SLAP........626100
C       |51  0 53  0 55|                                                 SLAP........626200
C                                                                        SLAP........626300
C       =================== S L A P Column format ==================     SLAP........626400
C                                                                        SLAP........626500
C       This routine  requires that  the matrix A  be stored in  the     SLAP........626600
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........626700
C       counting down columns (except for  the diagonal entry, which     SLAP........626800
C       must appear first in each  "column")  and are stored  in the     SLAP........626900
C       double precision array A.   In other words,  for each column     SLAP........627000
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........627100
C       other non-zero  elements going down  the column (except  the     SLAP........627200
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........627300
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........627400
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........627500
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........627600
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........627700
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........627800
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........627900
C       number of columns in  the matrix and NELT  is the number  of     SLAP........628000
C       non-zeros in the matrix.                                         SLAP........628100
C                                                                        SLAP........628200
C       Here is an example of the  SLAP Column  storage format for a     SLAP........628300
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........628400
C       column):                                                         SLAP........628500
C                                                                        SLAP........628600
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........628700
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........628800
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........628900
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........629000
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........629100
C       | 0  0  0 44  0|                                                 SLAP........629200
C       |51  0 53  0 55|                                                 SLAP........629300
C                                                                        SLAP........629400
C *Side Effects:                                                         SLAP........629500
C       The SLAP Triad format (IA, JA,  A) is modified internally to     SLAP........629600
C       be the SLAP Column format.  See above.                           SLAP........629700
C                                                                        SLAP........629800
C *Cautions:                                                             SLAP........629900
C     This routine will attempt to write to the Fortran logical output   SLAP........630000
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........630100
C     this logical unit is attached to a file or terminal before calling SLAP........630200
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........630300
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........630400
C                                                                        SLAP........630500
C***SEE ALSO  DOMN, DSDOMN                                               SLAP........630600
C***REFERENCES  (NONE)                                                   SLAP........630700
C***ROUTINES CALLED  DCHKW, DOMN, DS2Y, DSILUS, DSLUI, DSMV              SLAP........630800
C***REVISION HISTORY  (YYMMDD)                                           SLAP........630900
C   890404  DATE WRITTEN                                                 SLAP........631000
C   890404  Previous REVISION DATE                                       SLAP........631100
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........631200
C   890921  Removed TeX from comments.  (FNF)                            SLAP........631300
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........631400
C           standard.  (FNF)                                             SLAP........631500
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........631600
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........631700
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........631800
C   920511  Added complete declaration section.  (WRB)                   SLAP........631900
C   921019  Corrected NEL to NL.  (FNF)                                  SLAP........632000
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........632100
C***END PROLOGUE  DSLUOM                                                 SLAP........632200
C     .. Parameters ..                                                   SLAP........632300
      INTEGER LOCRB, LOCIB                                               SLAP........632400
      PARAMETER (LOCRB=1, LOCIB=11)                                      SLAP........632500
C     .. Scalar Arguments ..                                             SLAP........632600
      DOUBLE PRECISION ERR, TOL                                          SLAP........632700
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, LENIW, LENW, N,      SLAP........632800
     +        NELT, NSAVE                                                SLAP........632900
C     .. Array Arguments ..                                              SLAP........633000
      DOUBLE PRECISION A(N), B(N), RWORK(LENW), X(N)                     SLAP........633100
      INTEGER IA(NELT), IWORK(LENIW), JA(NELT)                           SLAP........633200
C     .. Local Scalars ..                                                SLAP........633300
      INTEGER ICOL, J, JBGN, JEND, LOCAP, LOCCSA, LOCDIN, LOCDZ, LOCEMA, SLAP........633400
     +        LOCIL, LOCIU, LOCIW, LOCJL, LOCJU, LOCL, LOCNC, LOCNR,     SLAP........633500
     +        LOCP, LOCR, LOCU, LOCW, LOCZ, NL, NU                       SLAP........633600
C     .. External Subroutines ..                                         SLAP........633700
      EXTERNAL DCHKW, DOMN, DS2Y, DSILUS, DSLUI, DSMV                    SLAP........633800
C***FIRST EXECUTABLE STATEMENT  DSLUOM                                   SLAP........633900
C                                                                        SLAP........634000
      IERR = 0                                                           SLAP........634100
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN                                   SLAP........634200
         IERR = 3                                                        SLAP........634300
         RETURN                                                          SLAP........634400
      ENDIF                                                              SLAP........634500
C                                                                        SLAP........634600
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.  SLAP........634700
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )                              SLAP........634800
C                                                                        SLAP........634900
C         Count number of Non-Zero elements preconditioner ILU matrix.   SLAP........635000
C         Then set up the work arrays.                                   SLAP........635100
      NL = 0                                                             SLAP........635200
      NU = 0                                                             SLAP........635300
      DO 20 ICOL = 1, N                                                  SLAP........635400
C         Don't count diagonal.                                          SLAP........635500
         JBGN = JA(ICOL)+1                                               SLAP........635600
         JEND = JA(ICOL+1)-1                                             SLAP........635700
         IF( JBGN.LE.JEND ) THEN                                         SLAP........635800
CVD$ NOVECTOR                                                            SLAP........635900
            DO 10 J = JBGN, JEND                                         SLAP........636000
               IF( IA(J).GT.ICOL ) THEN                                  SLAP........636100
                  NL = NL + 1                                            SLAP........636200
                  IF( ISYM.NE.0 ) NU = NU + 1                            SLAP........636300
               ELSE                                                      SLAP........636400
                  NU = NU + 1                                            SLAP........636500
               ENDIF                                                     SLAP........636600
 10         CONTINUE                                                     SLAP........636700
         ENDIF                                                           SLAP........636800
 20   CONTINUE                                                           SLAP........636900
C                                                                        SLAP........637000
      LOCIL = LOCIB                                                      SLAP........637100
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........637200
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........637300
C        LOCJL = LOCIL + N+1                                             SLAP........637400
      LOCJL = LOCIL + NL                                                 SLAP........637500
      LOCIU = LOCJL + NL                                                 SLAP........637600
      LOCJU = LOCIU + NU                                                 SLAP........637700
C.....THE NEXT LINE OF CODE INCORPORATES A BUG FIX MADE DURING           SLAP........637800
C        INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:            SLAP........637900
C        LOCNR = LOCJU + N+1                                             SLAP........638000
      LOCNR = LOCJU + NU                                                 SLAP........638100
      LOCNC = LOCNR + N                                                  SLAP........638200
      LOCIW = LOCNC + N                                                  SLAP........638300
C                                                                        SLAP........638400
      LOCL   = LOCRB                                                     SLAP........638500
      LOCDIN = LOCL + NL                                                 SLAP........638600
      LOCU   = LOCDIN + N                                                SLAP........638700
      LOCR   = LOCU + NU                                                 SLAP........638800
      LOCZ   = LOCR + N                                                  SLAP........638900
      LOCP   = LOCZ + N                                                  SLAP........639000
      LOCAP  = LOCP + N*(NSAVE+1)                                        SLAP........639100
      LOCEMA = LOCAP + N*(NSAVE+1)                                       SLAP........639200
      LOCDZ  = LOCEMA + N*(NSAVE+1)                                      SLAP........639300
      LOCCSA = LOCDZ + N                                                 SLAP........639400
      LOCW   = LOCCSA + NSAVE                                            SLAP........639500
C                                                                        SLAP........639600
C         Check the workspace allocations.                               SLAP........639700
      CALL DCHKW( 'DSLUOM', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )  SLAP........639800
      IF( IERR.NE.0 ) RETURN                                             SLAP........639900
C                                                                        SLAP........640000
      IWORK(1) = LOCIL                                                   SLAP........640100
      IWORK(2) = LOCJL                                                   SLAP........640200
      IWORK(3) = LOCIU                                                   SLAP........640300
      IWORK(4) = LOCJU                                                   SLAP........640400
      IWORK(5) = LOCL                                                    SLAP........640500
      IWORK(6) = LOCDIN                                                  SLAP........640600
      IWORK(7) = LOCU                                                    SLAP........640700
      IWORK(9) = LOCIW                                                   SLAP........640800
      IWORK(10) = LOCW                                                   SLAP........640900
C                                                                        SLAP........641000
C         Compute the Incomplete LU decomposition.                       SLAP........641100
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),           SLAP........641200
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),   SLAP........641300
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )       SLAP........641400
C                                                                        SLAP........641500
C         Perform the incomplete LU preconditioned OrthoMin algorithm.   SLAP........641600
      CALL DOMN(N, B, X, NELT, IA, JA, A, ISYM, DSMV,                    SLAP........641700
     $     DSLUI, NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,       SLAP........641800
     $     RWORK(LOCR), RWORK(LOCZ), RWORK(LOCP), RWORK(LOCAP),          SLAP........641900
     $     RWORK(LOCEMA), RWORK(LOCDZ), RWORK(LOCCSA),                   SLAP........642000
     $     RWORK, IWORK )                                                SLAP........642100
      RETURN                                                             SLAP........642200
      END                                                                SLAP........642300
*DECK DSMV                                                               SLAP........642400
      SUBROUTINE DSMV (N, X, Y, NELT, IA, JA, A, ISYM)                   SLAP........642500
C***BEGIN PROLOGUE  DSMV                                                 SLAP........642600
C***PURPOSE  SLAP Column Format Sparse Matrix Vector Product.            SLAP........642700
C            Routine to calculate the sparse matrix vector product:      SLAP........642800
C            Y = A*X.                                                    SLAP........642900
C***LIBRARY   SLATEC (SLAP)                                              SLAP........643000
C***CATEGORY  D1B4                                                       SLAP........643100
C***TYPE      DOUBLE PRECISION (SSMV-S, DSMV-D)                          SLAP........643200
C***KEYWORDS  MATRIX VECTOR MULTIPLY, SLAP, SPARSE                       SLAP........643300
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........643400
C           Seager, Mark K., (LLNL)                                      SLAP........643500
C             Lawrence Livermore National Laboratory                     SLAP........643600
C             PO BOX 808, L-60                                           SLAP........643700
C             Livermore, CA 94550 (510) 423-3141                         SLAP........643800
C             seager@llnl.gov                                            SLAP........643900
C***DESCRIPTION                                                          SLAP........644000
C                                                                        SLAP........644100
C *Usage:                                                                SLAP........644200
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM                         SLAP........644300
C     DOUBLE PRECISION X(N), Y(N), A(NELT)                               SLAP........644400
C                                                                        SLAP........644500
C     CALL DSMV(N, X, Y, NELT, IA, JA, A, ISYM )                         SLAP........644600
C                                                                        SLAP........644700
C *Arguments:                                                            SLAP........644800
C N      :IN       Integer.                                              SLAP........644900
C         Order of the Matrix.                                           SLAP........645000
C X      :IN       Double Precision X(N).                                SLAP........645100
C         The vector that should be multiplied by the matrix.            SLAP........645200
C Y      :OUT      Double Precision Y(N).                                SLAP........645300
C         The product of the matrix and the vector.                      SLAP........645400
C NELT   :IN       Integer.                                              SLAP........645500
C         Number of Non-Zeros stored in A.                               SLAP........645600
C IA     :IN       Integer IA(NELT).                                     SLAP........645700
C JA     :IN       Integer JA(NELT).                                     SLAP........645800
C A      :IN       Double Precision A(NELT).                             SLAP........645900
C         These arrays should hold the matrix A in the SLAP Column       SLAP........646000
C         format.  See "Description", below.                             SLAP........646100
C ISYM   :IN       Integer.                                              SLAP........646200
C         Flag to indicate symmetric storage format.                     SLAP........646300
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........646400
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........646500
C         or lower triangle of the matrix is stored.                     SLAP........646600
C                                                                        SLAP........646700
C *Description                                                           SLAP........646800
C       =================== S L A P Column format ==================     SLAP........646900
C       This routine  requires that  the matrix A  be stored in  the     SLAP........647000
C       SLAP Column format.  In this format the non-zeros are stored     SLAP........647100
C       counting down columns (except for  the diagonal entry, which     SLAP........647200
C       must appear first in each  "column")  and are stored  in the     SLAP........647300
C       double precision array A.   In other words,  for each column     SLAP........647400
C       in the matrix put the diagonal entry in  A.  Then put in the     SLAP........647500
C       other non-zero  elements going down  the column (except  the     SLAP........647600
C       diagonal) in order.   The  IA array holds the  row index for     SLAP........647700
C       each non-zero.  The JA array holds the offsets  into the IA,     SLAP........647800
C       A arrays  for  the  beginning  of each   column.   That  is,     SLAP........647900
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the     SLAP........648000
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),     SLAP........648100
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.     SLAP........648200
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the     SLAP........648300
C       number of columns in  the matrix and NELT  is the number  of     SLAP........648400
C       non-zeros in the matrix.                                         SLAP........648500
C                                                                        SLAP........648600
C       Here is an example of the  SLAP Column  storage format for a     SLAP........648700
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a     SLAP........648800
C       column):                                                         SLAP........648900
C                                                                        SLAP........649000
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.   SLAP........649100
C                              1  2  3    4  5    6  7    8    9 10 11   SLAP........649200
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35   SLAP........649300
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3   SLAP........649400
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12                      SLAP........649500
C       | 0  0  0 44  0|                                                 SLAP........649600
C       |51  0 53  0 55|                                                 SLAP........649700
C                                                                        SLAP........649800
C       With  the SLAP  format  the "inner  loops" of  this  routine     SLAP........649900
C       should vectorize   on machines with   hardware  support  for     SLAP........650000
C       vector gather/scatter operations.  Your compiler may require     SLAP........650100
C       a  compiler directive  to  convince   it that there  are  no     SLAP........650200
C       implicit vector  dependencies.  Compiler directives  for the     SLAP........650300
C       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied     SLAP........650400
C       with the standard SLAP distribution.                             SLAP........650500
C                                                                        SLAP........650600
C *Cautions:                                                             SLAP........650700
C     This   routine   assumes  that  the matrix A is stored in SLAP     SLAP........650800
C     Column format.  It does not check  for  this (for  speed)  and     SLAP........650900
C     evil, ugly, ornery and nasty things  will happen if the matrix     SLAP........651000
C     data  structure  is,  in fact, not SLAP Column.  Beware of the     SLAP........651100
C     wrong data structure!!!                                            SLAP........651200
C                                                                        SLAP........651300
C***SEE ALSO  DSMTV                                                      SLAP........651400
C***REFERENCES  (NONE)                                                   SLAP........651500
C***ROUTINES CALLED  (NONE)                                              SLAP........651600
C***REVISION HISTORY  (YYMMDD)                                           SLAP........651700
C   871119  DATE WRITTEN                                                 SLAP........651800
C   881213  Previous REVISION DATE                                       SLAP........651900
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........652000
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........652100
C           standard.  (FNF)                                             SLAP........652200
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........652300
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........652400
C   920511  Added complete declaration section.  (WRB)                   SLAP........652500
C   930701  Updated CATEGORY section.  (FNF, WRB)                        SLAP........652600
C***END PROLOGUE  DSMV                                                   SLAP........652700
C     .. Scalar Arguments ..                                             SLAP........652800
      INTEGER ISYM, N, NELT                                              SLAP........652900
C     .. Array Arguments ..                                              SLAP........653000
      DOUBLE PRECISION A(NELT), X(N), Y(N)                               SLAP........653100
      INTEGER IA(NELT), JA(NELT)                                         SLAP........653200
C     .. Local Scalars ..                                                SLAP........653300
      INTEGER I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND                   SLAP........653400
C***FIRST EXECUTABLE STATEMENT  DSMV                                     SLAP........653500
C                                                                        SLAP........653600
C         Zero out the result vector.                                    SLAP........653700
C                                                                        SLAP........653800
      DO 10 I = 1, N                                                     SLAP........653900
         Y(I) = 0                                                        SLAP........654000
 10   CONTINUE                                                           SLAP........654100
C                                                                        SLAP........654200
C         Multiply by A.                                                 SLAP........654300
C                                                                        SLAP........654400
CVD$R NOCONCUR                                                           SLAP........654500
      DO 30 ICOL = 1, N                                                  SLAP........654600
         IBGN = JA(ICOL)                                                 SLAP........654700
         IEND = JA(ICOL+1)-1                                             SLAP........654800
CLLL. OPTION ASSERT (NOHAZARD)                                           SLAP........654900
CDIR$ IVDEP                                                              SLAP........655000
CVD$ NODEPCHK                                                            SLAP........655100
         DO 20 I = IBGN, IEND                                            SLAP........655200
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)                           SLAP........655300
 20      CONTINUE                                                        SLAP........655400
 30   CONTINUE                                                           SLAP........655500
C                                                                        SLAP........655600
      IF( ISYM.EQ.1 ) THEN                                               SLAP........655700
C                                                                        SLAP........655800
C         The matrix is non-symmetric.  Need to get the other half in... SLAP........655900
C         This loops assumes that the diagonal is the first entry in     SLAP........656000
C         each column.                                                   SLAP........656100
C                                                                        SLAP........656200
         DO 50 IROW = 1, N                                               SLAP........656300
            JBGN = JA(IROW)+1                                            SLAP........656400
            JEND = JA(IROW+1)-1                                          SLAP........656500
            IF( JBGN.GT.JEND ) GOTO 50                                   SLAP........656600
            DO 40 J = JBGN, JEND                                         SLAP........656700
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))                         SLAP........656800
 40         CONTINUE                                                     SLAP........656900
 50      CONTINUE                                                        SLAP........657000
      ENDIF                                                              SLAP........657100
      RETURN                                                             SLAP........657200
C------------- LAST LINE OF DSMV FOLLOWS ----------------------------    SLAP........657300
      END                                                                SLAP........657400
*DECK DXLCAL                                                             SLAP........657500
      SUBROUTINE DXLCAL (N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,   SLAP........657600
     +   WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, NELT, IA, JA, A, SLAP........657700
     +   ISYM)                                                           SLAP........657800
C***BEGIN PROLOGUE  DXLCAL                                               SLAP........657900
C***SUBSIDIARY                                                           SLAP........658000
C***PURPOSE  Internal routine for DGMRES.                                SLAP........658100
C***LIBRARY   SLATEC (SLAP)                                              SLAP........658200
C***CATEGORY  D2A4, D2B4                                                 SLAP........658300
C***TYPE      DOUBLE PRECISION (SXLCAL-S, DXLCAL-D)                      SLAP........658400
C***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,      SLAP........658500
C             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE                  SLAP........658600
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........658700
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........658800
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........658900
C             Lawrence Livermore National Laboratory                     SLAP........659000
C             PO Box 808, L-60                                           SLAP........659100
C             Livermore, CA 94550 (510) 423-3141                         SLAP........659200
C***DESCRIPTION                                                          SLAP........659300
C        This  routine computes the solution  XL,  the current DGMRES    SLAP........659400
C        iterate, given the  V(I)'s and  the  QR factorization of the    SLAP........659500
C        Hessenberg  matrix HES.   This routine  is  only called when    SLAP........659600
C        ITOL=11.                                                        SLAP........659700
C                                                                        SLAP........659800
C *Usage:                                                                SLAP........659900
C      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED)    SLAP........660000
C      INTEGER NELT, IA(NELT), JA(NELT), ISYM                            SLAP........660100
C      DOUBLE PRECISION X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL), SLAP........660200
C     $                 V(N,MAXLP1), R0NRM, WK(N), SZ(N),                SLAP........660300
C     $                 RPAR(USER DEFINED), A(NELT)                      SLAP........660400
C      EXTERNAL MSOLVE                                                   SLAP........660500
C                                                                        SLAP........660600
C      CALL DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM,         SLAP........660700
C     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR,               SLAP........660800
C     $     NELT, IA, JA, A, ISYM)                                       SLAP........660900
C                                                                        SLAP........661000
C *Arguments:                                                            SLAP........661100
C N      :IN       Integer                                               SLAP........661200
C         The order of the matrix A, and the lengths                     SLAP........661300
C         of the vectors SR, SZ, R0 and Z.                               SLAP........661400
C LGMR   :IN       Integer                                               SLAP........661500
C         The number of iterations performed and                         SLAP........661600
C         the current order of the upper Hessenberg                      SLAP........661700
C         matrix HES.                                                    SLAP........661800
C X      :IN       Double Precision X(N)                                 SLAP........661900
C         The current approximate solution as of the last restart.       SLAP........662000
C XL     :OUT      Double Precision XL(N)                                SLAP........662100
C         An array of length N used to hold the approximate              SLAP........662200
C         solution X(L).                                                 SLAP........662300
C         Warning: XL and ZL are the same array in the calling routine.  SLAP........662400
C ZL     :IN       Double Precision ZL(N)                                SLAP........662500
C         An array of length N used to hold the approximate              SLAP........662600
C         solution Z(L).                                                 SLAP........662700
C HES    :IN       Double Precision HES(MAXLP1,MAXL)                     SLAP........662800
C         The upper triangular factor of the QR decomposition            SLAP........662900
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose          SLAP........663000
C         entries are the scaled inner-products of A*V(*,i) and V(*,k).  SLAP........663100
C MAXLP1 :IN       Integer                                               SLAP........663200
C         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.       SLAP........663300
C         MAXL is the maximum allowable order of the matrix HES.         SLAP........663400
C Q      :IN       Double Precision Q(2*MAXL)                            SLAP........663500
C         A double precision array of length 2*MAXL containing the       SLAP........663600
C         components of the Givens rotations used in the QR              SLAP........663700
C         decomposition of HES.  It is loaded in DHEQR.                  SLAP........663800
C V      :IN       Double Precision V(N,MAXLP1)                          SLAP........663900
C         The N by(LGMR+1) array containing the LGMR                     SLAP........664000
C         orthogonal vectors V(*,1) to V(*,LGMR).                        SLAP........664100
C R0NRM  :IN       Double Precision                                      SLAP........664200
C         The scaled norm of the initial residual for the                SLAP........664300
C         current call to DPIGMR.                                        SLAP........664400
C WK     :IN       Double Precision WK(N)                                SLAP........664500
C         A double precision work array of length N.                     SLAP........664600
C SZ     :IN       Double Precision SZ(N)                                SLAP........664700
C         A vector of length N containing the non-zero                   SLAP........664800
C         elements of the diagonal scaling matrix for Z.                 SLAP........664900
C JSCAL  :IN       Integer                                               SLAP........665000
C         A flag indicating whether arrays SR and SZ are used.           SLAP........665100
C         JSCAL=0 means SR and SZ are not used and the                   SLAP........665200
C                 algorithm will perform as if all                       SLAP........665300
C                 SR(i) = 1 and SZ(i) = 1.                               SLAP........665400
C         JSCAL=1 means only SZ is used, and the algorithm               SLAP........665500
C                 performs as if all SR(i) = 1.                          SLAP........665600
C         JSCAL=2 means only SR is used, and the algorithm               SLAP........665700
C                 performs as if all SZ(i) = 1.                          SLAP........665800
C         JSCAL=3 means both SR and SZ are used.                         SLAP........665900
C JPRE   :IN       Integer                                               SLAP........666000
C         The preconditioner type flag.                                  SLAP........666100
C MSOLVE :EXT      External.                                             SLAP........666200
C         Name of the routine which solves a linear system Mz = r for    SLAP........666300
C         z given r with the preconditioning matrix M (M is supplied via SLAP........666400
C         RPAR and IPAR arrays.  The name of the MSOLVE routine must     SLAP........666500
C         be declared external in the calling program.  The calling      SLAP........666600
C         sequence to MSOLVE is:                                         SLAP........666700
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)    SLAP........666800
C         Where N is the number of unknowns, R is the right-hand side    SLAP........666900
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........667000
C         ISYM are defined as below.  RPAR is a double precision array   SLAP........667100
C         that can be used to pass necessary preconditioning information SLAP........667200
C         and/or workspace to MSOLVE.  IPAR is an integer work array     SLAP........667300
C         for the same purpose as RPAR.                                  SLAP........667400
C NMSL   :IN       Integer                                               SLAP........667500
C         The number of calls to MSOLVE.                                 SLAP........667600
C RPAR   :IN       Double Precision RPAR(USER DEFINED)                   SLAP........667700
C         Double Precision workspace passed directly to the MSOLVE       SLAP........667800
C         routine.                                                       SLAP........667900
C IPAR   :IN       Integer IPAR(USER DEFINED)                            SLAP........668000
C         Integer workspace passed directly to the MSOLVE routine.       SLAP........668100
C NELT   :IN       Integer                                               SLAP........668200
C         The length of arrays IA, JA and A.                             SLAP........668300
C IA     :IN       Integer IA(NELT)                                      SLAP........668400
C         An integer array of length NELT containing matrix data.        SLAP........668500
C         It is passed directly to the MATVEC and MSOLVE routines.       SLAP........668600
C JA     :IN       Integer JA(NELT)                                      SLAP........668700
C         An integer array of length NELT containing matrix data.        SLAP........668800
C         It is passed directly to the MATVEC and MSOLVE routines.       SLAP........668900
C A      :IN       Double Precision A(NELT)                              SLAP........669000
C         A double precision array of length NELT containing matrix      SLAP........669100
C         data.                                                          SLAP........669200
C         It is passed directly to the MATVEC and MSOLVE routines.       SLAP........669300
C ISYM   :IN       Integer                                               SLAP........669400
C         A flag to indicate symmetric matrix storage.                   SLAP........669500
C         If ISYM=0, all non-zero entries of the matrix are              SLAP........669600
C         stored.  If ISYM=1, the matrix is symmetric and                SLAP........669700
C         only the upper or lower triangular part is stored.             SLAP........669800
C                                                                        SLAP........669900
C***SEE ALSO  DGMRES                                                     SLAP........670000
C***ROUTINES CALLED  DAXPY, DCOPY, DHELS                                 SLAP........670100
C***REVISION HISTORY  (YYMMDD)                                           SLAP........670200
C   890404  DATE WRITTEN                                                 SLAP........670300
C   890404  Previous REVISION DATE                                       SLAP........670400
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........670500
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........670600
C           standard.  (FNF)                                             SLAP........670700
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........670800
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........670900
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)             SLAP........671000
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........671100
C   920511  Added complete declaration section.  (WRB)                   SLAP........671200
C***END PROLOGUE  DXLCAL                                                 SLAP........671300
C         The following is for optimized compilation on LLNL/LTSS Crays. SLAP........671400
CLLL. OPTIMIZE                                                           SLAP........671500
C     .. Scalar Arguments ..                                             SLAP........671600
      DOUBLE PRECISION R0NRM                                             SLAP........671700
      INTEGER ISYM, JPRE, JSCAL, LGMR, MAXLP1, N, NELT, NMSL             SLAP........671800
C     .. Array Arguments ..                                              SLAP........671900
      DOUBLE PRECISION A(NELT), HES(MAXLP1,*), Q(*), RPAR(*), SZ(*),     SLAP........672000
     +                 V(N,*), WK(N), X(N), XL(N), ZL(N)                 SLAP........672100
      INTEGER IA(NELT), IPAR(*), JA(NELT)                                SLAP........672200
C     .. Subroutine Arguments ..                                         SLAP........672300
      EXTERNAL MSOLVE                                                    SLAP........672400
C     .. Local Scalars ..                                                SLAP........672500
      INTEGER I, K, LL, LLP1                                             SLAP........672600
C     .. External Subroutines ..                                         SLAP........672700
      EXTERNAL DAXPY, DCOPY, DHELS                                       SLAP........672800
C***FIRST EXECUTABLE STATEMENT  DXLCAL                                   SLAP........672900
      LL = LGMR                                                          SLAP........673000
      LLP1 = LL + 1                                                      SLAP........673100
      DO 10 K = 1,LLP1                                                   SLAP........673200
         WK(K) = 0                                                       SLAP........673300
 10   CONTINUE                                                           SLAP........673400
      WK(1) = R0NRM                                                      SLAP........673500
      CALL DHELS(HES, MAXLP1, LL, Q, WK)                                 SLAP........673600
      DO 20 K = 1,N                                                      SLAP........673700
         ZL(K) = 0                                                       SLAP........673800
 20   CONTINUE                                                           SLAP........673900
      DO 30 I = 1,LL                                                     SLAP........674000
         CALL DAXPY(N, WK(I), V(1,I), 1, ZL, 1)                          SLAP........674100
 30   CONTINUE                                                           SLAP........674200
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN                        SLAP........674300
         DO 40 K = 1,N                                                   SLAP........674400
            ZL(K) = ZL(K)/SZ(K)                                          SLAP........674500
 40      CONTINUE                                                        SLAP........674600
      ENDIF                                                              SLAP........674700
      IF (JPRE .GT. 0) THEN                                              SLAP........674800
         CALL DCOPY(N, ZL, 1, WK, 1)                                     SLAP........674900
         CALL MSOLVE(N, WK, ZL, NELT, IA, JA, A, ISYM, RPAR, IPAR)       SLAP........675000
         NMSL = NMSL + 1                                                 SLAP........675100
      ENDIF                                                              SLAP........675200
C         calculate XL from X and ZL.                                    SLAP........675300
      DO 50 K = 1,N                                                      SLAP........675400
         XL(K) = X(K) + ZL(K)                                            SLAP........675500
 50   CONTINUE                                                           SLAP........675600
      RETURN                                                             SLAP........675700
C------------- LAST LINE OF DXLCAL FOLLOWS ----------------------------  SLAP........675800
      END                                                                SLAP........675900
*DECK FDUMP                                                              SLAP........676000
      SUBROUTINE FDUMP                                                   SLAP........676100
C***BEGIN PROLOGUE  FDUMP                                                SLAP........676200
C***PURPOSE  Symbolic dump (should be locally written).                  SLAP........676300
C***LIBRARY   SLATEC (XERROR)                                            SLAP........676400
C***CATEGORY  R3                                                         SLAP........676500
C***TYPE      ALL (FDUMP-A)                                              SLAP........676600
C***KEYWORDS  ERROR, XERMSG                                              SLAP........676700
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........676800
C***DESCRIPTION                                                          SLAP........676900
C                                                                        SLAP........677000
C        ***Note*** Machine Dependent Routine                            SLAP........677100
C        FDUMP is intended to be replaced by a locally written           SLAP........677200
C        version which produces a symbolic dump.  Failing this,          SLAP........677300
C        it should be replaced by a version which prints the             SLAP........677400
C        subprogram nesting list.  Note that this dump must be           SLAP........677500
C        printed on each of up to five files, as indicated by the        SLAP........677600
C        XGETUA routine.  See XSETUA and XGETUA for details.             SLAP........677700
C                                                                        SLAP........677800
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee SLAP........677900
C                                                                        SLAP........678000
C***REFERENCES  (NONE)                                                   SLAP........678100
C***ROUTINES CALLED  (NONE)                                              SLAP........678200
C***REVISION HISTORY  (YYMMDD)                                           SLAP........678300
C   790801  DATE WRITTEN                                                 SLAP........678400
C   861211  REVISION DATE from Version 3.2                               SLAP........678500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........678600
C***END PROLOGUE  FDUMP                                                  SLAP........678700
C***FIRST EXECUTABLE STATEMENT  FDUMP                                    SLAP........678800
      RETURN                                                             SLAP........678900
      END                                                                SLAP........679000
*DECK I1MACH                                                             SLAP........679100
      INTEGER FUNCTION I1MACH (I)                                        SLAP........679200
C***BEGIN PROLOGUE  I1MACH                                               SLAP........679300
C***PURPOSE  Return integer machine dependent constants.                 SLAP........679400
C***LIBRARY   SLATEC                                                     SLAP........679500
C***CATEGORY  R1                                                         SLAP........679600
C***TYPE      INTEGER (I1MACH-I)                                         SLAP........679700
C***KEYWORDS  MACHINE CONSTANTS                                          SLAP........679800
C***AUTHOR  Fox, P. A., (Bell Labs)                                      SLAP........679900
C           Hall, A. D., (Bell Labs)                                     SLAP........680000
C           Schryer, N. L., (Bell Labs)                                  SLAP........680100
C***DESCRIPTION                                                          SLAP........680200
C                                                                        SLAP........680300
C   I1MACH can be used to obtain machine-dependent parameters for the    SLAP........680400
C   local machine environment.  It is a function subprogram with one     SLAP........680500
C   (input) argument and can be referenced as follows:                   SLAP........680600
C                                                                        SLAP........680700
C        K = I1MACH(I)                                                   SLAP........680800
C                                                                        SLAP........680900
C   where I=1,...,16.  The (output) value of K above is determined by    SLAP........681000
C   the (input) value of I.  The results for various values of I are     SLAP........681100
C   discussed below.                                                     SLAP........681200
C                                                                        SLAP........681300
C   I/O unit numbers:                                                    SLAP........681400
C     I1MACH( 1) = the standard input unit.                              SLAP........681500
C     I1MACH( 2) = the standard output unit.                             SLAP........681600
C     I1MACH( 3) = the standard punch unit.                              SLAP........681700
C     I1MACH( 4) = the standard error message unit.                      SLAP........681800
C                                                                        SLAP........681900
C   Words:                                                               SLAP........682000
C     I1MACH( 5) = the number of bits per integer storage unit.          SLAP........682100
C     I1MACH( 6) = the number of characters per integer storage unit.    SLAP........682200
C                                                                        SLAP........682300
C   Integers:                                                            SLAP........682400
C     assume integers are represented in the S-digit, base-A form        SLAP........682500
C                                                                        SLAP........682600
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )          SLAP........682700
C                                                                        SLAP........682800
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.               SLAP........682900
C     I1MACH( 7) = A, the base.                                          SLAP........683000
C     I1MACH( 8) = S, the number of base-A digits.                       SLAP........683100
C     I1MACH( 9) = A**S - 1, the largest magnitude.                      SLAP........683200
C                                                                        SLAP........683300
C   Floating-Point Numbers:                                              SLAP........683400
C     Assume floating-point numbers are represented in the T-digit,      SLAP........683500
C     base-B form                                                        SLAP........683600
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )            SLAP........683700
C                                                                        SLAP........683800
C                where 0 .LE. X(I) .LT. B for I=1,...,T,                 SLAP........683900
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.                 SLAP........684000
C     I1MACH(10) = B, the base.                                          SLAP........684100
C                                                                        SLAP........684200
C   Single-Precision:                                                    SLAP........684300
C     I1MACH(11) = T, the number of base-B digits.                       SLAP........684400
C     I1MACH(12) = EMIN, the smallest exponent E.                        SLAP........684500
C     I1MACH(13) = EMAX, the largest exponent E.                         SLAP........684600
C                                                                        SLAP........684700
C   Double-Precision:                                                    SLAP........684800
C     I1MACH(14) = T, the number of base-B digits.                       SLAP........684900
C     I1MACH(15) = EMIN, the smallest exponent E.                        SLAP........685000
C     I1MACH(16) = EMAX, the largest exponent E.                         SLAP........685100
C                                                                        SLAP........685200
C   To alter this function for a particular environment, the desired     SLAP........685300
C   set of DATA statements should be activated by removing the C from    SLAP........685400
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be       SLAP........685500
C   checked for consistency with the local operating system.             SLAP........685600
C                                                                        SLAP........685700
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for   SLAP........685800
C                 a portable library, ACM Transactions on Mathematical   SLAP........685900
C                 Software 4, 2 (June 1978), pp. 177-188.                SLAP........686000
C***ROUTINES CALLED  (NONE)                                              SLAP........686100
C***REVISION HISTORY  (YYMMDD)                                           SLAP........686200
C   750101  DATE WRITTEN                                                 SLAP........686300
C   891012  Added VAX G-floating constants.  (WRB)                       SLAP........686400
C   891012  REVISION DATE from Version 3.2                               SLAP........686500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........686600
C   900618  Added DEC RISC constants.  (WRB)                             SLAP........686700
C   900723  Added IBM RS 6000 constants.  (WRB)                          SLAP........686800
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.    SLAP........686900
C           (RWC)                                                        SLAP........687000
C   910710  Added HP 730 constants.  (SMR)                               SLAP........687100
C   911114  Added Convex IEEE constants.  (WRB)                          SLAP........687200
C   920121  Added SUN -r8 compiler option constants.  (WRB)              SLAP........687300
C   920229  Added Touchstone Delta i860 constants.  (WRB)                SLAP........687400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........687500
C   920625  Added Convex -p8 and -pd8 compiler option constants.         SLAP........687600
C           (BKS, WRB)                                                   SLAP........687700
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)            SLAP........687800
C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler         SLAP........687900
C           options.  (DWL, RWC and WRB).                                SLAP........688000
C***END PROLOGUE  I1MACH                                                 SLAP........688100
C                                                                        SLAP........688200
      INTEGER IMACH(16),OUTPUT                                           SLAP........688300
      SAVE IMACH                                                         SLAP........688400
      EQUIVALENCE (IMACH(4),OUTPUT)                                      SLAP........688500
C                                                                        SLAP........688600
C     MACHINE CONSTANTS FOR THE AMIGA                                    SLAP........688700
C     ABSOFT COMPILER                                                    SLAP........688800
C                                                                        SLAP........688900
C     DATA IMACH( 1) /          5 /                                      SLAP........689000
C     DATA IMACH( 2) /          6 /                                      SLAP........689100
C     DATA IMACH( 3) /          5 /                                      SLAP........689200
C     DATA IMACH( 4) /          6 /                                      SLAP........689300
C     DATA IMACH( 5) /         32 /                                      SLAP........689400
C     DATA IMACH( 6) /          4 /                                      SLAP........689500
C     DATA IMACH( 7) /          2 /                                      SLAP........689600
C     DATA IMACH( 8) /         31 /                                      SLAP........689700
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........689800
C     DATA IMACH(10) /          2 /                                      SLAP........689900
C     DATA IMACH(11) /         24 /                                      SLAP........690000
C     DATA IMACH(12) /       -126 /                                      SLAP........690100
C     DATA IMACH(13) /        127 /                                      SLAP........690200
C     DATA IMACH(14) /         53 /                                      SLAP........690300
C     DATA IMACH(15) /      -1022 /                                      SLAP........690400
C     DATA IMACH(16) /       1023 /                                      SLAP........690500
C                                                                        SLAP........690600
C     MACHINE CONSTANTS FOR THE APOLLO                                   SLAP........690700
C                                                                        SLAP........690800
C     DATA IMACH( 1) /          5 /                                      SLAP........690900
C     DATA IMACH( 2) /          6 /                                      SLAP........691000
C     DATA IMACH( 3) /          6 /                                      SLAP........691100
C     DATA IMACH( 4) /          6 /                                      SLAP........691200
C     DATA IMACH( 5) /         32 /                                      SLAP........691300
C     DATA IMACH( 6) /          4 /                                      SLAP........691400
C     DATA IMACH( 7) /          2 /                                      SLAP........691500
C     DATA IMACH( 8) /         31 /                                      SLAP........691600
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........691700
C     DATA IMACH(10) /          2 /                                      SLAP........691800
C     DATA IMACH(11) /         24 /                                      SLAP........691900
C     DATA IMACH(12) /       -125 /                                      SLAP........692000
C     DATA IMACH(13) /        129 /                                      SLAP........692100
C     DATA IMACH(14) /         53 /                                      SLAP........692200
C     DATA IMACH(15) /      -1021 /                                      SLAP........692300
C     DATA IMACH(16) /       1025 /                                      SLAP........692400
C                                                                        SLAP........692500
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM                    SLAP........692600
C                                                                        SLAP........692700
C     DATA IMACH( 1) /          7 /                                      SLAP........692800
C     DATA IMACH( 2) /          2 /                                      SLAP........692900
C     DATA IMACH( 3) /          2 /                                      SLAP........693000
C     DATA IMACH( 4) /          2 /                                      SLAP........693100
C     DATA IMACH( 5) /         36 /                                      SLAP........693200
C     DATA IMACH( 6) /          4 /                                      SLAP........693300
C     DATA IMACH( 7) /          2 /                                      SLAP........693400
C     DATA IMACH( 8) /         33 /                                      SLAP........693500
C     DATA IMACH( 9) / Z1FFFFFFFF /                                      SLAP........693600
C     DATA IMACH(10) /          2 /                                      SLAP........693700
C     DATA IMACH(11) /         24 /                                      SLAP........693800
C     DATA IMACH(12) /       -256 /                                      SLAP........693900
C     DATA IMACH(13) /        255 /                                      SLAP........694000
C     DATA IMACH(14) /         60 /                                      SLAP........694100
C     DATA IMACH(15) /       -256 /                                      SLAP........694200
C     DATA IMACH(16) /        255 /                                      SLAP........694300
C                                                                        SLAP........694400
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM                    SLAP........694500
C                                                                        SLAP........694600
C     DATA IMACH( 1) /          5 /                                      SLAP........694700
C     DATA IMACH( 2) /          6 /                                      SLAP........694800
C     DATA IMACH( 3) /          7 /                                      SLAP........694900
C     DATA IMACH( 4) /          6 /                                      SLAP........695000
C     DATA IMACH( 5) /         48 /                                      SLAP........695100
C     DATA IMACH( 6) /          6 /                                      SLAP........695200
C     DATA IMACH( 7) /          2 /                                      SLAP........695300
C     DATA IMACH( 8) /         39 /                                      SLAP........695400
C     DATA IMACH( 9) / O0007777777777777 /                               SLAP........695500
C     DATA IMACH(10) /          8 /                                      SLAP........695600
C     DATA IMACH(11) /         13 /                                      SLAP........695700
C     DATA IMACH(12) /        -50 /                                      SLAP........695800
C     DATA IMACH(13) /         76 /                                      SLAP........695900
C     DATA IMACH(14) /         26 /                                      SLAP........696000
C     DATA IMACH(15) /        -50 /                                      SLAP........696100
C     DATA IMACH(16) /         76 /                                      SLAP........696200
C                                                                        SLAP........696300
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS              SLAP........696400
C                                                                        SLAP........696500
C     DATA IMACH( 1) /          5 /                                      SLAP........696600
C     DATA IMACH( 2) /          6 /                                      SLAP........696700
C     DATA IMACH( 3) /          7 /                                      SLAP........696800
C     DATA IMACH( 4) /          6 /                                      SLAP........696900
C     DATA IMACH( 5) /         48 /                                      SLAP........697000
C     DATA IMACH( 6) /          6 /                                      SLAP........697100
C     DATA IMACH( 7) /          2 /                                      SLAP........697200
C     DATA IMACH( 8) /         39 /                                      SLAP........697300
C     DATA IMACH( 9) / O0007777777777777 /                               SLAP........697400
C     DATA IMACH(10) /          8 /                                      SLAP........697500
C     DATA IMACH(11) /         13 /                                      SLAP........697600
C     DATA IMACH(12) /        -50 /                                      SLAP........697700
C     DATA IMACH(13) /         76 /                                      SLAP........697800
C     DATA IMACH(14) /         26 /                                      SLAP........697900
C     DATA IMACH(15) /     -32754 /                                      SLAP........698000
C     DATA IMACH(16) /      32780 /                                      SLAP........698100
C                                                                        SLAP........698200
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE          SLAP........698300
C                                                                        SLAP........698400
C     DATA IMACH( 1) /          5 /                                      SLAP........698500
C     DATA IMACH( 2) /          6 /                                      SLAP........698600
C     DATA IMACH( 3) /          7 /                                      SLAP........698700
C     DATA IMACH( 4) /          6 /                                      SLAP........698800
C     DATA IMACH( 5) /         64 /                                      SLAP........698900
C     DATA IMACH( 6) /          8 /                                      SLAP........699000
C     DATA IMACH( 7) /          2 /                                      SLAP........699100
C     DATA IMACH( 8) /         63 /                                      SLAP........699200
C     DATA IMACH( 9) / 9223372036854775807 /                             SLAP........699300
C     DATA IMACH(10) /          2 /                                      SLAP........699400
C     DATA IMACH(11) /         47 /                                      SLAP........699500
C     DATA IMACH(12) /      -4095 /                                      SLAP........699600
C     DATA IMACH(13) /       4094 /                                      SLAP........699700
C     DATA IMACH(14) /         94 /                                      SLAP........699800
C     DATA IMACH(15) /      -4095 /                                      SLAP........699900
C     DATA IMACH(16) /       4094 /                                      SLAP........700000
C                                                                        SLAP........700100
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES                     SLAP........700200
C                                                                        SLAP........700300
C     DATA IMACH( 1) /          5 /                                      SLAP........700400
C     DATA IMACH( 2) /          6 /                                      SLAP........700500
C     DATA IMACH( 3) /          7 /                                      SLAP........700600
C     DATA IMACH( 4) /    6LOUTPUT/                                      SLAP........700700
C     DATA IMACH( 5) /         60 /                                      SLAP........700800
C     DATA IMACH( 6) /         10 /                                      SLAP........700900
C     DATA IMACH( 7) /          2 /                                      SLAP........701000
C     DATA IMACH( 8) /         48 /                                      SLAP........701100
C     DATA IMACH( 9) / 00007777777777777777B /                           SLAP........701200
C     DATA IMACH(10) /          2 /                                      SLAP........701300
C     DATA IMACH(11) /         47 /                                      SLAP........701400
C     DATA IMACH(12) /       -929 /                                      SLAP........701500
C     DATA IMACH(13) /       1070 /                                      SLAP........701600
C     DATA IMACH(14) /         94 /                                      SLAP........701700
C     DATA IMACH(15) /       -929 /                                      SLAP........701800
C     DATA IMACH(16) /       1069 /                                      SLAP........701900
C                                                                        SLAP........702000
C     MACHINE CONSTANTS FOR THE CELERITY C1260                           SLAP........702100
C                                                                        SLAP........702200
C     DATA IMACH( 1) /          5 /                                      SLAP........702300
C     DATA IMACH( 2) /          6 /                                      SLAP........702400
C     DATA IMACH( 3) /          6 /                                      SLAP........702500
C     DATA IMACH( 4) /          0 /                                      SLAP........702600
C     DATA IMACH( 5) /         32 /                                      SLAP........702700
C     DATA IMACH( 6) /          4 /                                      SLAP........702800
C     DATA IMACH( 7) /          2 /                                      SLAP........702900
C     DATA IMACH( 8) /         31 /                                      SLAP........703000
C     DATA IMACH( 9) / Z'7FFFFFFF' /                                     SLAP........703100
C     DATA IMACH(10) /          2 /                                      SLAP........703200
C     DATA IMACH(11) /         24 /                                      SLAP........703300
C     DATA IMACH(12) /       -126 /                                      SLAP........703400
C     DATA IMACH(13) /        127 /                                      SLAP........703500
C     DATA IMACH(14) /         53 /                                      SLAP........703600
C     DATA IMACH(15) /      -1022 /                                      SLAP........703700
C     DATA IMACH(16) /       1023 /                                      SLAP........703800
C                                                                        SLAP........703900
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP........704000
C     USING THE -fn COMPILER OPTION                                      SLAP........704100
C                                                                        SLAP........704200
C     DATA IMACH( 1) /          5 /                                      SLAP........704300
C     DATA IMACH( 2) /          6 /                                      SLAP........704400
C     DATA IMACH( 3) /          7 /                                      SLAP........704500
C     DATA IMACH( 4) /          6 /                                      SLAP........704600
C     DATA IMACH( 5) /         32 /                                      SLAP........704700
C     DATA IMACH( 6) /          4 /                                      SLAP........704800
C     DATA IMACH( 7) /          2 /                                      SLAP........704900
C     DATA IMACH( 8) /         31 /                                      SLAP........705000
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........705100
C     DATA IMACH(10) /          2 /                                      SLAP........705200
C     DATA IMACH(11) /         24 /                                      SLAP........705300
C     DATA IMACH(12) /       -127 /                                      SLAP........705400
C     DATA IMACH(13) /        127 /                                      SLAP........705500
C     DATA IMACH(14) /         53 /                                      SLAP........705600
C     DATA IMACH(15) /      -1023 /                                      SLAP........705700
C     DATA IMACH(16) /       1023 /                                      SLAP........705800
C                                                                        SLAP........705900
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP........706000
C     USING THE -fi COMPILER OPTION                                      SLAP........706100
C                                                                        SLAP........706200
C     DATA IMACH( 1) /          5 /                                      SLAP........706300
C     DATA IMACH( 2) /          6 /                                      SLAP........706400
C     DATA IMACH( 3) /          7 /                                      SLAP........706500
C     DATA IMACH( 4) /          6 /                                      SLAP........706600
C     DATA IMACH( 5) /         32 /                                      SLAP........706700
C     DATA IMACH( 6) /          4 /                                      SLAP........706800
C     DATA IMACH( 7) /          2 /                                      SLAP........706900
C     DATA IMACH( 8) /         31 /                                      SLAP........707000
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........707100
C     DATA IMACH(10) /          2 /                                      SLAP........707200
C     DATA IMACH(11) /         24 /                                      SLAP........707300
C     DATA IMACH(12) /       -125 /                                      SLAP........707400
C     DATA IMACH(13) /        128 /                                      SLAP........707500
C     DATA IMACH(14) /         53 /                                      SLAP........707600
C     DATA IMACH(15) /      -1021 /                                      SLAP........707700
C     DATA IMACH(16) /       1024 /                                      SLAP........707800
C                                                                        SLAP........707900
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP........708000
C     USING THE -p8 COMPILER OPTION                                      SLAP........708100
C                                                                        SLAP........708200
C     DATA IMACH( 1) /          5 /                                      SLAP........708300
C     DATA IMACH( 2) /          6 /                                      SLAP........708400
C     DATA IMACH( 3) /          7 /                                      SLAP........708500
C     DATA IMACH( 4) /          6 /                                      SLAP........708600
C     DATA IMACH( 5) /         64 /                                      SLAP........708700
C     DATA IMACH( 6) /          4 /                                      SLAP........708800
C     DATA IMACH( 7) /          2 /                                      SLAP........708900
C     DATA IMACH( 8) /         63 /                                      SLAP........709000
C     DATA IMACH( 9) / 9223372036854775807 /                             SLAP........709100
C     DATA IMACH(10) /          2 /                                      SLAP........709200
C     DATA IMACH(11) /         53 /                                      SLAP........709300
C     DATA IMACH(12) /      -1023 /                                      SLAP........709400
C     DATA IMACH(13) /       1023 /                                      SLAP........709500
C     DATA IMACH(14) /        113 /                                      SLAP........709600
C     DATA IMACH(15) /     -16383 /                                      SLAP........709700
C     DATA IMACH(16) /      16383 /                                      SLAP........709800
C                                                                        SLAP........709900
C     MACHINE CONSTANTS FOR THE CONVEX                                   SLAP........710000
C     USING THE -pd8 COMPILER OPTION                                     SLAP........710100
C                                                                        SLAP........710200
C     DATA IMACH( 1) /          5 /                                      SLAP........710300
C     DATA IMACH( 2) /          6 /                                      SLAP........710400
C     DATA IMACH( 3) /          7 /                                      SLAP........710500
C     DATA IMACH( 4) /          6 /                                      SLAP........710600
C     DATA IMACH( 5) /         64 /                                      SLAP........710700
C     DATA IMACH( 6) /          4 /                                      SLAP........710800
C     DATA IMACH( 7) /          2 /                                      SLAP........710900
C     DATA IMACH( 8) /         63 /                                      SLAP........711000
C     DATA IMACH( 9) / 9223372036854775807 /                             SLAP........711100
C     DATA IMACH(10) /          2 /                                      SLAP........711200
C     DATA IMACH(11) /         53 /                                      SLAP........711300
C     DATA IMACH(12) /      -1023 /                                      SLAP........711400
C     DATA IMACH(13) /       1023 /                                      SLAP........711500
C     DATA IMACH(14) /         53 /                                      SLAP........711600
C     DATA IMACH(15) /      -1023 /                                      SLAP........711700
C     DATA IMACH(16) /       1023 /                                      SLAP........711800
C                                                                        SLAP........711900
C     MACHINE CONSTANTS FOR THE CRAY                                     SLAP........712000
C     USING THE 46 BIT INTEGER COMPILER OPTION                           SLAP........712100
C                                                                        SLAP........712200
C     DATA IMACH( 1) /        100 /                                      SLAP........712300
C     DATA IMACH( 2) /        101 /                                      SLAP........712400
C     DATA IMACH( 3) /        102 /                                      SLAP........712500
C     DATA IMACH( 4) /        101 /                                      SLAP........712600
C     DATA IMACH( 5) /         64 /                                      SLAP........712700
C     DATA IMACH( 6) /          8 /                                      SLAP........712800
C     DATA IMACH( 7) /          2 /                                      SLAP........712900
C     DATA IMACH( 8) /         46 /                                      SLAP........713000
C     DATA IMACH( 9) / 1777777777777777B /                               SLAP........713100
C     DATA IMACH(10) /          2 /                                      SLAP........713200
C     DATA IMACH(11) /         47 /                                      SLAP........713300
C     DATA IMACH(12) /      -8189 /                                      SLAP........713400
C     DATA IMACH(13) /       8190 /                                      SLAP........713500
C     DATA IMACH(14) /         94 /                                      SLAP........713600
C     DATA IMACH(15) /      -8099 /                                      SLAP........713700
C     DATA IMACH(16) /       8190 /                                      SLAP........713800
C                                                                        SLAP........713900
C     MACHINE CONSTANTS FOR THE CRAY                                     SLAP........714000
C     USING THE 64 BIT INTEGER COMPILER OPTION                           SLAP........714100
C                                                                        SLAP........714200
C     DATA IMACH( 1) /        100 /                                      SLAP........714300
C     DATA IMACH( 2) /        101 /                                      SLAP........714400
C     DATA IMACH( 3) /        102 /                                      SLAP........714500
C     DATA IMACH( 4) /        101 /                                      SLAP........714600
C     DATA IMACH( 5) /         64 /                                      SLAP........714700
C     DATA IMACH( 6) /          8 /                                      SLAP........714800
C     DATA IMACH( 7) /          2 /                                      SLAP........714900
C     DATA IMACH( 8) /         63 /                                      SLAP........715000
C     DATA IMACH( 9) / 777777777777777777777B /                          SLAP........715100
C     DATA IMACH(10) /          2 /                                      SLAP........715200
C     DATA IMACH(11) /         47 /                                      SLAP........715300
C     DATA IMACH(12) /      -8189 /                                      SLAP........715400
C     DATA IMACH(13) /       8190 /                                      SLAP........715500
C     DATA IMACH(14) /         94 /                                      SLAP........715600
C     DATA IMACH(15) /      -8099 /                                      SLAP........715700
C     DATA IMACH(16) /       8190 /                                      SLAP........715800
C                                                                        SLAP........715900
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200               SLAP........716000
C                                                                        SLAP........716100
C     DATA IMACH( 1) /         11 /                                      SLAP........716200
C     DATA IMACH( 2) /         12 /                                      SLAP........716300
C     DATA IMACH( 3) /          8 /                                      SLAP........716400
C     DATA IMACH( 4) /         10 /                                      SLAP........716500
C     DATA IMACH( 5) /         16 /                                      SLAP........716600
C     DATA IMACH( 6) /          2 /                                      SLAP........716700
C     DATA IMACH( 7) /          2 /                                      SLAP........716800
C     DATA IMACH( 8) /         15 /                                      SLAP........716900
C     DATA IMACH( 9) /      32767 /                                      SLAP........717000
C     DATA IMACH(10) /         16 /                                      SLAP........717100
C     DATA IMACH(11) /          6 /                                      SLAP........717200
C     DATA IMACH(12) /        -64 /                                      SLAP........717300
C     DATA IMACH(13) /         63 /                                      SLAP........717400
C     DATA IMACH(14) /         14 /                                      SLAP........717500
C     DATA IMACH(15) /        -64 /                                      SLAP........717600
C     DATA IMACH(16) /         63 /                                      SLAP........717700
C                                                                        SLAP........717800
C     MACHINE CONSTANTS FOR THE DEC ALPHA                                SLAP........717900
C     USING G_FLOAT                                                      SLAP........718000
C                                                                        SLAP........718100
C     DATA IMACH( 1) /          5 /                                      SLAP........718200
C     DATA IMACH( 2) /          6 /                                      SLAP........718300
C     DATA IMACH( 3) /          5 /                                      SLAP........718400
C     DATA IMACH( 4) /          6 /                                      SLAP........718500
C     DATA IMACH( 5) /         32 /                                      SLAP........718600
C     DATA IMACH( 6) /          4 /                                      SLAP........718700
C     DATA IMACH( 7) /          2 /                                      SLAP........718800
C     DATA IMACH( 8) /         31 /                                      SLAP........718900
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........719000
C     DATA IMACH(10) /          2 /                                      SLAP........719100
C     DATA IMACH(11) /         24 /                                      SLAP........719200
C     DATA IMACH(12) /       -127 /                                      SLAP........719300
C     DATA IMACH(13) /        127 /                                      SLAP........719400
C     DATA IMACH(14) /         53 /                                      SLAP........719500
C     DATA IMACH(15) /      -1023 /                                      SLAP........719600
C     DATA IMACH(16) /       1023 /                                      SLAP........719700
C                                                                        SLAP........719800
C     MACHINE CONSTANTS FOR THE DEC ALPHA                                SLAP........719900
C     USING IEEE_FLOAT                                                   SLAP........720000
C                                                                        SLAP........720100
C     DATA IMACH( 1) /          5 /                                      SLAP........720200
C     DATA IMACH( 2) /          6 /                                      SLAP........720300
C     DATA IMACH( 3) /          6 /                                      SLAP........720400
C     DATA IMACH( 4) /          6 /                                      SLAP........720500
C     DATA IMACH( 5) /         32 /                                      SLAP........720600
C     DATA IMACH( 6) /          4 /                                      SLAP........720700
C     DATA IMACH( 7) /          2 /                                      SLAP........720800
C     DATA IMACH( 8) /         31 /                                      SLAP........720900
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........721000
C     DATA IMACH(10) /          2 /                                      SLAP........721100
C     DATA IMACH(11) /         24 /                                      SLAP........721200
C     DATA IMACH(12) /       -125 /                                      SLAP........721300
C     DATA IMACH(13) /        128 /                                      SLAP........721400
C     DATA IMACH(14) /         53 /                                      SLAP........721500
C     DATA IMACH(15) /      -1021 /                                      SLAP........721600
C     DATA IMACH(16) /       1024 /                                      SLAP........721700
C                                                                        SLAP........721800
C     MACHINE CONSTANTS FOR THE DEC RISC                                 SLAP........721900
C                                                                        SLAP........722000
C     DATA IMACH( 1) /          5 /                                      SLAP........722100
C     DATA IMACH( 2) /          6 /                                      SLAP........722200
C     DATA IMACH( 3) /          6 /                                      SLAP........722300
C     DATA IMACH( 4) /          6 /                                      SLAP........722400
C     DATA IMACH( 5) /         32 /                                      SLAP........722500
C     DATA IMACH( 6) /          4 /                                      SLAP........722600
C     DATA IMACH( 7) /          2 /                                      SLAP........722700
C     DATA IMACH( 8) /         31 /                                      SLAP........722800
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........722900
C     DATA IMACH(10) /          2 /                                      SLAP........723000
C     DATA IMACH(11) /         24 /                                      SLAP........723100
C     DATA IMACH(12) /       -125 /                                      SLAP........723200
C     DATA IMACH(13) /        128 /                                      SLAP........723300
C     DATA IMACH(14) /         53 /                                      SLAP........723400
C     DATA IMACH(15) /      -1021 /                                      SLAP........723500
C     DATA IMACH(16) /       1024 /                                      SLAP........723600
C                                                                        SLAP........723700
C     MACHINE CONSTANTS FOR THE DEC VAX                                  SLAP........723800
C     USING D_FLOATING                                                   SLAP........723900
C                                                                        SLAP........724000
C     DATA IMACH( 1) /          5 /                                      SLAP........724100
C     DATA IMACH( 2) /          6 /                                      SLAP........724200
C     DATA IMACH( 3) /          5 /                                      SLAP........724300
C     DATA IMACH( 4) /          6 /                                      SLAP........724400
C     DATA IMACH( 5) /         32 /                                      SLAP........724500
C     DATA IMACH( 6) /          4 /                                      SLAP........724600
C     DATA IMACH( 7) /          2 /                                      SLAP........724700
C     DATA IMACH( 8) /         31 /                                      SLAP........724800
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........724900
C     DATA IMACH(10) /          2 /                                      SLAP........725000
C     DATA IMACH(11) /         24 /                                      SLAP........725100
C     DATA IMACH(12) /       -127 /                                      SLAP........725200
C     DATA IMACH(13) /        127 /                                      SLAP........725300
C     DATA IMACH(14) /         56 /                                      SLAP........725400
C     DATA IMACH(15) /       -127 /                                      SLAP........725500
C     DATA IMACH(16) /        127 /                                      SLAP........725600
C                                                                        SLAP........725700
C     MACHINE CONSTANTS FOR THE DEC VAX                                  SLAP........725800
C     USING G_FLOATING                                                   SLAP........725900
C                                                                        SLAP........726000
C     DATA IMACH( 1) /          5 /                                      SLAP........726100
C     DATA IMACH( 2) /          6 /                                      SLAP........726200
C     DATA IMACH( 3) /          5 /                                      SLAP........726300
C     DATA IMACH( 4) /          6 /                                      SLAP........726400
C     DATA IMACH( 5) /         32 /                                      SLAP........726500
C     DATA IMACH( 6) /          4 /                                      SLAP........726600
C     DATA IMACH( 7) /          2 /                                      SLAP........726700
C     DATA IMACH( 8) /         31 /                                      SLAP........726800
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........726900
C     DATA IMACH(10) /          2 /                                      SLAP........727000
C     DATA IMACH(11) /         24 /                                      SLAP........727100
C     DATA IMACH(12) /       -127 /                                      SLAP........727200
C     DATA IMACH(13) /        127 /                                      SLAP........727300
C     DATA IMACH(14) /         53 /                                      SLAP........727400
C     DATA IMACH(15) /      -1023 /                                      SLAP........727500
C     DATA IMACH(16) /       1023 /                                      SLAP........727600
C                                                                        SLAP........727700
C     MACHINE CONSTANTS FOR THE ELXSI 6400                               SLAP........727800
C                                                                        SLAP........727900
C     DATA IMACH( 1) /          5 /                                      SLAP........728000
C     DATA IMACH( 2) /          6 /                                      SLAP........728100
C     DATA IMACH( 3) /          6 /                                      SLAP........728200
C     DATA IMACH( 4) /          6 /                                      SLAP........728300
C     DATA IMACH( 5) /         32 /                                      SLAP........728400
C     DATA IMACH( 6) /          4 /                                      SLAP........728500
C     DATA IMACH( 7) /          2 /                                      SLAP........728600
C     DATA IMACH( 8) /         32 /                                      SLAP........728700
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........728800
C     DATA IMACH(10) /          2 /                                      SLAP........728900
C     DATA IMACH(11) /         24 /                                      SLAP........729000
C     DATA IMACH(12) /       -126 /                                      SLAP........729100
C     DATA IMACH(13) /        127 /                                      SLAP........729200
C     DATA IMACH(14) /         53 /                                      SLAP........729300
C     DATA IMACH(15) /      -1022 /                                      SLAP........729400
C     DATA IMACH(16) /       1023 /                                      SLAP........729500
C                                                                        SLAP........729600
C     MACHINE CONSTANTS FOR THE HARRIS 220                               SLAP........729700
C                                                                        SLAP........729800
C     DATA IMACH( 1) /          5 /                                      SLAP........729900
C     DATA IMACH( 2) /          6 /                                      SLAP........730000
C     DATA IMACH( 3) /          0 /                                      SLAP........730100
C     DATA IMACH( 4) /          6 /                                      SLAP........730200
C     DATA IMACH( 5) /         24 /                                      SLAP........730300
C     DATA IMACH( 6) /          3 /                                      SLAP........730400
C     DATA IMACH( 7) /          2 /                                      SLAP........730500
C     DATA IMACH( 8) /         23 /                                      SLAP........730600
C     DATA IMACH( 9) /    8388607 /                                      SLAP........730700
C     DATA IMACH(10) /          2 /                                      SLAP........730800
C     DATA IMACH(11) /         23 /                                      SLAP........730900
C     DATA IMACH(12) /       -127 /                                      SLAP........731000
C     DATA IMACH(13) /        127 /                                      SLAP........731100
C     DATA IMACH(14) /         38 /                                      SLAP........731200
C     DATA IMACH(15) /       -127 /                                      SLAP........731300
C     DATA IMACH(16) /        127 /                                      SLAP........731400
C                                                                        SLAP........731500
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES                SLAP........731600
C                                                                        SLAP........731700
C     DATA IMACH( 1) /          5 /                                      SLAP........731800
C     DATA IMACH( 2) /          6 /                                      SLAP........731900
C     DATA IMACH( 3) /         43 /                                      SLAP........732000
C     DATA IMACH( 4) /          6 /                                      SLAP........732100
C     DATA IMACH( 5) /         36 /                                      SLAP........732200
C     DATA IMACH( 6) /          6 /                                      SLAP........732300
C     DATA IMACH( 7) /          2 /                                      SLAP........732400
C     DATA IMACH( 8) /         35 /                                      SLAP........732500
C     DATA IMACH( 9) / O377777777777 /                                   SLAP........732600
C     DATA IMACH(10) /          2 /                                      SLAP........732700
C     DATA IMACH(11) /         27 /                                      SLAP........732800
C     DATA IMACH(12) /       -127 /                                      SLAP........732900
C     DATA IMACH(13) /        127 /                                      SLAP........733000
C     DATA IMACH(14) /         63 /                                      SLAP........733100
C     DATA IMACH(15) /       -127 /                                      SLAP........733200
C     DATA IMACH(16) /        127 /                                      SLAP........733300
C                                                                        SLAP........733400
C     MACHINE CONSTANTS FOR THE HP 730                                   SLAP........733500
C                                                                        SLAP........733600
C     DATA IMACH( 1) /          5 /                                      SLAP........733700
C     DATA IMACH( 2) /          6 /                                      SLAP........733800
C     DATA IMACH( 3) /          6 /                                      SLAP........733900
C     DATA IMACH( 4) /          6 /                                      SLAP........734000
C     DATA IMACH( 5) /         32 /                                      SLAP........734100
C     DATA IMACH( 6) /          4 /                                      SLAP........734200
C     DATA IMACH( 7) /          2 /                                      SLAP........734300
C     DATA IMACH( 8) /         31 /                                      SLAP........734400
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........734500
C     DATA IMACH(10) /          2 /                                      SLAP........734600
C     DATA IMACH(11) /         24 /                                      SLAP........734700
C     DATA IMACH(12) /       -125 /                                      SLAP........734800
C     DATA IMACH(13) /        128 /                                      SLAP........734900
C     DATA IMACH(14) /         53 /                                      SLAP........735000
C     DATA IMACH(15) /      -1021 /                                      SLAP........735100
C     DATA IMACH(16) /       1024 /                                      SLAP........735200
C                                                                        SLAP........735300
C     MACHINE CONSTANTS FOR THE HP 2100                                  SLAP........735400
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4                           SLAP........735500
C                                                                        SLAP........735600
C     DATA IMACH( 1) /          5 /                                      SLAP........735700
C     DATA IMACH( 2) /          6 /                                      SLAP........735800
C     DATA IMACH( 3) /          4 /                                      SLAP........735900
C     DATA IMACH( 4) /          1 /                                      SLAP........736000
C     DATA IMACH( 5) /         16 /                                      SLAP........736100
C     DATA IMACH( 6) /          2 /                                      SLAP........736200
C     DATA IMACH( 7) /          2 /                                      SLAP........736300
C     DATA IMACH( 8) /         15 /                                      SLAP........736400
C     DATA IMACH( 9) /      32767 /                                      SLAP........736500
C     DATA IMACH(10) /          2 /                                      SLAP........736600
C     DATA IMACH(11) /         23 /                                      SLAP........736700
C     DATA IMACH(12) /       -128 /                                      SLAP........736800
C     DATA IMACH(13) /        127 /                                      SLAP........736900
C     DATA IMACH(14) /         39 /                                      SLAP........737000
C     DATA IMACH(15) /       -128 /                                      SLAP........737100
C     DATA IMACH(16) /        127 /                                      SLAP........737200
C                                                                        SLAP........737300
C     MACHINE CONSTANTS FOR THE HP 2100                                  SLAP........737400
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4                           SLAP........737500
C                                                                        SLAP........737600
C     DATA IMACH( 1) /          5 /                                      SLAP........737700
C     DATA IMACH( 2) /          6 /                                      SLAP........737800
C     DATA IMACH( 3) /          4 /                                      SLAP........737900
C     DATA IMACH( 4) /          1 /                                      SLAP........738000
C     DATA IMACH( 5) /         16 /                                      SLAP........738100
C     DATA IMACH( 6) /          2 /                                      SLAP........738200
C     DATA IMACH( 7) /          2 /                                      SLAP........738300
C     DATA IMACH( 8) /         15 /                                      SLAP........738400
C     DATA IMACH( 9) /      32767 /                                      SLAP........738500
C     DATA IMACH(10) /          2 /                                      SLAP........738600
C     DATA IMACH(11) /         23 /                                      SLAP........738700
C     DATA IMACH(12) /       -128 /                                      SLAP........738800
C     DATA IMACH(13) /        127 /                                      SLAP........738900
C     DATA IMACH(14) /         55 /                                      SLAP........739000
C     DATA IMACH(15) /       -128 /                                      SLAP........739100
C     DATA IMACH(16) /        127 /                                      SLAP........739200
C                                                                        SLAP........739300
C     MACHINE CONSTANTS FOR THE HP 9000                                  SLAP........739400
C                                                                        SLAP........739500
C     DATA IMACH( 1) /          5 /                                      SLAP........739600
C     DATA IMACH( 2) /          6 /                                      SLAP........739700
C     DATA IMACH( 3) /          6 /                                      SLAP........739800
C     DATA IMACH( 4) /          7 /                                      SLAP........739900
C     DATA IMACH( 5) /         32 /                                      SLAP........740000
C     DATA IMACH( 6) /          4 /                                      SLAP........740100
C     DATA IMACH( 7) /          2 /                                      SLAP........740200
C     DATA IMACH( 8) /         32 /                                      SLAP........740300
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........740400
C     DATA IMACH(10) /          2 /                                      SLAP........740500
C     DATA IMACH(11) /         24 /                                      SLAP........740600
C     DATA IMACH(12) /       -126 /                                      SLAP........740700
C     DATA IMACH(13) /        127 /                                      SLAP........740800
C     DATA IMACH(14) /         53 /                                      SLAP........740900
C     DATA IMACH(15) /      -1015 /                                      SLAP........741000
C     DATA IMACH(16) /       1017 /                                      SLAP........741100
C                                                                        SLAP........741200
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                      SLAP........741300
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND                  SLAP........741400
C     THE PERKIN ELMER (INTERDATA) 7/32.                                 SLAP........741500
C                                                                        SLAP........741600
C     DATA IMACH( 1) /          5 /                                      SLAP........741700
C     DATA IMACH( 2) /          6 /                                      SLAP........741800
C     DATA IMACH( 3) /          7 /                                      SLAP........741900
C     DATA IMACH( 4) /          6 /                                      SLAP........742000
C     DATA IMACH( 5) /         32 /                                      SLAP........742100
C     DATA IMACH( 6) /          4 /                                      SLAP........742200
C     DATA IMACH( 7) /          2 /                                      SLAP........742300
C     DATA IMACH( 8) /         31 /                                      SLAP........742400
C     DATA IMACH( 9) /  Z7FFFFFFF /                                      SLAP........742500
C     DATA IMACH(10) /         16 /                                      SLAP........742600
C     DATA IMACH(11) /          6 /                                      SLAP........742700
C     DATA IMACH(12) /        -64 /                                      SLAP........742800
C     DATA IMACH(13) /         63 /                                      SLAP........742900
C     DATA IMACH(14) /         14 /                                      SLAP........743000
C     DATA IMACH(15) /        -64 /                                      SLAP........743100
C     DATA IMACH(16) /         63 /                                      SLAP........743200
C                                                                        SLAP........743300
C     MACHINE CONSTANTS FOR THE IBM PC                                   SLAP........743400
C                                                                        SLAP........743500
C     DATA IMACH( 1) /          5 /                                      SLAP........743600
C     DATA IMACH( 2) /          6 /                                      SLAP........743700
C     DATA IMACH( 3) /          0 /                                      SLAP........743800
C     DATA IMACH( 4) /          0 /                                      SLAP........743900
C     DATA IMACH( 5) /         32 /                                      SLAP........744000
C     DATA IMACH( 6) /          4 /                                      SLAP........744100
C     DATA IMACH( 7) /          2 /                                      SLAP........744200
C     DATA IMACH( 8) /         31 /                                      SLAP........744300
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........744400
C     DATA IMACH(10) /          2 /                                      SLAP........744500
C     DATA IMACH(11) /         24 /                                      SLAP........744600
C     DATA IMACH(12) /       -125 /                                      SLAP........744700
C     DATA IMACH(13) /        127 /                                      SLAP........744800
C     DATA IMACH(14) /         53 /                                      SLAP........744900
C     DATA IMACH(15) /      -1021 /                                      SLAP........745000
C     DATA IMACH(16) /       1023 /                                      SLAP........745100
C                                                                        SLAP........745200
C     MACHINE CONSTANTS FOR THE IBM RS 6000                              SLAP........745300
C                                                                        SLAP........745400
C     DATA IMACH( 1) /          5 /                                      SLAP........745500
C     DATA IMACH( 2) /          6 /                                      SLAP........745600
C     DATA IMACH( 3) /          6 /                                      SLAP........745700
C     DATA IMACH( 4) /          0 /                                      SLAP........745800
C     DATA IMACH( 5) /         32 /                                      SLAP........745900
C     DATA IMACH( 6) /          4 /                                      SLAP........746000
C     DATA IMACH( 7) /          2 /                                      SLAP........746100
C     DATA IMACH( 8) /         31 /                                      SLAP........746200
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........746300
C     DATA IMACH(10) /          2 /                                      SLAP........746400
C     DATA IMACH(11) /         24 /                                      SLAP........746500
C     DATA IMACH(12) /       -125 /                                      SLAP........746600
C     DATA IMACH(13) /        128 /                                      SLAP........746700
C     DATA IMACH(14) /         53 /                                      SLAP........746800
C     DATA IMACH(15) /      -1021 /                                      SLAP........746900
C     DATA IMACH(16) /       1024 /                                      SLAP........747000
C                                                                        SLAP........747100
C     MACHINE CONSTANTS FOR THE INTEL i860                               SLAP........747200
C                                                                        SLAP........747300
C     DATA IMACH( 1) /          5 /                                      SLAP........747400
C     DATA IMACH( 2) /          6 /                                      SLAP........747500
C     DATA IMACH( 3) /          6 /                                      SLAP........747600
C     DATA IMACH( 4) /          6 /                                      SLAP........747700
C     DATA IMACH( 5) /         32 /                                      SLAP........747800
C     DATA IMACH( 6) /          4 /                                      SLAP........747900
C     DATA IMACH( 7) /          2 /                                      SLAP........748000
C     DATA IMACH( 8) /         31 /                                      SLAP........748100
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........748200
C     DATA IMACH(10) /          2 /                                      SLAP........748300
C     DATA IMACH(11) /         24 /                                      SLAP........748400
C     DATA IMACH(12) /       -125 /                                      SLAP........748500
C     DATA IMACH(13) /        128 /                                      SLAP........748600
C     DATA IMACH(14) /         53 /                                      SLAP........748700
C     DATA IMACH(15) /      -1021 /                                      SLAP........748800
C     DATA IMACH(16) /       1024 /                                      SLAP........748900
C                                                                        SLAP........749000
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)                    SLAP........749100
C                                                                        SLAP........749200
C     DATA IMACH( 1) /          5 /                                      SLAP........749300
C     DATA IMACH( 2) /          6 /                                      SLAP........749400
C     DATA IMACH( 3) /          5 /                                      SLAP........749500
C     DATA IMACH( 4) /          6 /                                      SLAP........749600
C     DATA IMACH( 5) /         36 /                                      SLAP........749700
C     DATA IMACH( 6) /          5 /                                      SLAP........749800
C     DATA IMACH( 7) /          2 /                                      SLAP........749900
C     DATA IMACH( 8) /         35 /                                      SLAP........750000
C     DATA IMACH( 9) / "377777777777 /                                   SLAP........750100
C     DATA IMACH(10) /          2 /                                      SLAP........750200
C     DATA IMACH(11) /         27 /                                      SLAP........750300
C     DATA IMACH(12) /       -128 /                                      SLAP........750400
C     DATA IMACH(13) /        127 /                                      SLAP........750500
C     DATA IMACH(14) /         54 /                                      SLAP........750600
C     DATA IMACH(15) /       -101 /                                      SLAP........750700
C     DATA IMACH(16) /        127 /                                      SLAP........750800
C                                                                        SLAP........750900
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)                    SLAP........751000
C                                                                        SLAP........751100
C     DATA IMACH( 1) /          5 /                                      SLAP........751200
C     DATA IMACH( 2) /          6 /                                      SLAP........751300
C     DATA IMACH( 3) /          5 /                                      SLAP........751400
C     DATA IMACH( 4) /          6 /                                      SLAP........751500
C     DATA IMACH( 5) /         36 /                                      SLAP........751600
C     DATA IMACH( 6) /          5 /                                      SLAP........751700
C     DATA IMACH( 7) /          2 /                                      SLAP........751800
C     DATA IMACH( 8) /         35 /                                      SLAP........751900
C     DATA IMACH( 9) / "377777777777 /                                   SLAP........752000
C     DATA IMACH(10) /          2 /                                      SLAP........752100
C     DATA IMACH(11) /         27 /                                      SLAP........752200
C     DATA IMACH(12) /       -128 /                                      SLAP........752300
C     DATA IMACH(13) /        127 /                                      SLAP........752400
C     DATA IMACH(14) /         62 /                                      SLAP........752500
C     DATA IMACH(15) /       -128 /                                      SLAP........752600
C     DATA IMACH(16) /        127 /                                      SLAP........752700
C                                                                        SLAP........752800
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                    SLAP........752900
C     32-BIT INTEGER ARITHMETIC.                                         SLAP........753000
C                                                                        SLAP........753100
C     DATA IMACH( 1) /          5 /                                      SLAP........753200
C     DATA IMACH( 2) /          6 /                                      SLAP........753300
C     DATA IMACH( 3) /          5 /                                      SLAP........753400
C     DATA IMACH( 4) /          6 /                                      SLAP........753500
C     DATA IMACH( 5) /         32 /                                      SLAP........753600
C     DATA IMACH( 6) /          4 /                                      SLAP........753700
C     DATA IMACH( 7) /          2 /                                      SLAP........753800
C     DATA IMACH( 8) /         31 /                                      SLAP........753900
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........754000
C     DATA IMACH(10) /          2 /                                      SLAP........754100
C     DATA IMACH(11) /         24 /                                      SLAP........754200
C     DATA IMACH(12) /       -127 /                                      SLAP........754300
C     DATA IMACH(13) /        127 /                                      SLAP........754400
C     DATA IMACH(14) /         56 /                                      SLAP........754500
C     DATA IMACH(15) /       -127 /                                      SLAP........754600
C     DATA IMACH(16) /        127 /                                      SLAP........754700
C                                                                        SLAP........754800
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING                    SLAP........754900
C     16-BIT INTEGER ARITHMETIC.                                         SLAP........755000
C                                                                        SLAP........755100
C     DATA IMACH( 1) /          5 /                                      SLAP........755200
C     DATA IMACH( 2) /          6 /                                      SLAP........755300
C     DATA IMACH( 3) /          5 /                                      SLAP........755400
C     DATA IMACH( 4) /          6 /                                      SLAP........755500
C     DATA IMACH( 5) /         16 /                                      SLAP........755600
C     DATA IMACH( 6) /          2 /                                      SLAP........755700
C     DATA IMACH( 7) /          2 /                                      SLAP........755800
C     DATA IMACH( 8) /         15 /                                      SLAP........755900
C     DATA IMACH( 9) /      32767 /                                      SLAP........756000
C     DATA IMACH(10) /          2 /                                      SLAP........756100
C     DATA IMACH(11) /         24 /                                      SLAP........756200
C     DATA IMACH(12) /       -127 /                                      SLAP........756300
C     DATA IMACH(13) /        127 /                                      SLAP........756400
C     DATA IMACH(14) /         56 /                                      SLAP........756500
C     DATA IMACH(15) /       -127 /                                      SLAP........756600
C     DATA IMACH(16) /        127 /                                      SLAP........756700
C                                                                        SLAP........756800
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS                         SLAP........756900
C                                                                        SLAP........757000
C     DATA IMACH( 1) /          5 /                                      SLAP........757100
C     DATA IMACH( 2) /          6 /                                      SLAP........757200
C     DATA IMACH( 3) /          6 /                                      SLAP........757300
C     DATA IMACH( 4) /          6 /                                      SLAP........757400
C     DATA IMACH( 5) /         32 /                                      SLAP........757500
C     DATA IMACH( 6) /          4 /                                      SLAP........757600
C     DATA IMACH( 7) /          2 /                                      SLAP........757700
C     DATA IMACH( 8) /         31 /                                      SLAP........757800
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........757900
C     DATA IMACH(10) /          2 /                                      SLAP........758000
C     DATA IMACH(11) /         24 /                                      SLAP........758100
C     DATA IMACH(12) /       -125 /                                      SLAP........758200
C     DATA IMACH(13) /        128 /                                      SLAP........758300
C     DATA IMACH(14) /         53 /                                      SLAP........758400
C     DATA IMACH(15) /      -1021 /                                      SLAP........758500
C     DATA IMACH(16) /       1024 /                                      SLAP........758600
C                                                                        SLAP........758700
C     MACHINE CONSTANTS FOR THE SUN                                      SLAP........758800
C                                                                        SLAP........758900
C     DATA IMACH( 1) /          5 /                                      SLAP........759000
C     DATA IMACH( 2) /          6 /                                      SLAP........759100
C     DATA IMACH( 3) /          6 /                                      SLAP........759200
C     DATA IMACH( 4) /          6 /                                      SLAP........759300
C     DATA IMACH( 5) /         32 /                                      SLAP........759400
C     DATA IMACH( 6) /          4 /                                      SLAP........759500
C     DATA IMACH( 7) /          2 /                                      SLAP........759600
C     DATA IMACH( 8) /         31 /                                      SLAP........759700
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........759800
C     DATA IMACH(10) /          2 /                                      SLAP........759900
C     DATA IMACH(11) /         24 /                                      SLAP........760000
C     DATA IMACH(12) /       -125 /                                      SLAP........760100
C     DATA IMACH(13) /        128 /                                      SLAP........760200
C     DATA IMACH(14) /         53 /                                      SLAP........760300
C     DATA IMACH(15) /      -1021 /                                      SLAP........760400
C     DATA IMACH(16) /       1024 /                                      SLAP........760500
C                                                                        SLAP........760600
C     MACHINE CONSTANTS FOR THE SUN                                      SLAP........760700
C     USING THE -r8 COMPILER OPTION                                      SLAP........760800
C                                                                        SLAP........760900
C     DATA IMACH( 1) /          5 /                                      SLAP........761000
C     DATA IMACH( 2) /          6 /                                      SLAP........761100
C     DATA IMACH( 3) /          6 /                                      SLAP........761200
C     DATA IMACH( 4) /          6 /                                      SLAP........761300
C     DATA IMACH( 5) /         32 /                                      SLAP........761400
C     DATA IMACH( 6) /          4 /                                      SLAP........761500
C     DATA IMACH( 7) /          2 /                                      SLAP........761600
C     DATA IMACH( 8) /         31 /                                      SLAP........761700
C     DATA IMACH( 9) / 2147483647 /                                      SLAP........761800
C     DATA IMACH(10) /          2 /                                      SLAP........761900
C     DATA IMACH(11) /         53 /                                      SLAP........762000
C     DATA IMACH(12) /      -1021 /                                      SLAP........762100
C     DATA IMACH(13) /       1024 /                                      SLAP........762200
C     DATA IMACH(14) /        113 /                                      SLAP........762300
C     DATA IMACH(15) /     -16381 /                                      SLAP........762400
C     DATA IMACH(16) /      16384 /                                      SLAP........762500
C                                                                        SLAP........762600
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER          SLAP........762700
C                                                                        SLAP........762800
C     DATA IMACH( 1) /          5 /                                      SLAP........762900
C     DATA IMACH( 2) /          6 /                                      SLAP........763000
C     DATA IMACH( 3) /          1 /                                      SLAP........763100
C     DATA IMACH( 4) /          6 /                                      SLAP........763200
C     DATA IMACH( 5) /         36 /                                      SLAP........763300
C     DATA IMACH( 6) /          4 /                                      SLAP........763400
C     DATA IMACH( 7) /          2 /                                      SLAP........763500
C     DATA IMACH( 8) /         35 /                                      SLAP........763600
C     DATA IMACH( 9) / O377777777777 /                                   SLAP........763700
C     DATA IMACH(10) /          2 /                                      SLAP........763800
C     DATA IMACH(11) /         27 /                                      SLAP........763900
C     DATA IMACH(12) /       -128 /                                      SLAP........764000
C     DATA IMACH(13) /        127 /                                      SLAP........764100
C     DATA IMACH(14) /         60 /                                      SLAP........764200
C     DATA IMACH(15) /      -1024 /                                      SLAP........764300
C     DATA IMACH(16) /       1023 /                                      SLAP........764400
C                                                                        SLAP........764500
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR                       SLAP........764600
C                                                                        SLAP........764700
C     DATA IMACH( 1) /          1 /                                      SLAP........764800
C     DATA IMACH( 2) /          1 /                                      SLAP........764900
C     DATA IMACH( 3) /          0 /                                      SLAP........765000
C     DATA IMACH( 4) /          1 /                                      SLAP........765100
C     DATA IMACH( 5) /         16 /                                      SLAP........765200
C     DATA IMACH( 6) /          2 /                                      SLAP........765300
C     DATA IMACH( 7) /          2 /                                      SLAP........765400
C     DATA IMACH( 8) /         15 /                                      SLAP........765500
C     DATA IMACH( 9) /      32767 /                                      SLAP........765600
C     DATA IMACH(10) /          2 /                                      SLAP........765700
C     DATA IMACH(11) /         24 /                                      SLAP........765800
C     DATA IMACH(12) /       -127 /                                      SLAP........765900
C     DATA IMACH(13) /        127 /                                      SLAP........766000
C     DATA IMACH(14) /         56 /                                      SLAP........766100
C     DATA IMACH(15) /       -127 /                                      SLAP........766200
C     DATA IMACH(16) /        127 /                                      SLAP........766300
C                                                                        SLAP........766400
C.....GENERAL-PURPOSE VALUES INSERTED DURING INTEGRATION OF SLAP         SLAP........766500
C        WITH SUTRA.  (ONLY SPECIFIED VALUES USED BY SLAP ROUTINES.)     SLAP........766600
C                                                                        SLAP........766700
      DATA IMACH( 4) /          6 /                                      SLAP........766800
C                                                                        SLAP........766900
C***FIRST EXECUTABLE STATEMENT  I1MACH                                   SLAP........767000
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10                            SLAP........767100
C                                                                        SLAP........767200
      I1MACH = IMACH(I)                                                  SLAP........767300
      RETURN                                                             SLAP........767400
C                                                                        SLAP........767500
   10 CONTINUE                                                           SLAP........767600
      WRITE (UNIT = OUTPUT, FMT = 9000)                                  SLAP........767700
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')                 SLAP........767800
C                                                                        SLAP........767900
C     CALL FDUMP                                                         SLAP........768000
C                                                                        SLAP........768100
      STOP                                                               SLAP........768200
      END                                                                SLAP........768300
*DECK ISDCG                                                              SLAP........768400
      INTEGER FUNCTION ISDCG (N, B, X, NELT, IA, JA, A, ISYM, MSOLVE,    SLAP........768500
     +   ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,   SLAP........768600
     +   IWORK, AK, BK, BNRM, SOLNRM)                                    SLAP........768700
C***BEGIN PROLOGUE  ISDCG                                                SLAP........768800
C***SUBSIDIARY                                                           SLAP........768900
C***PURPOSE  Preconditioned Conjugate Gradient Stop Test.                SLAP........769000
C            This routine calculates the stop test for the Conjugate     SLAP........769100
C            Gradient iteration scheme.  It returns a non-zero if the    SLAP........769200
C            error estimate (the type of which is determined by ITOL)    SLAP........769300
C            is less than the user specified tolerance TOL.              SLAP........769400
C***LIBRARY   SLATEC (SLAP)                                              SLAP........769500
C***CATEGORY  D2B4                                                       SLAP........769600
C***TYPE      DOUBLE PRECISION (ISSCG-S, ISDCG-D)                        SLAP........769700
C***KEYWORDS  LINEAR SYSTEM, SLAP, SPARSE, STOP TEST                     SLAP........769800
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........769900
C           Seager, Mark K., (LLNL)                                      SLAP........770000
C             Lawrence Livermore National Laboratory                     SLAP........770100
C             PO BOX 808, L-60                                           SLAP........770200
C             Livermore, CA 94550 (510) 423-3141                         SLAP........770300
C             seager@llnl.gov                                            SLAP........770400
C***DESCRIPTION                                                          SLAP........770500
C                                                                        SLAP........770600
C *Usage:                                                                SLAP........770700
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER       SLAP........770800
C     INTEGER IERR, IUNIT, IWORK(USER DEFINED)                           SLAP........770900
C     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N)            SLAP........771000
C     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINED), AK, BK          SLAP........771100
C     DOUBLE PRECISION BNRM, SOLNRM                                      SLAP........771200
C     EXTERNAL MSOLVE                                                    SLAP........771300
C                                                                        SLAP........771400
C     IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,       SLAP........771500
C    $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK, IWORK,     SLAP........771600
C    $     AK, BK, BNRM, SOLNRM) .NE. 0 ) THEN ITERATION DONE            SLAP........771700
C                                                                        SLAP........771800
C *Arguments:                                                            SLAP........771900
C N      :IN       Integer.                                              SLAP........772000
C         Order of the Matrix.                                           SLAP........772100
C B      :IN       Double Precision B(N).                                SLAP........772200
C         Right-hand side vector.                                        SLAP........772300
C X      :IN       Double Precision X(N).                                SLAP........772400
C         The current approximate solution vector.                       SLAP........772500
C NELT   :IN       Integer.                                              SLAP........772600
C         Number of Non-Zeros stored in A.                               SLAP........772700
C IA     :IN       Integer IA(NELT).                                     SLAP........772800
C JA     :IN       Integer JA(NELT).                                     SLAP........772900
C A      :IN       Double Precision A(NELT).                             SLAP........773000
C         These arrays should hold the matrix A in either the SLAP       SLAP........773100
C         Triad format or the SLAP Column format.  See "Description"     SLAP........773200
C         in the DCG, DSDCG or DSICCG routines.                          SLAP........773300
C ISYM   :IN       Integer.                                              SLAP........773400
C         Flag to indicate symmetric storage format.                     SLAP........773500
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........773600
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........773700
C         or lower triangle of the matrix is stored.                     SLAP........773800
C MSOLVE :EXT      External.                                             SLAP........773900
C         Name of a routine which solves a linear system MZ = R for      SLAP........774000
C         Z given R with the preconditioning matrix M (M is supplied via SLAP........774100
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  SLAP........774200
C         be declared external in the calling program.  The calling      SLAP........774300
C         sequence to MSOLVE is:                                         SLAP........774400
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........774500
C         Where N is the number of unknowns, R is the right-hand side    SLAP........774600
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........774700
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........774800
C         that can be used to pass necessary preconditioning information SLAP........774900
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........775000
C         for the same purpose as RWORK.                                 SLAP........775100
C ITOL   :IN       Integer.                                              SLAP........775200
C         Flag to indicate type of convergence criterion.                SLAP........775300
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........775400
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........775500
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........775600
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........775700
C         side is less than TOL, where M-inv is the inverse of the       SLAP........775800
C         diagonal of A.                                                 SLAP........775900
C         ITOL=11 is often useful for checking and comparing different   SLAP........776000
C         routines.  For this case, the user must supply the "exact"     SLAP........776100
C         solution or a very accurate approximation (one with an error   SLAP........776200
C         much less than TOL) through a common block,                    SLAP........776300
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........776400
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........776500
C         between the iterative approximation and the user-supplied      SLAP........776600
C         solution divided by the 2-norm of the user-supplied solution   SLAP........776700
C         is less than TOL.  Note that this requires the user to set up  SLAP........776800
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........776900
C         The routine with this declaration should be loaded before the  SLAP........777000
C         stop test so that the correct length is used by the loader.    SLAP........777100
C         This procedure is not standard Fortran and may not work        SLAP........777200
C         correctly on your system (although it has worked on every      SLAP........777300
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........777400
C         common block is indeed standard Fortran.                       SLAP........777500
C TOL    :IN       Double Precision.                                     SLAP........777600
C         Convergence criterion, as described above.                     SLAP........777700
C ITMAX  :IN       Integer.                                              SLAP........777800
C         Maximum number of iterations.                                  SLAP........777900
C ITER   :IN       Integer.                                              SLAP........778000
C         Current iteration count.  (Must be zero on first call.)        SLAP........778100
C ERR    :OUT      Double Precision.                                     SLAP........778200
C         Error estimate of error in the X(N) approximate solution, as   SLAP........778300
C         defined by ITOL.                                               SLAP........778400
C IERR   :OUT      Integer.                                              SLAP........778500
C         Error flag.  IERR is set to 3 if ITOL is not one of the        SLAP........778600
C         acceptable values, see above.                                  SLAP........778700
C IUNIT  :IN       Integer.                                              SLAP........778800
C         Unit number on which to write the error at each iteration,     SLAP........778900
C         if this is desired for monitoring convergence.  If unit        SLAP........779000
C         number is 0, no writing will occur.                            SLAP........779100
C R      :IN       Double Precision R(N).                                SLAP........779200
C         The residual R = B-AX.                                         SLAP........779300
C Z      :WORK     Double Precision Z(N).                                SLAP........779400
C         Workspace used to hold the pseudo-residual M Z = R.            SLAP........779500
C P      :IN       Double Precision P(N).                                SLAP........779600
C         The conjugate direction vector.                                SLAP........779700
C DZ     :WORK     Double Precision DZ(N).                               SLAP........779800
C         Workspace used to hold temporary vector(s).                    SLAP........779900
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........780000
C         Double Precision array that can be used by MSOLVE.             SLAP........780100
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........780200
C         Integer array that can be used by MSOLVE.                      SLAP........780300
C AK     :IN       Double Precision.                                     SLAP........780400
C BK     :IN       Double Precision.                                     SLAP........780500
C         Current conjugate gradient parameters alpha and beta.          SLAP........780600
C BNRM   :INOUT    Double Precision.                                     SLAP........780700
C         Norm of the right hand side.  Type of norm depends on ITOL.    SLAP........780800
C         Calculated only on the first call.                             SLAP........780900
C SOLNRM :INOUT    Double Precision.                                     SLAP........781000
C         2-Norm of the true solution, SOLN.  Only computed and used     SLAP........781100
C         if ITOL = 11.                                                  SLAP........781200
C                                                                        SLAP........781300
C *Function Return Values:                                               SLAP........781400
C       0 : Error estimate (determined by ITOL) is *NOT* less than the   SLAP........781500
C           specified tolerance, TOL.  The iteration must continue.      SLAP........781600
C       1 : Error estimate (determined by ITOL) is less than the         SLAP........781700
C           specified tolerance, TOL.  The iteration can be considered   SLAP........781800
C           complete.                                                    SLAP........781900
C                                                                        SLAP........782000
C *Cautions:                                                             SLAP........782100
C     This routine will attempt to write to the Fortran logical output   SLAP........782200
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........782300
C     this logical unit is attached to a file or terminal before calling SLAP........782400
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........782500
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........782600
C                                                                        SLAP........782700
C***SEE ALSO  DCG, DSDCG, DSICCG                                         SLAP........782800
C***ROUTINES CALLED  D1MACH, DNRM2                                       SLAP........782900
C***COMMON BLOCKS    DSLBLK                                              SLAP........783000
C***REVISION HISTORY  (YYMMDD)                                           SLAP........783100
C   890404  DATE WRITTEN                                                 SLAP........783200
C   890404  Previous REVISION DATE                                       SLAP........783300
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........783400
C   890921  Removed TeX from comments.  (FNF)                            SLAP........783500
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........783600
C           standard.  (FNF)                                             SLAP........783700
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........783800
C   891003  Removed C***REFER TO line, per MKS.                          SLAP........783900
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........784000
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)             SLAP........784100
C   910506  Made subsidiary to DCG.  (FNF)                               SLAP........784200
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........784300
C   920511  Added complete declaration section.  (WRB)                   SLAP........784400
C   920930  Corrected to not print AK,BK when ITER=0.  (FNF)             SLAP........784500
C   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in          SLAP........784600
C           output format.  (FNF)                                        SLAP........784700
C***END PROLOGUE  ISDCG                                                  SLAP........784800
C     .. Scalar Arguments ..                                             SLAP........784900
      DOUBLE PRECISION AK, BK, BNRM, ERR, SOLNRM, TOL                    SLAP........785000
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT              SLAP........785100
C     .. Array Arguments ..                                              SLAP........785200
      DOUBLE PRECISION A(NELT), B(N), DZ(N), P(N), R(N), RWORK(*), X(N), SLAP........785300
     +                 Z(N)                                              SLAP........785400
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........785500
C     .. Subroutine Arguments ..                                         SLAP........785600
      EXTERNAL MSOLVE                                                    SLAP........785700
C     .. Arrays in Common ..                                             SLAP........785800
      DOUBLE PRECISION SOLN(1)                                           SLAP........785900
C     .. Local Scalars ..                                                SLAP........786000
      INTEGER I                                                          SLAP........786100
C     .. External Functions ..                                           SLAP........786200
      DOUBLE PRECISION D1MACH, DNRM2                                     SLAP........786300
      EXTERNAL D1MACH, DNRM2                                             SLAP........786400
C     .. Common blocks ..                                                SLAP........786500
      COMMON /DSLBLK/ SOLN                                               SLAP........786600
C***FIRST EXECUTABLE STATEMENT  ISDCG                                    SLAP........786700
      ISDCG = 0                                                          SLAP........786800
C                                                                        SLAP........786900
      IF( ITOL.EQ.1 ) THEN                                               SLAP........787000
C         err = ||Residual||/||RightHandSide|| (2-Norms).                SLAP........787100
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)                           SLAP........787200
         ERR = DNRM2(N, R, 1)/BNRM                                       SLAP........787300
      ELSE IF( ITOL.EQ.2 ) THEN                                          SLAP........787400
C                  -1              -1                                    SLAP........787500
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).          SLAP........787600
         IF(ITER .EQ. 0) THEN                                            SLAP........787700
C...........THE NEXT LINE OF CODE REFLECTS A BUG FIX MADE DURING         SLAP........787800
C              INTEGRATION OF SLAP WITH SUTRA.  IT ORIGINALLY READ:      SLAP........787900
C            CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)   SLAP........788000
            CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)   SLAP........788100
            BNRM = DNRM2(N, DZ, 1)                                       SLAP........788200
         ENDIF                                                           SLAP........788300
         ERR = DNRM2(N, Z, 1)/BNRM                                       SLAP........788400
      ELSE IF( ITOL.EQ.11 ) THEN                                         SLAP........788500
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).           SLAP........788600
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)                      SLAP........788700
         DO 10 I = 1, N                                                  SLAP........788800
            DZ(I) = X(I) - SOLN(I)                                       SLAP........788900
 10      CONTINUE                                                        SLAP........789000
         ERR = DNRM2(N, DZ, 1)/SOLNRM                                    SLAP........789100
      ELSE                                                               SLAP........789200
C                                                                        SLAP........789300
C         If we get here ITOL is not one of the acceptable values.       SLAP........789400
         ERR = D1MACH(2)                                                 SLAP........789500
         IERR = 3                                                        SLAP........789600
      ENDIF                                                              SLAP........789700
C                                                                        SLAP........789800
      IF(IUNIT .NE. 0) THEN                                              SLAP........789900
         IF( ITER.EQ.0 ) THEN                                            SLAP........790000
C...........THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION        SLAP........790100
C              OF SLAP WITH SUTRA.  IT ORIGINALLY READ:                  SLAP........790200
C              WRITE(IUNIT,1000) N, ITOL                                 SLAP........790300
            WRITE(IUNIT,1000)                                            SLAP........790400
            WRITE(IUNIT,1010) ITER, ERR                                  SLAP........790500
         ELSE                                                            SLAP........790600
C...........THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION        SLAP........790700
C              OF SLAP WITH SUTRA.  IT ORIGINALLY READ:                  SLAP........790800
C              WRITE(IUNIT,1010) ITER, ERR, AK, BK                       SLAP........790900
            WRITE(IUNIT,1010) ITER, ERR                                  SLAP........791000
         ENDIF                                                           SLAP........791100
      ENDIF                                                              SLAP........791200
      IF(ERR .LE. TOL) ISDCG = 1                                         SLAP........791300
      RETURN                                                             SLAP........791400
C.....THE NEXT TWO LINES OF CODE REFLECT MODIFICATIONS MADE DURING       SLAP........791500
C        INTEGRATION OF SLAP WITH SUTRA.  THEY ORIGINALLY READ:          SLAP........791600
C        1000 FORMAT(' Preconditioned Conjugate Gradient for ',          SLAP........791700
C            $     'N, ITOL = ',I5, I5,                                  SLAP........791800
C            $     /' ITER','   Error Estimate','            Alpha',     SLAP........791900
C            $     '             Beta')                                  SLAP........792000
C        1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7,1X,D16.7)                   SLAP........792100
1000  FORMAT(1X,6X,'Solver Iteration','   Error Estimate')               SLAP........792200
1010  FORMAT(1X,6X,I16,1X,1PE16.7)                                       SLAP........792300
C------------- LAST LINE OF ISDCG FOLLOWS ------------------------------ SLAP........792400
      END                                                                SLAP........792500
*DECK ISDGMR                                                             SLAP........792600
      INTEGER FUNCTION ISDGMR (N, B, X, XL, NELT, IA, JA, A, ISYM,       SLAP........792700
     +   MSOLVE, NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,     SLAP........792800
     +   RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL,       SLAP........792900
     +   MAXLP1, V, Q, SNORMW, PROD, R0NRM, HES, JPRE)                   SLAP........793000
C***BEGIN PROLOGUE  ISDGMR                                               SLAP........793100
C***SUBSIDIARY                                                           SLAP........793200
C***PURPOSE  Generalized Minimum Residual Stop Test.                     SLAP........793300
C            This routine calculates the stop test for the Generalized   SLAP........793400
C            Minimum RESidual (GMRES) iteration scheme.  It returns a    SLAP........793500
C            non-zero if the error estimate (the type of which is        SLAP........793600
C            determined by ITOL) is less than the user specified         SLAP........793700
C            tolerance TOL.                                              SLAP........793800
C***LIBRARY   SLATEC (SLAP)                                              SLAP........793900
C***CATEGORY  D2A4, D2B4                                                 SLAP........794000
C***TYPE      DOUBLE PRECISION (ISSGMR-S, ISDGMR-D)                      SLAP........794100
C***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST              SLAP........794200
C***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov                       SLAP........794300
C           Hindmarsh, Alan, (LLNL), alanh@llnl.gov                      SLAP........794400
C           Seager, Mark K., (LLNL), seager@llnl.gov                     SLAP........794500
C             Lawrence Livermore National Laboratory                     SLAP........794600
C             PO Box 808, L-60                                           SLAP........794700
C             Livermore, CA 94550 (510) 423-3141                         SLAP........794800
C***DESCRIPTION                                                          SLAP........794900
C                                                                        SLAP........795000
C *Usage:                                                                SLAP........795100
C      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL             SLAP........795200
C      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL            SLAP........795300
C      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE                             SLAP........795400
C      DOUBLE PRECISION B(N), X(N), XL(MAXL), A(NELT), TOL, ERR,         SLAP........795500
C     $                 R(N), Z(N), DZ(N), RWORK(USER DEFINED),          SLAP........795600
C     $                 RNRM, BNRM, SB(N), SX(N), V(N,MAXLP1),           SLAP........795700
C     $                 Q(2*MAXL), SNORMW, PROD, R0NRM,                  SLAP........795800
C     $                 HES(MAXLP1,MAXL)                                 SLAP........795900
C      EXTERNAL MSOLVE                                                   SLAP........796000
C                                                                        SLAP........796100
C      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,            SLAP........796200
C     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,          SLAP........796300
C     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,                     SLAP........796400
C     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,          SLAP........796500
C     $     HES, JPRE) .NE. 0) THEN ITERATION DONE                       SLAP........796600
C                                                                        SLAP........796700
C *Arguments:                                                            SLAP........796800
C N      :IN       Integer.                                              SLAP........796900
C         Order of the Matrix.                                           SLAP........797000
C B      :IN       Double Precision B(N).                                SLAP........797100
C         Right-hand-side vector.                                        SLAP........797200
C X      :IN       Double Precision X(N).                                SLAP........797300
C         Approximate solution vector as of the last restart.            SLAP........797400
C XL     :OUT      Double Precision XL(N)                                SLAP........797500
C         An array of length N used to hold the approximate              SLAP........797600
C         solution as of the current iteration.  Only computed by        SLAP........797700
C         this routine when ITOL=11.                                     SLAP........797800
C NELT   :IN       Integer.                                              SLAP........797900
C         Number of Non-Zeros stored in A.                               SLAP........798000
C IA     :IN       Integer IA(NELT).                                     SLAP........798100
C JA     :IN       Integer JA(NELT).                                     SLAP........798200
C A      :IN       Double Precision A(NELT).                             SLAP........798300
C         These arrays contain the matrix data structure for A.          SLAP........798400
C         It could take any form.  See "Description", in the DGMRES,     SLAP........798500
C         DSLUGM and DSDGMR routines for more details.                   SLAP........798600
C ISYM   :IN       Integer.                                              SLAP........798700
C         Flag to indicate symmetric storage format.                     SLAP........798800
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........798900
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........799000
C         or lower triangle of the matrix is stored.                     SLAP........799100
C MSOLVE :EXT      External.                                             SLAP........799200
C         Name of a routine which solves a linear system Mz = r for  z   SLAP........799300
C         given r with the preconditioning matrix M (M is supplied via   SLAP........799400
C         RWORK and IWORK arrays.  The name of the MSOLVE routine must   SLAP........799500
C         be declared external in the calling program.  The calling      SLAP........799600
C         sequence to MSOLVE is:                                         SLAP........799700
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........799800
C         Where N is the number of unknowns, R is the right-hand side    SLAP........799900
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........800000
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........800100
C         that can be used to pass necessary preconditioning information SLAP........800200
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........800300
C         for the same purpose as RWORK.                                 SLAP........800400
C NMSL   :INOUT    Integer.                                              SLAP........800500
C         A counter for the number of calls to MSOLVE.                   SLAP........800600
C ITOL   :IN       Integer.                                              SLAP........800700
C         Flag to indicate the type of convergence criterion used.       SLAP........800800
C         ITOL=0  Means the  iteration stops when the test described     SLAP........800900
C                 below on  the  residual RL  is satisfied.  This is     SLAP........801000
C                 the  "Natural Stopping Criteria" for this routine.     SLAP........801100
C                 Other values  of   ITOL  cause  extra,   otherwise     SLAP........801200
C                 unnecessary, computation per iteration and     are     SLAP........801300
C                 therefore much less efficient.                         SLAP........801400
C         ITOL=1  Means   the  iteration stops   when the first test     SLAP........801500
C                 described below on  the residual RL  is satisfied,     SLAP........801600
C                 and there  is either right  or  no preconditioning     SLAP........801700
C                 being used.                                            SLAP........801800
C         ITOL=2  Implies     that   the  user    is   using    left     SLAP........801900
C                 preconditioning, and the second stopping criterion     SLAP........802000
C                 below is used.                                         SLAP........802100
C         ITOL=3  Means the  iteration stops   when  the  third test     SLAP........802200
C                 described below on Minv*Residual is satisfied, and     SLAP........802300
C                 there is either left  or no  preconditioning begin     SLAP........802400
C                 used.                                                  SLAP........802500
C         ITOL=11 is    often  useful  for   checking  and comparing     SLAP........802600
C                 different routines.  For this case, the  user must     SLAP........802700
C                 supply  the  "exact" solution or  a  very accurate     SLAP........802800
C                 approximation (one with  an  error much less  than     SLAP........802900
C                 TOL) through a common block,                           SLAP........803000
C                     COMMON /DSLBLK/ SOLN( )                            SLAP........803100
C                 If ITOL=11, iteration stops when the 2-norm of the     SLAP........803200
C                 difference between the iterative approximation and     SLAP........803300
C                 the user-supplied solution  divided by the  2-norm     SLAP........803400
C                 of the  user-supplied solution  is  less than TOL.     SLAP........803500
C                 Note that this requires  the  user to  set up  the     SLAP........803600
C                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling     SLAP........803700
C                 routine.  The routine with this declaration should     SLAP........803800
C                 be loaded before the stop test so that the correct     SLAP........803900
C                 length is used by  the loader.  This procedure  is     SLAP........804000
C                 not standard Fortran and may not work correctly on     SLAP........804100
C                 your   system (although  it  has  worked  on every     SLAP........804200
C                 system the authors have tried).  If ITOL is not 11     SLAP........804300
C                 then this common block is indeed standard Fortran.     SLAP........804400
C TOL    :IN       Double Precision.                                     SLAP........804500
C         Convergence criterion, as described above.                     SLAP........804600
C ITMAX  :IN       Integer.                                              SLAP........804700
C         Maximum number of iterations.                                  SLAP........804800
C ITER   :IN       Integer.                                              SLAP........804900
C         The iteration for which to check for convergence.              SLAP........805000
C ERR    :OUT      Double Precision.                                     SLAP........805100
C         Error estimate of error in final approximate solution, as      SLAP........805200
C         defined by ITOL.  Letting norm() denote the Euclidean          SLAP........805300
C         norm, ERR is defined as follows..                              SLAP........805400
C                                                                        SLAP........805500
C         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........805600
C                               for right or no preconditioning, and     SLAP........805700
C                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........805800
C                                norm(SB*(M-inverse)*B),                 SLAP........805900
C                               for left preconditioning.                SLAP........806000
C         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),          SLAP........806100
C                               since right or no preconditioning        SLAP........806200
C                               being used.                              SLAP........806300
C         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/         SLAP........806400
C                                norm(SB*(M-inverse)*B),                 SLAP........806500
C                               since left preconditioning is being      SLAP........806600
C                               used.                                    SLAP........806700
C         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|        SLAP........806800
C                               i=1,n                                    SLAP........806900
C         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).     SLAP........807000
C IUNIT  :IN       Integer.                                              SLAP........807100
C         Unit number on which to write the error at each iteration,     SLAP........807200
C         if this is desired for monitoring convergence.  If unit        SLAP........807300
C         number is 0, no writing will occur.                            SLAP........807400
C R      :INOUT    Double Precision R(N).                                SLAP........807500
C         Work array used in calling routine.  It contains               SLAP........807600
C         information necessary to compute the residual RL = B-A*XL.     SLAP........807700
C Z      :WORK     Double Precision Z(N).                                SLAP........807800
C         Workspace used to hold the pseudo-residual M z = r.            SLAP........807900
C DZ     :WORK     Double Precision DZ(N).                               SLAP........808000
C         Workspace used to hold temporary vector(s).                    SLAP........808100
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........808200
C         Double Precision array that can be used by MSOLVE.             SLAP........808300
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........808400
C         Integer array that can be used by MSOLVE.                      SLAP........808500
C RNRM   :IN       Double Precision.                                     SLAP........808600
C         Norm of the current residual.  Type of norm depends on ITOL.   SLAP........808700
C BNRM   :IN       Double Precision.                                     SLAP........808800
C         Norm of the right hand side.  Type of norm depends on ITOL.    SLAP........808900
C SB     :IN       Double Precision SB(N).                               SLAP........809000
C         Scaling vector for B.                                          SLAP........809100
C SX     :IN       Double Precision SX(N).                               SLAP........809200
C         Scaling vector for X.                                          SLAP........809300
C JSCAL  :IN       Integer.                                              SLAP........809400
C         Flag indicating if scaling arrays SB and SX are being          SLAP........809500
C         used in the calling routine DPIGMR.                            SLAP........809600
C         JSCAL=0 means SB and SX are not used and the                   SLAP........809700
C                 algorithm will perform as if all                       SLAP........809800
C                 SB(i) = 1 and SX(i) = 1.                               SLAP........809900
C         JSCAL=1 means only SX is used, and the algorithm               SLAP........810000
C                 performs as if all SB(i) = 1.                          SLAP........810100
C         JSCAL=2 means only SB is used, and the algorithm               SLAP........810200
C                 performs as if all SX(i) = 1.                          SLAP........810300
C         JSCAL=3 means both SB and SX are used.                         SLAP........810400
C KMP    :IN       Integer                                               SLAP........810500
C         The number of previous vectors the new vector VNEW             SLAP........810600
C         must be made orthogonal to.  (KMP .le. MAXL)                   SLAP........810700
C LGMR   :IN       Integer                                               SLAP........810800
C         The number of GMRES iterations performed on the current call   SLAP........810900
C         to DPIGMR (i.e., # iterations since the last restart) and      SLAP........811000
C         the current order of the upper Hessenberg                      SLAP........811100
C         matrix HES.                                                    SLAP........811200
C MAXL   :IN       Integer                                               SLAP........811300
C         The maximum allowable order of the matrix H.                   SLAP........811400
C MAXLP1 :IN       Integer                                               SLAP........811500
C         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.       SLAP........811600
C V      :IN       Double Precision V(N,MAXLP1)                          SLAP........811700
C         The N by (LGMR+1) array containing the LGMR                    SLAP........811800
C         orthogonal vectors V(*,1) to V(*,LGMR).                        SLAP........811900
C Q      :IN       Double Precision Q(2*MAXL)                            SLAP........812000
C         A double precision array of length 2*MAXL containing the       SLAP........812100
C         components of the Givens rotations used in the QR              SLAP........812200
C         decomposition of HES.                                          SLAP........812300
C SNORMW :IN       Double Precision                                      SLAP........812400
C         A scalar containing the scaled norm of VNEW before it          SLAP........812500
C         is renormalized in DPIGMR.                                     SLAP........812600
C PROD   :IN       Double Precision                                      SLAP........812700
C         The product s1*s2*...*sl = the product of the sines of the     SLAP........812800
C         Givens rotations used in the QR factorization of the           SLAP........812900
C         Hessenberg matrix HES.                                         SLAP........813000
C R0NRM  :IN       Double Precision                                      SLAP........813100
C         The scaled norm of initial residual R0.                        SLAP........813200
C HES    :IN       Double Precision HES(MAXLP1,MAXL)                     SLAP........813300
C         The upper triangular factor of the QR decomposition            SLAP........813400
C         of the (LGMR+1) by LGMR upper Hessenberg matrix whose          SLAP........813500
C         entries are the scaled inner-products of A*V(*,I)              SLAP........813600
C         and V(*,K).                                                    SLAP........813700
C JPRE   :IN       Integer                                               SLAP........813800
C         Preconditioner type flag.                                      SLAP........813900
C         (See description of IGWK(4) in DGMRES.)                        SLAP........814000
C                                                                        SLAP........814100
C *Description                                                           SLAP........814200
C       When using the GMRES solver,  the preferred value  for ITOL      SLAP........814300
C       is 0.  This is due to the fact that when ITOL=0 the norm of      SLAP........814400
C       the residual required in the stopping test is  obtained for      SLAP........814500
C       free, since this value is already  calculated  in the GMRES      SLAP........814600
C       algorithm.   The  variable  RNRM contains the   appropriate      SLAP........814700
C       norm, which is equal to norm(SB*(RL - A*XL))  when right or      SLAP........814800
C       no   preconditioning is  being  performed,   and equal   to      SLAP........814900
C       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.      SLAP........815000
C       Here, norm() is the Euclidean norm.  Nonzero values of ITOL      SLAP........815100
C       require  additional work  to  calculate the  actual  scaled      SLAP........815200
C       residual  or its scaled/preconditioned  form,  and/or   the      SLAP........815300
C       approximate solution XL.  Hence, these values of  ITOL will      SLAP........815400
C       not be as efficient as ITOL=0.                                   SLAP........815500
C                                                                        SLAP........815600
C *Cautions:                                                             SLAP........815700
C     This routine will attempt to write to the Fortran logical output   SLAP........815800
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........815900
C     this logical unit is attached to a file or terminal before calling SLAP........816000
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........816100
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........816200
C                                                                        SLAP........816300
C     This routine does not verify that ITOL has a valid value.          SLAP........816400
C     The calling routine should make such a test before calling         SLAP........816500
C     ISDGMR, as is done in DGMRES.                                      SLAP........816600
C                                                                        SLAP........816700
C***SEE ALSO  DGMRES                                                     SLAP........816800
C***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DRLCAL, DSCAL, DXLCAL         SLAP........816900
C***COMMON BLOCKS    DSLBLK                                              SLAP........817000
C***REVISION HISTORY  (YYMMDD)                                           SLAP........817100
C   890404  DATE WRITTEN                                                 SLAP........817200
C   890404  Previous REVISION DATE                                       SLAP........817300
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........817400
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........817500
C           standard.  (FNF)                                             SLAP........817600
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........817700
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........817800
C   910502  Corrected conversion errors, etc.  (FNF)                     SLAP........817900
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)             SLAP........818000
C   910506  Made subsidiary to DGMRES.  (FNF)                            SLAP........818100
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........818200
C   920511  Added complete declaration section.  (WRB)                   SLAP........818300
C   921026  Corrected D to E in output format.  (FNF)                    SLAP........818400
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........818500
C***END PROLOGUE  ISDGMR                                                 SLAP........818600
C     .. Scalar Arguments ..                                             SLAP........818700
      DOUBLE PRECISION BNRM, ERR, PROD, R0NRM, RNRM, SNORMW, TOL         SLAP........818800
      INTEGER ISYM, ITER, ITMAX, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR,    SLAP........818900
     +        MAXL, MAXLP1, N, NELT, NMSL                                SLAP........819000
C     .. Array Arguments ..                                              SLAP........819100
      DOUBLE PRECISION A(*), B(*), DZ(*), HES(MAXLP1, MAXL), Q(*), R(*), SLAP........819200
     +                 RWORK(*), SB(*), SX(*), V(N,*), X(*), XL(*), Z(*) SLAP........819300
      INTEGER IA(*), IWORK(*), JA(*)                                     SLAP........819400
C     .. Subroutine Arguments ..                                         SLAP........819500
      EXTERNAL MSOLVE                                                    SLAP........819600
C     .. Arrays in Common ..                                             SLAP........819700
      DOUBLE PRECISION SOLN(1)                                           SLAP........819800
C     .. Local Scalars ..                                                SLAP........819900
      DOUBLE PRECISION DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM             SLAP........820000
      INTEGER I, IELMAX                                                  SLAP........820100
C     .. External Functions ..                                           SLAP........820200
      DOUBLE PRECISION D1MACH, DNRM2                                     SLAP........820300
      EXTERNAL D1MACH, DNRM2                                             SLAP........820400
C     .. External Subroutines ..                                         SLAP........820500
      EXTERNAL DCOPY, DRLCAL, DSCAL, DXLCAL                              SLAP........820600
C     .. Intrinsic Functions ..                                          SLAP........820700
      INTRINSIC ABS, MAX, SQRT                                           SLAP........820800
C     .. Common blocks ..                                                SLAP........820900
      COMMON /DSLBLK/ SOLN                                               SLAP........821000
C     .. Save statement ..                                               SLAP........821100
      SAVE SOLNRM                                                        SLAP........821200
C***FIRST EXECUTABLE STATEMENT  ISDGMR                                   SLAP........821300
      ISDGMR = 0                                                         SLAP........821400
      IF ( ITOL.EQ.0 ) THEN                                              SLAP........821500
C                                                                        SLAP........821600
C       Use input from DPIGMR to determine if stop conditions are met.   SLAP........821700
C                                                                        SLAP........821800
         ERR = RNRM/BNRM                                                 SLAP........821900
      ENDIF                                                              SLAP........822000
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN                          SLAP........822100
C                                                                        SLAP........822200
C       Use DRLCAL to calculate the scaled residual vector.              SLAP........822300
C       Store answer in R.                                               SLAP........822400
C                                                                        SLAP........822500
         IF ( LGMR.NE.0 ) CALL DRLCAL(N, KMP, LGMR, MAXL, V, Q, R,       SLAP........822600
     $                                SNORMW, PROD, R0NRM)               SLAP........822700
         IF ( ITOL.LE.2 ) THEN                                           SLAP........822800
C         err = ||Residual||/||RightHandSide||(2-Norms).                 SLAP........822900
            ERR = DNRM2(N, R, 1)/BNRM                                    SLAP........823000
C                                                                        SLAP........823100
C         Unscale R by R0NRM*PROD when KMP < MAXL.                       SLAP........823200
C                                                                        SLAP........823300
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN                  SLAP........823400
               TEM = 1.0D0/(R0NRM*PROD)                                  SLAP........823500
               CALL DSCAL(N, TEM, R, 1)                                  SLAP........823600
            ENDIF                                                        SLAP........823700
         ELSEIF ( ITOL.EQ.3 ) THEN                                       SLAP........823800
C         err = Max |(Minv*Residual)(i)/x(i)|                            SLAP........823900
C         When JPRE .lt. 0, R already contains Minv*Residual.            SLAP........824000
            IF ( JPRE.GT.0 ) THEN                                        SLAP........824100
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK,       SLAP........824200
     $              IWORK)                                               SLAP........824300
               NMSL = NMSL + 1                                           SLAP........824400
            ENDIF                                                        SLAP........824500
C                                                                        SLAP........824600
C         Unscale R by R0NRM*PROD when KMP < MAXL.                       SLAP........824700
C                                                                        SLAP........824800
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN                  SLAP........824900
               TEM = 1.0D0/(R0NRM*PROD)                                  SLAP........825000
               CALL DSCAL(N, TEM, R, 1)                                  SLAP........825100
            ENDIF                                                        SLAP........825200
C                                                                        SLAP........825300
            FUZZ = D1MACH(1)                                             SLAP........825400
            IELMAX = 1                                                   SLAP........825500
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ)                      SLAP........825600
            DO 25 I = 2, N                                               SLAP........825700
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ)                      SLAP........825800
               IF( RAT.GT.RATMAX ) THEN                                  SLAP........825900
                  IELMAX = I                                             SLAP........826000
                  RATMAX = RAT                                           SLAP........826100
               ENDIF                                                     SLAP........826200
 25         CONTINUE                                                     SLAP........826300
            ERR = RATMAX                                                 SLAP........826400
            IF( RATMAX.LE.TOL ) ISDGMR = 1                               SLAP........826500
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX      SLAP........826600
            RETURN                                                       SLAP........826700
         ENDIF                                                           SLAP........826800
      ENDIF                                                              SLAP........826900
      IF ( ITOL.EQ.11 ) THEN                                             SLAP........827000
C                                                                        SLAP........827100
C       Use DXLCAL to calculate the approximate solution XL.             SLAP........827200
C                                                                        SLAP........827300
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN                       SLAP........827400
            CALL DXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM,    SLAP........827500
     $           DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK,        SLAP........827600
     $           NELT, IA, JA, A, ISYM)                                  SLAP........827700
         ELSEIF ( ITER.EQ.0 ) THEN                                       SLAP........827800
C         Copy X to XL to check if initial guess is good enough.         SLAP........827900
            CALL DCOPY(N, X, 1, XL, 1)                                   SLAP........828000
         ELSE                                                            SLAP........828100
C         Return since this is the first call to DPIGMR on a restart.    SLAP........828200
            RETURN                                                       SLAP........828300
         ENDIF                                                           SLAP........828400
C                                                                        SLAP........828500
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN                     SLAP........828600
C         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).            SLAP........828700
            IF ( ITER.EQ.0 ) SOLNRM = DNRM2(N, SOLN, 1)                  SLAP........828800
            DO 30 I = 1, N                                               SLAP........828900
               DZ(I) = XL(I) - SOLN(I)                                   SLAP........829000
 30         CONTINUE                                                     SLAP........829100
            ERR = DNRM2(N, DZ, 1)/SOLNRM                                 SLAP........829200
         ELSE                                                            SLAP........829300
            IF (ITER .EQ. 0) THEN                                        SLAP........829400
               SOLNRM = 0                                                SLAP........829500
               DO 40 I = 1,N                                             SLAP........829600
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2                   SLAP........829700
 40            CONTINUE                                                  SLAP........829800
               SOLNRM = SQRT(SOLNRM)                                     SLAP........829900
            ENDIF                                                        SLAP........830000
            DXNRM = 0                                                    SLAP........830100
            DO 50 I = 1,N                                                SLAP........830200
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2                SLAP........830300
 50         CONTINUE                                                     SLAP........830400
            DXNRM = SQRT(DXNRM)                                          SLAP........830500
C         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).   SLAP........830600
            ERR = DXNRM/SOLNRM                                           SLAP........830700
         ENDIF                                                           SLAP........830800
      ENDIF                                                              SLAP........830900
C                                                                        SLAP........831000
      IF( IUNIT.NE.0 ) THEN                                              SLAP........831100
         IF( ITER.EQ.0 ) THEN                                            SLAP........831200
C...........THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION        SLAP........831300
C              OF SLAP WITH SUTRA.  IT ORIGINALLY READ:                  SLAP........831400
C              WRITE(IUNIT,1000) N, ITOL, MAXL, KMP                      SLAP........831500
            WRITE(IUNIT,1000)                                            SLAP........831600
         ENDIF                                                           SLAP........831700
C........THE NEXT LINE OF CODE WAS MODIFIED DURING INTEGRATION           SLAP........831800
C           OF SLAP WITH SUTRA.  IT ORIGINALLY READ:                     SLAP........831900
C           WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR                       SLAP........832000
         WRITE(IUNIT,1010) ITER, ERR                                     SLAP........832100
      ENDIF                                                              SLAP........832200
      IF ( ERR.LE.TOL ) ISDGMR = 1                                       SLAP........832300
C                                                                        SLAP........832400
      RETURN                                                             SLAP........832500
C.....THE NEXT TWO LINES OF CODE REFLECT MODIFICATIONS MADE DURING       SLAP........832600
C        INTEGRATION FO SLAP WITH SUTRA.  THEY ORIGINALLY READ:          SLAP........832700
C        1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ',    SLAP........832800
C            $     'N, ITOL = ',I5, I5,                                  SLAP........832900
C            $     /' ITER','   Natural Err Est','   Error Estimate')    SLAP........833000
C        1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7)                            SLAP........833100
1000  FORMAT(1X,6X,'Solver Iteration','   Error Estimate')               SLAP........833200
1010  FORMAT(1X,6X,I16,1X,1PE16.7)                                       SLAP........833300
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5,                          SLAP........833400
     $     ' |R(IELMAX)/X(IELMAX)| = ',D12.5)                            SLAP........833500
C------------- LAST LINE OF ISDGMR FOLLOWS ----------------------------  SLAP........833600
      END                                                                SLAP........833700
*DECK ISDOMN                                                             SLAP........833800
      INTEGER FUNCTION ISDOMN (N, B, X, NELT, IA, JA, A, ISYM, MSOLVE,   SLAP........833900
     +   NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,   SLAP........834000
     +   EMAP, DZ, CSAV, RWORK, IWORK, AK, BNRM, SOLNRM)                 SLAP........834100
C***BEGIN PROLOGUE  ISDOMN                                               SLAP........834200
C***SUBSIDIARY                                                           SLAP........834300
C***PURPOSE  Preconditioned Orthomin Stop Test.                          SLAP........834400
C            This routine calculates the stop test for the Orthomin      SLAP........834500
C            iteration scheme.  It returns a non-zero if the error       SLAP........834600
C            estimate (the type of which is determined by ITOL) is       SLAP........834700
C            less than the user specified tolerance TOL.                 SLAP........834800
C***LIBRARY   SLATEC (SLAP)                                              SLAP........834900
C***CATEGORY  D2A4, D2B4                                                 SLAP........835000
C***TYPE      DOUBLE PRECISION (ISSOMN-S, ISDOMN-D)                      SLAP........835100
C***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM,       SLAP........835200
C             ORTHOMIN, SLAP, SPARSE, STOP TEST                          SLAP........835300
C***AUTHOR  Greenbaum, Anne, (Courant Institute)                         SLAP........835400
C           Seager, Mark K., (LLNL)                                      SLAP........835500
C             Lawrence Livermore National Laboratory                     SLAP........835600
C             PO BOX 808, L-60                                           SLAP........835700
C             Livermore, CA 94550 (510) 423-3141                         SLAP........835800
C             seager@llnl.gov                                            SLAP........835900
C***DESCRIPTION                                                          SLAP........836000
C                                                                        SLAP........836100
C *Usage:                                                                SLAP........836200
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX     SLAP........836300
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)                    SLAP........836400
C     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)         SLAP........836500
C     DOUBLE PRECISION P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE)      SLAP........836600
C     DOUBLE PRECISION DZ(N), CSAV(NSAVE), RWORK(USER DEFINED), AK       SLAP........836700
C     DOUBLE PRECISION BNRM, SOLNRM                                      SLAP........836800
C     EXTERNAL MSOLVE                                                    SLAP........836900
C                                                                        SLAP........837000
C     IF( ISDOMN(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, NSAVE,          SLAP........837100
C    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, AP,        SLAP........837200
C    $     EMAP, DZ, CSAV, RWORK, IWORK, AK, BNRM, SOLNRM)               SLAP........837300
C    $     .NE.0 ) THEN ITERATION CONVERGED                              SLAP........837400
C                                                                        SLAP........837500
C *Arguments:                                                            SLAP........837600
C N      :IN       Integer.                                              SLAP........837700
C         Order of the matrix.                                           SLAP........837800
C B      :IN       Double Precision B(N).                                SLAP........837900
C         Right-hand side vector.                                        SLAP........838000
C X      :IN       Double Precision X(N).                                SLAP........838100
C         On input X is your initial guess for solution vector.          SLAP........838200
C         On output X is the final approximate solution.                 SLAP........838300
C NELT   :IN       Integer.                                              SLAP........838400
C         Number of Non-Zeros stored in A.                               SLAP........838500
C IA     :IN       Integer IA(NELT).                                     SLAP........838600
C JA     :IN       Integer JA(NELT).                                     SLAP........838700
C A      :IN       Double Precision A(NELT).                             SLAP........838800
C         These arrays should hold the matrix A in either the SLAP       SLAP........838900
C         Triad format or the SLAP Column format.  See "Description"     SLAP........839000
C         in the DSDOMN or DSLUOM prologue.                              SLAP........839100
C ISYM   :IN       Integer.                                              SLAP........839200
C         Flag to indicate symmetric storage format.                     SLAP........839300
C         If ISYM=0, all non-zero entries of the matrix are stored.      SLAP........839400
C         If ISYM=1, the matrix is symmetric, and only the upper         SLAP........839500
C         or lower triangle of the matrix is stored.                     SLAP........839600
C MSOLVE :EXT      External.                                             SLAP........839700
C         Name of a routine which solves a linear system MZ = R for      SLAP........839800
C         Z given R with the preconditioning matrix M (M is supplied via SLAP........839900
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must  SLAP........840000
C         be declared external in the calling program.  The calling      SLAP........840100
C         sequence to MSOLVE is:                                         SLAP........840200
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)  SLAP........840300
C         Where N is the number of unknowns, R is the right-hand side    SLAP........840400
C         vector and Z is the solution upon return.  NELT, IA, JA, A and SLAP........840500
C         ISYM are defined as above.  RWORK is a double precision array  SLAP........840600
C         that can be used to pass necessary preconditioning information SLAP........840700
C         and/or workspace to MSOLVE.  IWORK is an integer work array    SLAP........840800
C         for the same purpose as RWORK.                                 SLAP........840900
C NSAVE  :IN       Integer.                                              SLAP........841000
C         Number of direction vectors to save and orthogonalize against. SLAP........841100
C ITOL   :IN       Integer.                                              SLAP........841200
C         Flag to indicate type of convergence criterion.                SLAP........841300
C         If ITOL=1, iteration stops when the 2-norm of the residual     SLAP........841400
C         divided by the 2-norm of the right-hand side is less than TOL. SLAP........841500
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the  SLAP........841600
C         residual divided by the 2-norm of M-inv times the right hand   SLAP........841700
C         side is less than TOL, where M-inv is the inverse of the       SLAP........841800
C         diagonal of A.                                                 SLAP........841900
C         ITOL=11 is often useful for checking and comparing different   SLAP........842000
C         routines.  For this case, the user must supply the "exact"     SLAP........842100
C         solution or a very accurate approximation (one with an error   SLAP........842200
C         much less than TOL) through a common block,                    SLAP........842300
C             COMMON /DSLBLK/ SOLN( )                                    SLAP........842400
C         If ITOL=11, iteration stops when the 2-norm of the difference  SLAP........842500
C         between the iterative approximation and the user-supplied      SLAP........842600
C         solution divided by the 2-norm of the user-supplied solution   SLAP........842700
C         is less than TOL.  Note that this requires the user to set up  SLAP........842800
C         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.     SLAP........842900
C         The routine with this declaration should be loaded before the  SLAP........843000
C         stop test so that the correct length is used by the loader.    SLAP........843100
C         This procedure is not standard Fortran and may not work        SLAP........843200
C         correctly on your system (although it has worked on every      SLAP........843300
C         system the authors have tried).  If ITOL is not 11 then this   SLAP........843400
C         common block is indeed standard Fortran.                       SLAP........843500
C TOL    :IN       Double Precision.                                     SLAP........843600
C         Convergence criterion, as described above.                     SLAP........843700
C ITMAX  :IN       Integer.                                              SLAP........843800
C         Maximum number of iterations.                                  SLAP........843900
C ITER   :IN       Integer.                                              SLAP........844000
C         Current iteration count.  (Must be zero on first call.)        SLAP........844100
C ERR    :OUT      Double Precision.                                     SLAP........844200
C         Error estimate of error in final approximate solution, as      SLAP........844300
C         defined by ITOL.                                               SLAP........844400
C IERR   :OUT      Integer.                                              SLAP........844500
C         Error flag.  IERR is set to 3 if ITOL is not one of the        SLAP........844600
C         acceptable values, see above.                                  SLAP........844700
C IUNIT  :IN       Integer.                                              SLAP........844800
C         Unit number on which to write the error at each iteration,     SLAP........844900
C         if this is desired for monitoring convergence.  If unit        SLAP........845000
C         number is 0, no writing will occur.                            SLAP........845100
C R      :IN       Double Precision R(N).                                SLAP........845200
C         The residual R = B-AX.                                         SLAP........845300
C Z      :WORK     Double Precision Z(N).                                SLAP........845400
C P      :IN       Double Precision P(N,0:NSAVE).                        SLAP........845500
C         Workspace used to hold the conjugate direction vector(s).      SLAP........845600
C AP     :IN       Double Precision AP(N,0:NSAVE).                       SLAP........845700
C         Workspace used to hold the matrix A times the P vector(s).     SLAP........845800
C EMAP   :IN       Double Precision EMAP(N,0:NSAVE).                     SLAP........845900
C         Workspace used to hold M-inv times the AP vector(s).           SLAP........846000
C DZ     :WORK     Double Precision DZ(N).                               SLAP........846100
C         Workspace.                                                     SLAP........846200
C CSAV   :DUMMY    Double Precision CSAV(NSAVE)                          SLAP........846300
C         Reserved for future use.                                       SLAP........846400
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).                 SLAP........846500
C         Double Precision array that can be used for workspace in       SLAP........846600
C         MSOLVE.                                                        SLAP........846700
C IWORK  :WORK     Integer IWORK(USER DEFINED).                          SLAP........846800
C         Integer array that can be used for workspace in MSOLVE.        SLAP........846900
C AK     :IN       Double Precision.                                     SLAP........847000
C         Current iterate Orthomin iteration parameter.                  SLAP........847100
C BNRM   :OUT      Double Precision.                                     SLAP........847200
C         Current solution B-norm, if ITOL = 1 or 2.                     SLAP........847300
C SOLNRM :OUT      Double Precision.                                     SLAP........847400
C         True solution norm, if ITOL = 11.                              SLAP........847500
C                                                                        SLAP........847600
C *Function Return Values:                                               SLAP........847700
C       0 : Error estimate (determined by ITOL) is *NOT* less than the   SLAP........847800
C           specified tolerance, TOL.  The iteration must continue.      SLAP........847900
C       1 : Error estimate (determined by ITOL) is less than the         SLAP........848000
C           specified tolerance, TOL.  The iteration can be considered   SLAP........848100
C           complete.                                                    SLAP........848200
C                                                                        SLAP........848300
C *Cautions:                                                             SLAP........848400
C     This routine will attempt to write to the Fortran logical output   SLAP........848500
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that   SLAP........848600
C     this logical unit is attached to a file or terminal before calling SLAP........848700
C     this routine with a non-zero value for IUNIT.  This routine does   SLAP........848800
C     not check for the validity of a non-zero IUNIT unit number.        SLAP........848900
C                                                                        SLAP........849000
C***SEE ALSO  DOMN, DSDOMN, DSLUOM                                       SLAP........849100
C***ROUTINES CALLED  D1MACH, DNRM2                                       SLAP........849200
C***COMMON BLOCKS    DSLBLK                                              SLAP........849300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........849400
C   890404  DATE WRITTEN                                                 SLAP........849500
C   890404  Previous REVISION DATE                                       SLAP........849600
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........849700
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........849800
C           standard.  (FNF)                                             SLAP........849900
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........850000
C   891003  Removed C***REFER TO line, per MKS.                          SLAP........850100
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........850200
C   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)             SLAP........850300
C   910506  Made subsidiary to DOMN.  (FNF)                              SLAP........850400
C   920407  COMMON BLOCK renamed DSLBLK.  (WRB)                          SLAP........850500
C   920511  Added complete declaration section.  (WRB)                   SLAP........850600
C   920930  Corrected to not print AK when ITER=0.  (FNF)                SLAP........850700
C   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in          SLAP........850800
C           output format.  (FNF)                                        SLAP........850900
C   921113  Corrected C***CATEGORY line.  (FNF)                          SLAP........851000
C***END PROLOGUE  ISDOMN                                                 SLAP........851100
C     .. Scalar Arguments ..                                             SLAP........851200
      DOUBLE PRECISION AK, BNRM, ERR, SOLNRM, TOL                        SLAP........851300
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT, NSAVE       SLAP........851400
C     .. Array Arguments ..                                              SLAP........851500
      DOUBLE PRECISION A(NELT), AP(N,0:NSAVE), B(N), CSAV(NSAVE),        SLAP........851600
     +                 DZ(N), EMAP(N,0:NSAVE), P(N,0:NSAVE), R(N),       SLAP........851700
     +                 RWORK(*), X(N), Z(N)                              SLAP........851800
      INTEGER IA(NELT), IWORK(*), JA(NELT)                               SLAP........851900
C     .. Subroutine Arguments ..                                         SLAP........852000
      EXTERNAL MSOLVE                                                    SLAP........852100
C     .. Arrays in Common ..                                             SLAP........852200
      DOUBLE PRECISION SOLN(1)                                           SLAP........852300
C     .. Local Scalars ..                                                SLAP........852400
      INTEGER I                                                          SLAP........852500
C     .. External Functions ..                                           SLAP........852600
      DOUBLE PRECISION D1MACH, DNRM2                                     SLAP........852700
      EXTERNAL D1MACH, DNRM2                                             SLAP........852800
C     .. Common blocks ..                                                SLAP........852900
      COMMON /DSLBLK/ SOLN                                               SLAP........853000
C***FIRST EXECUTABLE STATEMENT  ISDOMN                                   SLAP........853100
      ISDOMN = 0                                                         SLAP........853200
C                                                                        SLAP........853300
      IF( ITOL.EQ.1 ) THEN                                               SLAP........853400
C         err = ||Residual||/||RightHandSide|| (2-Norms).                SLAP........853500
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)                           SLAP........853600
         ERR = DNRM2(N, R, 1)/BNRM                                       SLAP........853700
      ELSE IF( ITOL.EQ.2 ) THEN                                          SLAP........853800
C                  -1              -1                                    SLAP........853900
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).          SLAP........854000
         IF(ITER .EQ. 0) THEN                                            SLAP........854100
            CALL MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)   SLAP........854200
            BNRM = DNRM2(N, DZ, 1)                                       SLAP........854300
         ENDIF                                                           SLAP........854400
         ERR = DNRM2(N, Z, 1)/BNRM                                       SLAP........854500
      ELSE IF( ITOL.EQ.11 ) THEN                                         SLAP........854600
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).           SLAP........854700
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)                      SLAP........854800
         DO 10 I = 1, N                                                  SLAP........854900
            DZ(I) = X(I) - SOLN(I)                                       SLAP........855000
 10      CONTINUE                                                        SLAP........855100
         ERR = DNRM2(N, DZ, 1)/SOLNRM                                    SLAP........855200
      ELSE                                                               SLAP........855300
C                                                                        SLAP........855400
C         If we get here ITOL is not one of the acceptable values.       SLAP........855500
         ERR = D1MACH(2)                                                 SLAP........855600
         IERR = 3                                                        SLAP........855700
      ENDIF                                                              SLAP........855800
C                                                                        SLAP........855900
C.....THE BLOCK IF STATEMENT THAT IMMEDIATELY FOLLOWS THIS COMMMENT      SLAP........856000
C        WAS MODIFIED DURING INTEGRATION OF SLAP WITH SUTRA.  IT         SLAP........856100
C        ORIGINALLY READ:                                                SLAP........856200
C        IF(IUNIT .NE. 0) THEN                                           SLAP........856300
C           IF( ITER.EQ.0 ) THEN                                         SLAP........856400
C              WRITE(IUNIT,1000) NSAVE, N, ITOL                          SLAP........856500
C              WRITE(IUNIT,1010) ITER, ERR                               SLAP........856600
C           ELSE                                                         SLAP........856700
C              WRITE(IUNIT,1010) ITER, ERR, AK                           SLAP........856800
C           ENDIF                                                        SLAP........856900
C        ENDIF                                                           SLAP........857000
      IF(IUNIT .NE. 0) THEN                                              SLAP........857100
         IF( ITER.EQ.0 ) WRITE(IUNIT,1000)                               SLAP........857200
         WRITE(IUNIT,1010) ITER, ERR                                     SLAP........857300
      ENDIF                                                              SLAP........857400
      IF(ERR .LE. TOL) ISDOMN = 1                                        SLAP........857500
C                                                                        SLAP........857600
      RETURN                                                             SLAP........857700
C.....THE NEXT TWO LINES OF CODE REFLECT MODIFICATIONS MADE DURING       SLAP........857800
C        INTEGRATION OF SLAP WITH SUTRA.  THEY ORIGINALLY READ:          SLAP........857900
C        1000 FORMAT(' Preconditioned Orthomin(',I3,') for ',            SLAP........858000
C            $     'N, ITOL = ',I5, I5,                                  SLAP........858100
C            $     /' ITER','   Error Estimate','            Alpha')     SLAP........858200
C        1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7)                            SLAP........858300
1000  FORMAT(1X,6X,'Solver Iteration','   Error Estimate')               SLAP........858400
1010  FORMAT(1X,6X,I16,1X,1PE16.7)                                       SLAP........858500
C------------- LAST LINE OF ISDOMN FOLLOWS ----------------------------  SLAP........858600
      END                                                                SLAP........858700
*DECK J4SAVE                                                             SLAP........858800
      FUNCTION J4SAVE (IWHICH, IVALUE, ISET)                             SLAP........858900
C***BEGIN PROLOGUE  J4SAVE                                               SLAP........859000
C***SUBSIDIARY                                                           SLAP........859100
C***PURPOSE  Save or recall global variables needed by error             SLAP........859200
C            handling routines.                                          SLAP........859300
C***LIBRARY   SLATEC (XERROR)                                            SLAP........859400
C***TYPE      INTEGER (J4SAVE-I)                                         SLAP........859500
C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR         SLAP........859600
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........859700
C***DESCRIPTION                                                          SLAP........859800
C                                                                        SLAP........859900
C     Abstract                                                           SLAP........860000
C        J4SAVE saves and recalls several global variables needed        SLAP........860100
C        by the library error handling routines.                         SLAP........860200
C                                                                        SLAP........860300
C     Description of Parameters                                          SLAP........860400
C      --Input--                                                         SLAP........860500
C        IWHICH - Index of item desired.                                 SLAP........860600
C                = 1 Refers to current error number.                     SLAP........860700
C                = 2 Refers to current error control flag.               SLAP........860800
C                = 3 Refers to current unit number to which error        SLAP........860900
C                    messages are to be sent.  (0 means use standard.)   SLAP........861000
C                = 4 Refers to the maximum number of times any           SLAP........861100
C                     message is to be printed (as set by XERMAX).       SLAP........861200
C                = 5 Refers to the total number of units to which        SLAP........861300
C                     each error message is to be written.               SLAP........861400
C                = 6 Refers to the 2nd unit for error messages           SLAP........861500
C                = 7 Refers to the 3rd unit for error messages           SLAP........861600
C                = 8 Refers to the 4th unit for error messages           SLAP........861700
C                = 9 Refers to the 5th unit for error messages           SLAP........861800
C        IVALUE - The value to be set for the IWHICH-th parameter,       SLAP........861900
C                 if ISET is .TRUE. .                                    SLAP........862000
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE        SLAP........862100
C                 given the value, IVALUE.  If ISET=.FALSE., the         SLAP........862200
C                 IWHICH-th parameter will be unchanged, and IVALUE      SLAP........862300
C                 is a dummy parameter.                                  SLAP........862400
C      --Output--                                                        SLAP........862500
C        The (old) value of the IWHICH-th parameter will be returned     SLAP........862600
C        in the function value, J4SAVE.                                  SLAP........862700
C                                                                        SLAP........862800
C***SEE ALSO  XERMSG                                                     SLAP........862900
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........863000
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........863100
C                 Laboratories, 1982.                                    SLAP........863200
C***ROUTINES CALLED  (NONE)                                              SLAP........863300
C***REVISION HISTORY  (YYMMDD)                                           SLAP........863400
C   790801  DATE WRITTEN                                                 SLAP........863500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........863600
C   900205  Minor modifications to prologue.  (WRB)                      SLAP........863700
C   900402  Added TYPE section.  (WRB)                                   SLAP........863800
C   910411  Added KEYWORDS section.  (WRB)                               SLAP........863900
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........864000
C***END PROLOGUE  J4SAVE                                                 SLAP........864100
      LOGICAL ISET                                                       SLAP........864200
      INTEGER IPARAM(9)                                                  SLAP........864300
      SAVE IPARAM                                                        SLAP........864400
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/             SLAP........864500
      DATA IPARAM(5)/1/                                                  SLAP........864600
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/              SLAP........864700
C***FIRST EXECUTABLE STATEMENT  J4SAVE                                   SLAP........864800
      J4SAVE = IPARAM(IWHICH)                                            SLAP........864900
      IF (ISET) IPARAM(IWHICH) = IVALUE                                  SLAP........865000
      RETURN                                                             SLAP........865100
      END                                                                SLAP........865200
*DECK QS2I1D                                                             SLAP........865300
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)                            SLAP........865400
C***BEGIN PROLOGUE  QS2I1D                                               SLAP........865500
C***SUBSIDIARY                                                           SLAP........865600
C***PURPOSE  Sort an integer array, moving an integer and DP array.      SLAP........865700
C            This routine sorts the integer array IA and makes the same  SLAP........865800
C            interchanges in the integer array JA and the double pre-    SLAP........865900
C            cision array A.  The array IA may be sorted in increasing   SLAP........866000
C            order or decreasing order.  A slightly modified QUICKSORT   SLAP........866100
C            algorithm is used.                                          SLAP........866200
C***LIBRARY   SLATEC (SLAP)                                              SLAP........866300
C***CATEGORY  N6A2A                                                      SLAP........866400
C***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)                      SLAP........866500
C***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING                   SLAP........866600
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........866700
C           Kahaner, D. K., (NBS)                                        SLAP........866800
C           Seager, M. K., (LLNL) seager@llnl.gov                        SLAP........866900
C           Wisniewski, J. A., (SNLA)                                    SLAP........867000
C***DESCRIPTION                                                          SLAP........867100
C     Written by Rondall E Jones                                         SLAP........867200
C     Modified by John A. Wisniewski to use the Singleton QUICKSORT      SLAP........867300
C     algorithm. date 18 November 1976.                                  SLAP........867400
C                                                                        SLAP........867500
C     Further modified by David K. Kahaner                               SLAP........867600
C     National Bureau of Standards                                       SLAP........867700
C     August, 1981                                                       SLAP........867800
C                                                                        SLAP........867900
C     Even further modification made to bring the code up to the         SLAP........868000
C     Fortran 77 level and make it more readable and to carry            SLAP........868100
C     along one integer array and one double precision array during      SLAP........868200
C     the sort by                                                        SLAP........868300
C     Mark K. Seager                                                     SLAP........868400
C     Lawrence Livermore National Laboratory                             SLAP........868500
C     November, 1987                                                     SLAP........868600
C     This routine was adapted from the ISORT routine.                   SLAP........868700
C                                                                        SLAP........868800
C     ABSTRACT                                                           SLAP........868900
C         This routine sorts an integer array IA and makes the same      SLAP........869000
C         interchanges in the integer array JA and the double precision  SLAP........869100
C         array A.                                                       SLAP........869200
C         The array IA may be sorted in increasing order or decreasing   SLAP........869300
C         order.  A slightly modified quicksort algorithm is used.       SLAP........869400
C                                                                        SLAP........869500
C     DESCRIPTION OF PARAMETERS                                          SLAP........869600
C        IA - Integer array of values to be sorted.                      SLAP........869700
C        JA - Integer array to be carried along.                         SLAP........869800
C         A - Double Precision array to be carried along.                SLAP........869900
C         N - Number of values in integer array IA to be sorted.         SLAP........870000
C     KFLAG - Control parameter                                          SLAP........870100
C           = 1 means sort IA in INCREASING order.                       SLAP........870200
C           =-1 means sort IA in DECREASING order.                       SLAP........870300
C                                                                        SLAP........870400
C***SEE ALSO  DS2Y                                                       SLAP........870500
C***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm   SLAP........870600
C                 for Sorting With Minimal Storage, Communications ACM   SLAP........870700
C                 12:3 (1969), pp.185-7.                                 SLAP........870800
C***ROUTINES CALLED  XERMSG                                              SLAP........870900
C***REVISION HISTORY  (YYMMDD)                                           SLAP........871000
C   761118  DATE WRITTEN                                                 SLAP........871100
C   890125  Previous REVISION DATE                                       SLAP........871200
C   890915  Made changes requested at July 1989 CML Meeting.  (MKS)      SLAP........871300
C   890922  Numerous changes to prologue to make closer to SLATEC        SLAP........871400
C           standard.  (FNF)                                             SLAP........871500
C   890929  Numerous changes to reduce SP/DP differences.  (FNF)         SLAP........871600
C   900805  Changed XERROR calls to calls to XERMSG.  (RWC)              SLAP........871700
C   910411  Prologue converted to Version 4.0 format.  (BAB)             SLAP........871800
C   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)      SLAP........871900
C   920511  Added complete declaration section.  (WRB)                   SLAP........872000
C   920929  Corrected format of reference.  (FNF)                        SLAP........872100
C   921012  Corrected all f.p. constants to double precision.  (FNF)     SLAP........872200
C***END PROLOGUE  QS2I1D                                                 SLAP........872300
CVD$R NOVECTOR                                                           SLAP........872400
CVD$R NOCONCUR                                                           SLAP........872500
C     .. Scalar Arguments ..                                             SLAP........872600
      INTEGER KFLAG, N                                                   SLAP........872700
C     .. Array Arguments ..                                              SLAP........872800
      DOUBLE PRECISION A(N)                                              SLAP........872900
      INTEGER IA(N), JA(N)                                               SLAP........873000
C     .. Local Scalars ..                                                SLAP........873100
      DOUBLE PRECISION R, TA, TTA                                        SLAP........873200
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN                SLAP........873300
C     .. Local Arrays ..                                                 SLAP........873400
      INTEGER IL(21), IU(21)                                             SLAP........873500
C     .. External Subroutines ..                                         SLAP........873600
      EXTERNAL XERMSG                                                    SLAP........873700
C     .. Intrinsic Functions ..                                          SLAP........873800
      INTRINSIC ABS, INT                                                 SLAP........873900
C***FIRST EXECUTABLE STATEMENT  QS2I1D                                   SLAP........874000
      NN = N                                                             SLAP........874100
      IF (NN.LT.1) THEN                                                  SLAP........874200
         CALL XERMSG ('SLATEC', 'QS2I1D',                                SLAP........874300
     $      'The number of values to be sorted was not positive.', 1, 1) SLAP........874400
         RETURN                                                          SLAP........874500
      ENDIF                                                              SLAP........874600
      IF( N.EQ.1 ) RETURN                                                SLAP........874700
      KK = ABS(KFLAG)                                                    SLAP........874800
      IF ( KK.NE.1 ) THEN                                                SLAP........874900
         CALL XERMSG ('SLATEC', 'QS2I1D',                                SLAP........875000
     $      'The sort control parameter, K, was not 1 or -1.', 2, 1)     SLAP........875100
         RETURN                                                          SLAP........875200
      ENDIF                                                              SLAP........875300
C                                                                        SLAP........875400
C     Alter array IA to get decreasing order if needed.                  SLAP........875500
C                                                                        SLAP........875600
      IF( KFLAG.LT.1 ) THEN                                              SLAP........875700
         DO 20 I=1,NN                                                    SLAP........875800
            IA(I) = -IA(I)                                               SLAP........875900
 20      CONTINUE                                                        SLAP........876000
      ENDIF                                                              SLAP........876100
C                                                                        SLAP........876200
C     Sort IA and carry JA and A along.                                  SLAP........876300
C     And now...Just a little black magic...                             SLAP........876400
      M = 1                                                              SLAP........876500
      I = 1                                                              SLAP........876600
      J = NN                                                             SLAP........876700
      R = .375D0                                                         SLAP........876800
 210  IF( R.LE.0.5898437D0 ) THEN                                        SLAP........876900
         R = R + 3.90625D-2                                              SLAP........877000
      ELSE                                                               SLAP........877100
         R = R-.21875D0                                                  SLAP........877200
      ENDIF                                                              SLAP........877300
 225  K = I                                                              SLAP........877400
C                                                                        SLAP........877500
C     Select a central element of the array and save it in location      SLAP........877600
C     it, jt, at.                                                        SLAP........877700
C                                                                        SLAP........877800
      IJ = I + INT ((J-I)*R)                                             SLAP........877900
      IT = IA(IJ)                                                        SLAP........878000
      JT = JA(IJ)                                                        SLAP........878100
      TA = A(IJ)                                                         SLAP........878200
C                                                                        SLAP........878300
C     If first element of array is greater than it, interchange with it. SLAP........878400
C                                                                        SLAP........878500
      IF( IA(I).GT.IT ) THEN                                             SLAP........878600
         IA(IJ) = IA(I)                                                  SLAP........878700
         IA(I)  = IT                                                     SLAP........878800
         IT     = IA(IJ)                                                 SLAP........878900
         JA(IJ) = JA(I)                                                  SLAP........879000
         JA(I)  = JT                                                     SLAP........879100
         JT     = JA(IJ)                                                 SLAP........879200
         A(IJ)  = A(I)                                                   SLAP........879300
         A(I)   = TA                                                     SLAP........879400
         TA     = A(IJ)                                                  SLAP........879500
      ENDIF                                                              SLAP........879600
      L=J                                                                SLAP........879700
C                                                                        SLAP........879800
C     If last element of array is less than it, swap with it.            SLAP........879900
C                                                                        SLAP........880000
      IF( IA(J).LT.IT ) THEN                                             SLAP........880100
         IA(IJ) = IA(J)                                                  SLAP........880200
         IA(J)  = IT                                                     SLAP........880300
         IT     = IA(IJ)                                                 SLAP........880400
         JA(IJ) = JA(J)                                                  SLAP........880500
         JA(J)  = JT                                                     SLAP........880600
         JT     = JA(IJ)                                                 SLAP........880700
         A(IJ)  = A(J)                                                   SLAP........880800
         A(J)   = TA                                                     SLAP........880900
         TA     = A(IJ)                                                  SLAP........881000
C                                                                        SLAP........881100
C     If first element of array is greater than it, swap with it.        SLAP........881200
C                                                                        SLAP........881300
         IF ( IA(I).GT.IT ) THEN                                         SLAP........881400
            IA(IJ) = IA(I)                                               SLAP........881500
            IA(I)  = IT                                                  SLAP........881600
            IT     = IA(IJ)                                              SLAP........881700
            JA(IJ) = JA(I)                                               SLAP........881800
            JA(I)  = JT                                                  SLAP........881900
            JT     = JA(IJ)                                              SLAP........882000
            A(IJ)  = A(I)                                                SLAP........882100
            A(I)   = TA                                                  SLAP........882200
            TA     = A(IJ)                                               SLAP........882300
         ENDIF                                                           SLAP........882400
      ENDIF                                                              SLAP........882500
C                                                                        SLAP........882600
C     Find an element in the second half of the array which is           SLAP........882700
C     smaller than it.                                                   SLAP........882800
C                                                                        SLAP........882900
  240 L=L-1                                                              SLAP........883000
      IF( IA(L).GT.IT ) GO TO 240                                        SLAP........883100
C                                                                        SLAP........883200
C     Find an element in the first half of the array which is            SLAP........883300
C     greater than it.                                                   SLAP........883400
C                                                                        SLAP........883500
  245 K=K+1                                                              SLAP........883600
      IF( IA(K).LT.IT ) GO TO 245                                        SLAP........883700
C                                                                        SLAP........883800
C     Interchange these elements.                                        SLAP........883900
C                                                                        SLAP........884000
      IF( K.LE.L ) THEN                                                  SLAP........884100
         IIT   = IA(L)                                                   SLAP........884200
         IA(L) = IA(K)                                                   SLAP........884300
         IA(K) = IIT                                                     SLAP........884400
         JJT   = JA(L)                                                   SLAP........884500
         JA(L) = JA(K)                                                   SLAP........884600
         JA(K) = JJT                                                     SLAP........884700
         TTA   = A(L)                                                    SLAP........884800
         A(L)  = A(K)                                                    SLAP........884900
         A(K)  = TTA                                                     SLAP........885000
         GOTO 240                                                        SLAP........885100
      ENDIF                                                              SLAP........885200
C                                                                        SLAP........885300
C     Save upper and lower subscripts of the array yet to be sorted.     SLAP........885400
C                                                                        SLAP........885500
      IF( L-I.GT.J-K ) THEN                                              SLAP........885600
         IL(M) = I                                                       SLAP........885700
         IU(M) = L                                                       SLAP........885800
         I = K                                                           SLAP........885900
         M = M+1                                                         SLAP........886000
      ELSE                                                               SLAP........886100
         IL(M) = K                                                       SLAP........886200
         IU(M) = J                                                       SLAP........886300
         J = L                                                           SLAP........886400
         M = M+1                                                         SLAP........886500
      ENDIF                                                              SLAP........886600
      GO TO 260                                                          SLAP........886700
C                                                                        SLAP........886800
C     Begin again on another portion of the unsorted array.              SLAP........886900
C                                                                        SLAP........887000
  255 M = M-1                                                            SLAP........887100
      IF( M.EQ.0 ) GO TO 300                                             SLAP........887200
      I = IL(M)                                                          SLAP........887300
      J = IU(M)                                                          SLAP........887400
  260 IF( J-I.GE.1 ) GO TO 225                                           SLAP........887500
      IF( I.EQ.J ) GO TO 255                                             SLAP........887600
      IF( I.EQ.1 ) GO TO 210                                             SLAP........887700
      I = I-1                                                            SLAP........887800
  265 I = I+1                                                            SLAP........887900
      IF( I.EQ.J ) GO TO 255                                             SLAP........888000
      IT = IA(I+1)                                                       SLAP........888100
      JT = JA(I+1)                                                       SLAP........888200
      TA =  A(I+1)                                                       SLAP........888300
      IF( IA(I).LE.IT ) GO TO 265                                        SLAP........888400
      K=I                                                                SLAP........888500
  270 IA(K+1) = IA(K)                                                    SLAP........888600
      JA(K+1) = JA(K)                                                    SLAP........888700
      A(K+1)  =  A(K)                                                    SLAP........888800
      K = K-1                                                            SLAP........888900
      IF( IT.LT.IA(K) ) GO TO 270                                        SLAP........889000
      IA(K+1) = IT                                                       SLAP........889100
      JA(K+1) = JT                                                       SLAP........889200
      A(K+1)  = TA                                                       SLAP........889300
      GO TO 265                                                          SLAP........889400
C                                                                        SLAP........889500
C     Clean up, if necessary.                                            SLAP........889600
C                                                                        SLAP........889700
  300 IF( KFLAG.LT.1 ) THEN                                              SLAP........889800
         DO 310 I=1,NN                                                   SLAP........889900
            IA(I) = -IA(I)                                               SLAP........890000
 310     CONTINUE                                                        SLAP........890100
      ENDIF                                                              SLAP........890200
      RETURN                                                             SLAP........890300
C------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------  SLAP........890400
      END                                                                SLAP........890500
*DECK XERCNT                                                             SLAP........890600
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)     SLAP........890700
C***BEGIN PROLOGUE  XERCNT                                               SLAP........890800
C***SUBSIDIARY                                                           SLAP........890900
C***PURPOSE  Allow user control over handling of errors.                 SLAP........891000
C***LIBRARY   SLATEC (XERROR)                                            SLAP........891100
C***CATEGORY  R3C                                                        SLAP........891200
C***TYPE      ALL (XERCNT-A)                                             SLAP........891300
C***KEYWORDS  ERROR, XERROR                                              SLAP........891400
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........891500
C***DESCRIPTION                                                          SLAP........891600
C                                                                        SLAP........891700
C     Abstract                                                           SLAP........891800
C        Allows user control over handling of individual errors.         SLAP........891900
C        Just after each message is recorded, but before it is           SLAP........892000
C        processed any further (i.e., before it is printed or            SLAP........892100
C        a decision to abort is made), a call is made to XERCNT.         SLAP........892200
C        If the user has provided his own version of XERCNT, he          SLAP........892300
C        can then override the value of KONTROL used in processing       SLAP........892400
C        this message by redefining its value.                           SLAP........892500
C        KONTRL may be set to any value from -2 to 2.                    SLAP........892600
C        The meanings for KONTRL are the same as in XSETF, except        SLAP........892700
C        that the value of KONTRL changes only for this message.         SLAP........892800
C        If KONTRL is set to a value outside the range from -2 to 2,     SLAP........892900
C        it will be moved back into that range.                          SLAP........893000
C                                                                        SLAP........893100
C     Description of Parameters                                          SLAP........893200
C                                                                        SLAP........893300
C      --Input--                                                         SLAP........893400
C        LIBRAR - the library that the routine is in.                    SLAP........893500
C        SUBROU - the subroutine that XERMSG is being called from        SLAP........893600
C        MESSG  - the first 20 characters of the error message.          SLAP........893700
C        NERR   - same as in the call to XERMSG.                         SLAP........893800
C        LEVEL  - same as in the call to XERMSG.                         SLAP........893900
C        KONTRL - the current value of the control flag as set           SLAP........894000
C                 by a call to XSETF.                                    SLAP........894100
C                                                                        SLAP........894200
C      --Output--                                                        SLAP........894300
C        KONTRL - the new value of KONTRL.  If KONTRL is not             SLAP........894400
C                 defined, it will remain at its original value.         SLAP........894500
C                 This changed value of control affects only             SLAP........894600
C                 the current occurrence of the current message.         SLAP........894700
C                                                                        SLAP........894800
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........894900
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........895000
C                 Laboratories, 1982.                                    SLAP........895100
C***ROUTINES CALLED  (NONE)                                              SLAP........895200
C***REVISION HISTORY  (YYMMDD)                                           SLAP........895300
C   790801  DATE WRITTEN                                                 SLAP........895400
C   861211  REVISION DATE from Version 3.2                               SLAP........895500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........895600
C   900206  Routine changed from user-callable to subsidiary.  (WRB)     SLAP........895700
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE   SLAP........895800
C           names, changed routine name from XERCTL to XERCNT.  (RWC)    SLAP........895900
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........896000
C***END PROLOGUE  XERCNT                                                 SLAP........896100
      CHARACTER*(*) LIBRAR, SUBROU, MESSG                                SLAP........896200
C***FIRST EXECUTABLE STATEMENT  XERCNT                                   SLAP........896300
C.....THE NEXT TWO NON-COMMENT LINES (ONE CONTINUED LINE) WERE ADDED     SLAP........896400
C        DURING INTEGRATION OF SLAP WITH SUTRA.  THEY ALLOW THE SOLUTION SLAP........896500
C        PROCESS TO CONTINUE IF THE IC FACTORIZATION FOR THE CG SOLVER   SLAP........896600
C        IS FUDGED.                                                      SLAP........896700
      IF ((LIBRAR.EQ.'SLATEC').AND.(SUBROU.EQ.'DSICCG').AND.             SLAP........896800
     1    (MESSG(1:20).EQ.'IC factorization bro')) KONTRL = 0            SLAP........896900
      RETURN                                                             SLAP........897000
      END                                                                SLAP........897100
*DECK XERHLT                                                             SLAP........897200
      SUBROUTINE XERHLT (MESSG)                                          SLAP........897300
C***BEGIN PROLOGUE  XERHLT                                               SLAP........897400
C***SUBSIDIARY                                                           SLAP........897500
C***PURPOSE  Abort program execution and print error message.            SLAP........897600
C***LIBRARY   SLATEC (XERROR)                                            SLAP........897700
C***CATEGORY  R3C                                                        SLAP........897800
C***TYPE      ALL (XERHLT-A)                                             SLAP........897900
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR                     SLAP........898000
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........898100
C***DESCRIPTION                                                          SLAP........898200
C                                                                        SLAP........898300
C     Abstract                                                           SLAP........898400
C        ***Note*** machine dependent routine                            SLAP........898500
C        XERHLT aborts the execution of the program.                     SLAP........898600
C        The error message causing the abort is given in the calling     SLAP........898700
C        sequence, in case one needs it for printing on a dayfile,       SLAP........898800
C        for example.                                                    SLAP........898900
C                                                                        SLAP........899000
C     Description of Parameters                                          SLAP........899100
C        MESSG is as in XERMSG.                                          SLAP........899200
C                                                                        SLAP........899300
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........899400
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........899500
C                 Laboratories, 1982.                                    SLAP........899600
C***ROUTINES CALLED  (NONE)                                              SLAP........899700
C***REVISION HISTORY  (YYMMDD)                                           SLAP........899800
C   790801  DATE WRITTEN                                                 SLAP........899900
C   861211  REVISION DATE from Version 3.2                               SLAP........900000
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........900100
C   900206  Routine changed from user-callable to subsidiary.  (WRB)     SLAP........900200
C   900510  Changed calling sequence to delete length of character       SLAP........900300
C           and changed routine name from XERABT to XERHLT.  (RWC)       SLAP........900400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........900500
C***END PROLOGUE  XERHLT                                                 SLAP........900600
      CHARACTER*(*) MESSG                                                SLAP........900700
C***FIRST EXECUTABLE STATEMENT  XERHLT                                   SLAP........900800
      STOP                                                               SLAP........900900
      END                                                                SLAP........901000
*DECK XERMSG                                                             SLAP........901100
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)             SLAP........901200
C***BEGIN PROLOGUE  XERMSG                                               SLAP........901300
C***PURPOSE  Process error messages for SLATEC and other libraries.      SLAP........901400
C***LIBRARY   SLATEC (XERROR)                                            SLAP........901500
C***CATEGORY  R3C                                                        SLAP........901600
C***TYPE      ALL (XERMSG-A)                                             SLAP........901700
C***KEYWORDS  ERROR MESSAGE, XERROR                                      SLAP........901800
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)                                SLAP........901900
C***DESCRIPTION                                                          SLAP........902000
C                                                                        SLAP........902100
C   XERMSG processes a diagnostic message in a manner determined by the  SLAP........902200
C   value of LEVEL and the current value of the library error control    SLAP........902300
C   flag, KONTRL.  See subroutine XSETF for details.                     SLAP........902400
C                                                                        SLAP........902500
C    LIBRAR   A character constant (or character variable) with the name SLAP........902600
C             of the library.  This will be 'SLATEC' for the SLATEC      SLAP........902700
C             Common Math Library.  The error handling package is        SLAP........902800
C             general enough to be used by many libraries                SLAP........902900
C             simultaneously, so it is desirable for the routine that    SLAP........903000
C             detects and reports an error to identify the library name  SLAP........903100
C             as well as the routine name.                               SLAP........903200
C                                                                        SLAP........903300
C    SUBROU   A character constant (or character variable) with the name SLAP........903400
C             of the routine that detected the error.  Usually it is the SLAP........903500
C             name of the routine that is calling XERMSG.  There are     SLAP........903600
C             some instances where a user callable library routine calls SLAP........903700
C             lower level subsidiary routines where the error is         SLAP........903800
C             detected.  In such cases it may be more informative to     SLAP........903900
C             supply the name of the routine the user called rather than SLAP........904000
C             the name of the subsidiary routine that detected the       SLAP........904100
C             error.                                                     SLAP........904200
C                                                                        SLAP........904300
C    MESSG    A character constant (or character variable) with the text SLAP........904400
C             of the error or warning message.  In the example below,    SLAP........904500
C             the message is a character constant that contains a        SLAP........904600
C             generic message.                                           SLAP........904700
C                                                                        SLAP........904800
C                   CALL XERMSG ('SLATEC', 'MMPY',                       SLAP........904900
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION', SLAP........905000
C                  *3, 1)                                                SLAP........905100
C                                                                        SLAP........905200
C             It is possible (and is sometimes desirable) to generate a  SLAP........905300
C             specific message--e.g., one that contains actual numeric   SLAP........905400
C             values.  Specific numeric values can be converted into     SLAP........905500
C             character strings using formatted WRITE statements into    SLAP........905600
C             character variables.  This is called standard Fortran      SLAP........905700
C             internal file I/O and is exemplified in the first three    SLAP........905800
C             lines of the following example.  You can also catenate     SLAP........905900
C             substrings of characters to construct the error message.   SLAP........906000
C             Here is an example showing the use of both writing to      SLAP........906100
C             an internal file and catenating character strings.         SLAP........906200
C                                                                        SLAP........906300
C                   CHARACTER*5 CHARN, CHARL                             SLAP........906400
C                   WRITE (CHARN,10) N                                   SLAP........906500
C                   WRITE (CHARL,10) LDA                                 SLAP........906600
C                10 FORMAT(I5)                                           SLAP........906700
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//  SLAP........906800
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//   SLAP........906900
C                  *   CHARL, 3, 1)                                      SLAP........907000
C                                                                        SLAP........907100
C             There are two subtleties worth mentioning.  One is that    SLAP........907200
C             the // for character catenation is used to construct the   SLAP........907300
C             error message so that no single character constant is      SLAP........907400
C             continued to the next line.  This avoids confusion as to   SLAP........907500
C             whether there are trailing blanks at the end of the line.  SLAP........907600
C             The second is that by catenating the parts of the message  SLAP........907700
C             as an actual argument rather than encoding the entire      SLAP........907800
C             message into one large character variable, we avoid        SLAP........907900
C             having to know how long the message will be in order to    SLAP........908000
C             declare an adequate length for that large character        SLAP........908100
C             variable.  XERMSG calls XERPRN to print the message using  SLAP........908200
C             multiple lines if necessary.  If the message is very long, SLAP........908300
C             XERPRN will break it into pieces of 72 characters (as      SLAP........908400
C             requested by XERMSG) for printing on multiple lines.       SLAP........908500
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '   SLAP........908600
C             so that the total line length could be 76 characters.      SLAP........908700
C             Note also that XERPRN scans the error message backwards    SLAP........908800
C             to ignore trailing blanks.  Another feature is that        SLAP........908900
C             the substring '$$' is treated as a new line sentinel       SLAP........909000
C             by XERPRN.  If you want to construct a multiline           SLAP........909100
C             message without having to count out multiples of 72        SLAP........909200
C             characters, just use '$$' as a separator.  '$$'            SLAP........909300
C             obviously must occur within 72 characters of the           SLAP........909400
C             start of each line to have its intended effect since       SLAP........909500
C             XERPRN is asked to wrap around at 72 characters in         SLAP........909600
C             addition to looking for '$$'.                              SLAP........909700
C                                                                        SLAP........909800
C    NERR     An integer value that is chosen by the library routine's   SLAP........909900
C             author.  It must be in the range -99 to 999 (three         SLAP........910000
C             printable digits).  Each distinct error should have its    SLAP........910100
C             own error number.  These error numbers should be described SLAP........910200
C             in the machine readable documentation for the routine.     SLAP........910300
C             The error numbers need be unique only within each routine, SLAP........910400
C             so it is reasonable for each routine to start enumerating  SLAP........910500
C             errors from 1 and proceeding to the next integer.          SLAP........910600
C                                                                        SLAP........910700
C    LEVEL    An integer value in the range 0 to 2 that indicates the    SLAP........910800
C             level (severity) of the error.  Their meanings are         SLAP........910900
C                                                                        SLAP........911000
C            -1  A warning message.  This is used if it is not clear     SLAP........911100
C                that there really is an error, but the user's attention SLAP........911200
C                may be needed.  An attempt is made to only print this   SLAP........911300
C                message once.                                           SLAP........911400
C                                                                        SLAP........911500
C             0  A warning message.  This is used if it is not clear     SLAP........911600
C                that there really is an error, but the user's attention SLAP........911700
C                may be needed.                                          SLAP........911800
C                                                                        SLAP........911900
C             1  A recoverable error.  This is used even if the error is SLAP........912000
C                so serious that the routine cannot return any useful    SLAP........912100
C                answer.  If the user has told the error package to      SLAP........912200
C                return after recoverable errors, then XERMSG will       SLAP........912300
C                return to the Library routine which can then return to  SLAP........912400
C                the user's routine.  The user may also permit the error SLAP........912500
C                package to terminate the program upon encountering a    SLAP........912600
C                recoverable error.                                      SLAP........912700
C                                                                        SLAP........912800
C             2  A fatal error.  XERMSG will not return to its caller    SLAP........912900
C                after it receives a fatal error.  This level should     SLAP........913000
C                hardly ever be used; it is much better to allow the     SLAP........913100
C                user a chance to recover.  An example of one of the few SLAP........913200
C                cases in which it is permissible to declare a level 2   SLAP........913300
C                error is a reverse communication Library routine that   SLAP........913400
C                is likely to be called repeatedly until it integrates   SLAP........913500
C                across some interval.  If there is a serious error in   SLAP........913600
C                the input such that another step cannot be taken and    SLAP........913700
C                the Library routine is called again without the input   SLAP........913800
C                error having been corrected by the caller, the Library  SLAP........913900
C                routine will probably be called forever with improper   SLAP........914000
C                input.  In this case, it is reasonable to declare the   SLAP........914100
C                error to be fatal.                                      SLAP........914200
C                                                                        SLAP........914300
C    Each of the arguments to XERMSG is input; none will be modified by  SLAP........914400
C    XERMSG.  A routine may make multiple calls to XERMSG with warning   SLAP........914500
C    level messages; however, after a call to XERMSG with a recoverable  SLAP........914600
C    error, the routine should return to the user.  Do not try to call   SLAP........914700
C    XERMSG with a second recoverable error after the first recoverable  SLAP........914800
C    error because the error package saves the error number.  The user   SLAP........914900
C    can retrieve this error number by calling another entry point in    SLAP........915000
C    the error handling package and then clear the error number when     SLAP........915100
C    recovering from the error.  Calling XERMSG in succession causes the SLAP........915200
C    old error number to be overwritten by the latest error number.      SLAP........915300
C    This is considered harmless for error numbers associated with       SLAP........915400
C    warning messages but must not be done for error numbers of serious  SLAP........915500
C    errors.  After a call to XERMSG with a recoverable error, the user  SLAP........915600
C    must be given a chance to call NUMXER or XERCLR to retrieve or      SLAP........915700
C    clear the error number.                                             SLAP........915800
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........915900
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........916000
C                 Laboratories, 1982.                                    SLAP........916100
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE       SLAP........916200
C***REVISION HISTORY  (YYMMDD)                                           SLAP........916300
C   880101  DATE WRITTEN                                                 SLAP........916400
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.  SLAP........916500
C           THERE ARE TWO BASIC CHANGES.                                 SLAP........916600
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO      SLAP........916700
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES   SLAP........916800
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS     SLAP........916900
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE        SLAP........917000
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER    SLAP........917100
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY      SLAP........917200
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE         SLAP........917300
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.                 SLAP........917400
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE  SLAP........917500
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE   SLAP........917600
C               OF LOWER CASE.                                           SLAP........917700
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.      SLAP........917800
C           THE PRINCIPAL CHANGES ARE                                    SLAP........917900
C           1.  CLARIFY COMMENTS IN THE PROLOGUES                        SLAP........918000
C           2.  RENAME XRPRNT TO XERPRN                                  SLAP........918100
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES  SLAP........918200
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /        SLAP........918300
C               CHARACTER FOR NEW RECORDS.                               SLAP........918400
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO     SLAP........918500
C           CLEAN UP THE CODING.                                         SLAP........918600
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN  SLAP........918700
C           PREFIX.                                                      SLAP........918800
C   891013  REVISED TO CORRECT COMMENTS.                                 SLAP........918900
C   891214  Prologue converted to Version 4.0 format.  (WRB)             SLAP........919000
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but   SLAP........919100
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added       SLAP........919200
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and       SLAP........919300
C           XERCTL to XERCNT.  (RWC)                                     SLAP........919400
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........919500
C***END PROLOGUE  XERMSG                                                 SLAP........919600
      CHARACTER*(*) LIBRAR, SUBROU, MESSG                                SLAP........919700
      CHARACTER*8 XLIBR, XSUBR                                           SLAP........919800
      CHARACTER*72  TEMP                                                 SLAP........919900
      CHARACTER*20  LFIRST                                               SLAP........920000
C***FIRST EXECUTABLE STATEMENT  XERMSG                                   SLAP........920100
      LKNTRL = J4SAVE (2, 0, .FALSE.)                                    SLAP........920200
      MAXMES = J4SAVE (4, 0, .FALSE.)                                    SLAP........920300
C                                                                        SLAP........920400
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.               SLAP........920500
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE     SLAP........920600
C          SHOULD BE PRINTED.                                            SLAP........920700
C                                                                        SLAP........920800
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN     SLAP........920900
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,         SLAP........921000
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.                      SLAP........921100
C                                                                        SLAP........921200
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.     SLAP........921300
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN                               SLAP........921400
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //              SLAP........921500
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//               SLAP........921600
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)                         SLAP........921700
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)                    SLAP........921800
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')                     SLAP........921900
         RETURN                                                          SLAP........922000
      ENDIF                                                              SLAP........922100
C                                                                        SLAP........922200
C       RECORD THE MESSAGE.                                              SLAP........922300
C                                                                        SLAP........922400
      I = J4SAVE (1, NERR, .TRUE.)                                       SLAP........922500
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)         SLAP........922600
C                                                                        SLAP........922700
C       HANDLE PRINT-ONCE WARNING MESSAGES.                              SLAP........922800
C                                                                        SLAP........922900
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN                           SLAP........923000
C                                                                        SLAP........923100
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.               SLAP........923200
C                                                                        SLAP........923300
      XLIBR  = LIBRAR                                                    SLAP........923400
      XSUBR  = SUBROU                                                    SLAP........923500
      LFIRST = MESSG                                                     SLAP........923600
      LERR   = NERR                                                      SLAP........923700
      LLEVEL = LEVEL                                                     SLAP........923800
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)           SLAP........923900
C                                                                        SLAP........924000
      LKNTRL = MAX(-2, MIN(2,LKNTRL))                                    SLAP........924100
      MKNTRL = ABS(LKNTRL)                                               SLAP........924200
C                                                                        SLAP........924300
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS    SLAP........924400
C       ZERO AND THE ERROR IS NOT FATAL.                                 SLAP........924500
C                                                                        SLAP........924600
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30                         SLAP........924700
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30                     SLAP........924800
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30   SLAP........924900
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30              SLAP........925000
C                                                                        SLAP........925100
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A   SLAP........925200
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS) SLAP........925300
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG       SLAP........925400
C       IS NOT ZERO.                                                     SLAP........925500
C                                                                        SLAP........925600
      IF (LKNTRL .NE. 0) THEN                                            SLAP........925700
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '                            SLAP........925800
         I = MIN(LEN(SUBROU), 16)                                        SLAP........925900
         TEMP(22:21+I) = SUBROU(1:I)                                     SLAP........926000
         TEMP(22+I:33+I) = ' IN LIBRARY '                                SLAP........926100
         LTEMP = 33 + I                                                  SLAP........926200
         I = MIN(LEN(LIBRAR), 16)                                        SLAP........926300
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)                            SLAP........926400
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'                                 SLAP........926500
         LTEMP = LTEMP + I + 1                                           SLAP........926600
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)                     SLAP........926700
      ENDIF                                                              SLAP........926800
C                                                                        SLAP........926900
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE         SLAP........927000
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE    SLAP........927100
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.                        SLAP........927200
C       1.  LEVEL OF THE MESSAGE                                         SLAP........927300
C              'INFORMATIVE MESSAGE'                                     SLAP........927400
C              'POTENTIALLY RECOVERABLE ERROR'                           SLAP........927500
C              'FATAL ERROR'                                             SLAP........927600
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE          SLAP........927700
C              'PROG CONTINUES'                                          SLAP........927800
C              'PROG ABORTED'                                            SLAP........927900
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK    SLAP........928000
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS     SLAP........928100
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)                 SLAP........928200
C              'TRACEBACK REQUESTED'                                     SLAP........928300
C              'TRACEBACK NOT REQUESTED'                                 SLAP........928400
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT   SLAP........928500
C       EXCEED 74 CHARACTERS.                                            SLAP........928600
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.   SLAP........928700
C                                                                        SLAP........928800
      IF (LKNTRL .GT. 0) THEN                                            SLAP........928900
C                                                                        SLAP........929000
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.             SLAP........929100
C                                                                        SLAP........929200
         IF (LEVEL .LE. 0) THEN                                          SLAP........929300
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'                          SLAP........929400
            LTEMP = 20                                                   SLAP........929500
         ELSEIF (LEVEL .EQ. 1) THEN                                      SLAP........929600
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'                SLAP........929700
            LTEMP = 30                                                   SLAP........929800
         ELSE                                                            SLAP........929900
            TEMP(1:12) = 'FATAL ERROR,'                                  SLAP........930000
            LTEMP = 12                                                   SLAP........930100
         ENDIF                                                           SLAP........930200
C                                                                        SLAP........930300
C       THEN WHETHER THE PROGRAM WILL CONTINUE.                          SLAP........930400
C                                                                        SLAP........930500
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.                         SLAP........930600
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN                        SLAP........930700
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'                    SLAP........930800
            LTEMP = LTEMP + 14                                           SLAP........930900
         ELSE                                                            SLAP........931000
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'                  SLAP........931100
            LTEMP = LTEMP + 16                                           SLAP........931200
         ENDIF                                                           SLAP........931300
C                                                                        SLAP........931400
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.                SLAP........931500
C                                                                        SLAP........931600
         IF (LKNTRL .GT. 0) THEN                                         SLAP........931700
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'              SLAP........931800
            LTEMP = LTEMP + 20                                           SLAP........931900
         ELSE                                                            SLAP........932000
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'          SLAP........932100
            LTEMP = LTEMP + 24                                           SLAP........932200
         ENDIF                                                           SLAP........932300
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)                     SLAP........932400
      ENDIF                                                              SLAP........932500
C                                                                        SLAP........932600
C       NOW SEND OUT THE MESSAGE.                                        SLAP........932700
C                                                                        SLAP........932800
      CALL XERPRN (' *  ', -1, MESSG, 72)                                SLAP........932900
C                                                                        SLAP........933000
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A      SLAP........933100
C          TRACEBACK.                                                    SLAP........933200
C                                                                        SLAP........933300
      IF (LKNTRL .GT. 0) THEN                                            SLAP........933400
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR                  SLAP........933500
         DO 10 I=16,22                                                   SLAP........933600
            IF (TEMP(I:I) .NE. ' ') GO TO 20                             SLAP........933700
   10    CONTINUE                                                        SLAP........933800
C                                                                        SLAP........933900
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)          SLAP........934000
         CALL FDUMP                                                      SLAP........934100
      ENDIF                                                              SLAP........934200
C                                                                        SLAP........934300
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE. SLAP........934400
C                                                                        SLAP........934500
      IF (LKNTRL .NE. 0) THEN                                            SLAP........934600
         CALL XERPRN (' *  ', -1, ' ', 72)                               SLAP........934700
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)                  SLAP........934800
         CALL XERPRN ('    ',  0, ' ', 72)                               SLAP........934900
      ENDIF                                                              SLAP........935000
C                                                                        SLAP........935100
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE    SLAP........935200
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.                   SLAP........935300
C                                                                        SLAP........935400
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN         SLAP........935500
C                                                                        SLAP........935600
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A     SLAP........935700
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR       SLAP........935800
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.  SLAP........935900
C                                                                        SLAP........936000
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN                 SLAP........936100
         IF (LEVEL .EQ. 1) THEN                                          SLAP........936200
            CALL XERPRN                                                  SLAP........936300
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)   SLAP........936400
         ELSE                                                            SLAP........936500
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72) SLAP........936600
         ENDIF                                                           SLAP........936700
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)                   SLAP........936800
         CALL XERHLT (' ')                                               SLAP........936900
      ELSE                                                               SLAP........937000
         CALL XERHLT (MESSG)                                             SLAP........937100
      ENDIF                                                              SLAP........937200
      RETURN                                                             SLAP........937300
      END                                                                SLAP........937400
*DECK XERPRN                                                             SLAP........937500
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)                    SLAP........937600
C***BEGIN PROLOGUE  XERPRN                                               SLAP........937700
C***SUBSIDIARY                                                           SLAP........937800
C***PURPOSE  Print error messages processed by XERMSG.                   SLAP........937900
C***LIBRARY   SLATEC (XERROR)                                            SLAP........938000
C***CATEGORY  R3C                                                        SLAP........938100
C***TYPE      ALL (XERPRN-A)                                             SLAP........938200
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR                           SLAP........938300
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)                                SLAP........938400
C***DESCRIPTION                                                          SLAP........938500
C                                                                        SLAP........938600
C This routine sends one or more lines to each of the (up to five)       SLAP........938700
C logical units to which error messages are to be sent.  This routine    SLAP........938800
C is called several times by XERMSG, sometimes with a single line to     SLAP........938900
C print and sometimes with a (potentially very long) message that may    SLAP........939000
C wrap around into multiple lines.                                       SLAP........939100
C                                                                        SLAP........939200
C PREFIX  Input argument of type CHARACTER.  This argument contains      SLAP........939300
C         characters to be put at the beginning of each line before      SLAP........939400
C         the body of the message.  No more than 16 characters of        SLAP........939500
C         PREFIX will be used.                                           SLAP........939600
C                                                                        SLAP........939700
C NPREF   Input argument of type INTEGER.  This argument is the number   SLAP........939800
C         of characters to use from PREFIX.  If it is negative, the      SLAP........939900
C         intrinsic function LEN is used to determine its length.  If    SLAP........940000
C         it is zero, PREFIX is not used.  If it exceeds 16 or if        SLAP........940100
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be   SLAP........940200
C         used.  If NPREF is positive and the length of PREFIX is less   SLAP........940300
C         than NPREF, a copy of PREFIX extended with blanks to length    SLAP........940400
C         NPREF will be used.                                            SLAP........940500
C                                                                        SLAP........940600
C MESSG   Input argument of type CHARACTER.  This is the text of a       SLAP........940700
C         message to be printed.  If it is a long message, it will be    SLAP........940800
C         broken into pieces for printing on multiple lines.  Each line  SLAP........940900
C         will start with the appropriate prefix and be followed by a    SLAP........941000
C         piece of the message.  NWRAP is the number of characters per   SLAP........941100
C         piece; that is, after each NWRAP characters, we break and      SLAP........941200
C         start a new line.  In addition the characters '$$' embedded    SLAP........941300
C         in MESSG are a sentinel for a new line.  The counting of       SLAP........941400
C         characters up to NWRAP starts over for each new line.  The     SLAP........941500
C         value of NWRAP typically used by XERMSG is 72 since many       SLAP........941600
C         older error messages in the SLATEC Library are laid out to     SLAP........941700
C         rely on wrap-around every 72 characters.                       SLAP........941800
C                                                                        SLAP........941900
C NWRAP   Input argument of type INTEGER.  This gives the maximum size   SLAP........942000
C         piece into which to break MESSG for printing on multiple       SLAP........942100
C         lines.  An embedded '$$' ends a line, and the count restarts   SLAP........942200
C         at the following character.  If a line break does not occur    SLAP........942300
C         on a blank (it would split a word) that word is moved to the   SLAP........942400
C         next line.  Values of NWRAP less than 16 will be treated as    SLAP........942500
C         16.  Values of NWRAP greater than 132 will be treated as 132.  SLAP........942600
C         The actual line length will be NPREF + NWRAP after NPREF has   SLAP........942700
C         been adjusted to fall between 0 and 16 and NWRAP has been      SLAP........942800
C         adjusted to fall between 16 and 132.                           SLAP........942900
C                                                                        SLAP........943000
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........943100
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........943200
C                 Laboratories, 1982.                                    SLAP........943300
C***ROUTINES CALLED  I1MACH, XGETUA                                      SLAP........943400
C***REVISION HISTORY  (YYMMDD)                                           SLAP........943500
C   880621  DATE WRITTEN                                                 SLAP........943600
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF         SLAP........943700
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK    SLAP........943800
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE     SLAP........943900
C           SLASH CHARACTER IN FORMAT STATEMENTS.                        SLAP........944000
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO     SLAP........944100
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK  SLAP........944200
C           LINES TO BE PRINTED.                                         SLAP........944300
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF     SLAP........944400
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.                 SLAP........944500
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.       SLAP........944600
C   891214  Prologue converted to Version 4.0 format.  (WRB)             SLAP........944700
C   900510  Added code to break messages between words.  (RWC)           SLAP........944800
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........944900
C***END PROLOGUE  XERPRN                                                 SLAP........945000
      CHARACTER*(*) PREFIX, MESSG                                        SLAP........945100
      INTEGER NPREF, NWRAP                                               SLAP........945200
      CHARACTER*148 CBUFF                                                SLAP........945300
      INTEGER IU(5), NUNIT                                               SLAP........945400
      CHARACTER*2 NEWLIN                                                 SLAP........945500
      PARAMETER (NEWLIN = '$$')                                          SLAP........945600
C***FIRST EXECUTABLE STATEMENT  XERPRN                                   SLAP........945700
      CALL XGETUA(IU,NUNIT)                                              SLAP........945800
C                                                                        SLAP........945900
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD SLAP........946000
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD    SLAP........946100
C       ERROR MESSAGE UNIT.                                              SLAP........946200
C                                                                        SLAP........946300
      N = I1MACH(4)                                                      SLAP........946400
      DO 10 I=1,NUNIT                                                    SLAP........946500
         IF (IU(I) .EQ. 0) IU(I) = N                                     SLAP........946600
   10 CONTINUE                                                           SLAP........946700
C                                                                        SLAP........946800
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE  SLAP........946900
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING  SLAP........947000
C       THE REST OF THIS ROUTINE.                                        SLAP........947100
C                                                                        SLAP........947200
      IF ( NPREF .LT. 0 ) THEN                                           SLAP........947300
         LPREF = LEN(PREFIX)                                             SLAP........947400
      ELSE                                                               SLAP........947500
         LPREF = NPREF                                                   SLAP........947600
      ENDIF                                                              SLAP........947700
      LPREF = MIN(16, LPREF)                                             SLAP........947800
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX                          SLAP........947900
C                                                                        SLAP........948000
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE SLAP........948100
C       TIME FROM MESSG TO PRINT ON ONE LINE.                            SLAP........948200
C                                                                        SLAP........948300
      LWRAP = MAX(16, MIN(132, NWRAP))                                   SLAP........948400
C                                                                        SLAP........948500
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.   SLAP........948600
C                                                                        SLAP........948700
      LENMSG = LEN(MESSG)                                                SLAP........948800
      N = LENMSG                                                         SLAP........948900
      DO 20 I=1,N                                                        SLAP........949000
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30                     SLAP........949100
         LENMSG = LENMSG - 1                                             SLAP........949200
   20 CONTINUE                                                           SLAP........949300
   30 CONTINUE                                                           SLAP........949400
C                                                                        SLAP........949500
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.         SLAP........949600
C                                                                        SLAP........949700
      IF (LENMSG .EQ. 0) THEN                                            SLAP........949800
         CBUFF(LPREF+1:LPREF+1) = ' '                                    SLAP........949900
         DO 40 I=1,NUNIT                                                 SLAP........950000
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)                         SLAP........950100
   40    CONTINUE                                                        SLAP........950200
         RETURN                                                          SLAP........950300
      ENDIF                                                              SLAP........950400
C                                                                        SLAP........950500
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING      SLAP........950600
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.   SLAP........950700
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.            SLAP........950800
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.     SLAP........950900
C                                                                        SLAP........951000
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE   SLAP........951100
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE  SLAP........951200
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH   SLAP........951300
C       OF THE SECOND ARGUMENT.                                          SLAP........951400
C                                                                        SLAP........951500
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE       SLAP........951600
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER  SLAP........951700
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT        SLAP........951800
C       POSITION NEXTC.                                                  SLAP........951900
C                                                                        SLAP........952000
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE      SLAP........952100
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE       SLAP........952200
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,        SLAP........952300
C                       WHICHEVER IS LESS.                               SLAP........952400
C                                                                        SLAP........952500
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:     SLAP........952600
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE      SLAP........952700
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY     SLAP........952800
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION   SLAP........952900
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF       SLAP........953000
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE     SLAP........953100
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC     SLAP........953200
C                       SHOULD BE INCREMENTED BY 2.                      SLAP........953300
C                                                                        SLAP........953400
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.                     SLAP........953500
C                                                                        SLAP........953600
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1  SLAP........953700
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS         SLAP........953800
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.  SLAP........953900
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY    SLAP........954000
C                       AT THE END OF A LINE.                            SLAP........954100
C                                                                        SLAP........954200
      NEXTC = 1                                                          SLAP........954300
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)                        SLAP........954400
      IF (LPIECE .EQ. 0) THEN                                            SLAP........954500
C                                                                        SLAP........954600
C       THERE WAS NO NEW LINE SENTINEL FOUND.                            SLAP........954700
C                                                                        SLAP........954800
         IDELTA = 0                                                      SLAP........954900
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)                             SLAP........955000
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN                            SLAP........955100
            DO 52 I=LPIECE+1,2,-1                                        SLAP........955200
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN             SLAP........955300
                  LPIECE = I-1                                           SLAP........955400
                  IDELTA = 1                                             SLAP........955500
                  GOTO 54                                                SLAP........955600
               ENDIF                                                     SLAP........955700
   52       CONTINUE                                                     SLAP........955800
         ENDIF                                                           SLAP........955900
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)       SLAP........956000
         NEXTC = NEXTC + LPIECE + IDELTA                                 SLAP........956100
      ELSEIF (LPIECE .EQ. 1) THEN                                        SLAP........956200
C                                                                        SLAP........956300
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).             SLAP........956400
C       DON'T PRINT A BLANK LINE.                                        SLAP........956500
C                                                                        SLAP........956600
         NEXTC = NEXTC + 2                                               SLAP........956700
         GO TO 50                                                        SLAP........956800
      ELSEIF (LPIECE .GT. LWRAP+1) THEN                                  SLAP........956900
C                                                                        SLAP........957000
C       LPIECE SHOULD BE SET DOWN TO LWRAP.                              SLAP........957100
C                                                                        SLAP........957200
         IDELTA = 0                                                      SLAP........957300
         LPIECE = LWRAP                                                  SLAP........957400
         DO 56 I=LPIECE+1,2,-1                                           SLAP........957500
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN                SLAP........957600
               LPIECE = I-1                                              SLAP........957700
               IDELTA = 1                                                SLAP........957800
               GOTO 58                                                   SLAP........957900
            ENDIF                                                        SLAP........958000
   56    CONTINUE                                                        SLAP........958100
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)       SLAP........958200
         NEXTC = NEXTC + LPIECE + IDELTA                                 SLAP........958300
      ELSE                                                               SLAP........958400
C                                                                        SLAP........958500
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.          SLAP........958600
C       WE SHOULD DECREMENT LPIECE BY ONE.                               SLAP........958700
C                                                                        SLAP........958800
         LPIECE = LPIECE - 1                                             SLAP........958900
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)       SLAP........959000
         NEXTC  = NEXTC + LPIECE + 2                                     SLAP........959100
      ENDIF                                                              SLAP........959200
C                                                                        SLAP........959300
C       PRINT                                                            SLAP........959400
C                                                                        SLAP........959500
      DO 60 I=1,NUNIT                                                    SLAP........959600
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)                       SLAP........959700
   60 CONTINUE                                                           SLAP........959800
C                                                                        SLAP........959900
      IF (NEXTC .LE. LENMSG) GO TO 50                                    SLAP........960000
      RETURN                                                             SLAP........960100
      END                                                                SLAP........960200
*DECK XERSVE                                                             SLAP........960300
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,      SLAP........960400
     +   ICOUNT)                                                         SLAP........960500
C***BEGIN PROLOGUE  XERSVE                                               SLAP........960600
C***SUBSIDIARY                                                           SLAP........960700
C***PURPOSE  Record that an error has occurred.                          SLAP........960800
C***LIBRARY   SLATEC (XERROR)                                            SLAP........960900
C***CATEGORY  R3                                                         SLAP........961000
C***TYPE      ALL (XERSVE-A)                                             SLAP........961100
C***KEYWORDS  ERROR, XERROR                                              SLAP........961200
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........961300
C***DESCRIPTION                                                          SLAP........961400
C                                                                        SLAP........961500
C *Usage:                                                                SLAP........961600
C                                                                        SLAP........961700
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT                             SLAP........961800
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG                         SLAP........961900
C                                                                        SLAP........962000
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT) SLAP........962100
C                                                                        SLAP........962200
C *Arguments:                                                            SLAP........962300
C                                                                        SLAP........962400
C        LIBRAR :IN    is the library that the message is from.          SLAP........962500
C        SUBROU :IN    is the subroutine that the message is from.       SLAP........962600
C        MESSG  :IN    is the message to be saved.                       SLAP........962700
C        KFLAG  :IN    indicates the action to be performed.             SLAP........962800
C                      when KFLAG > 0, the message in MESSG is saved.    SLAP........962900
C                      when KFLAG=0 the tables will be dumped and        SLAP........963000
C                      cleared.                                          SLAP........963100
C                      when KFLAG < 0, the tables will be dumped and     SLAP........963200
C                      not cleared.                                      SLAP........963300
C        NERR   :IN    is the error number.                              SLAP........963400
C        LEVEL  :IN    is the error severity.                            SLAP........963500
C        ICOUNT :OUT   the number of times this message has been seen,   SLAP........963600
C                      or zero if the table has overflowed and does not  SLAP........963700
C                      contain this message specifically.  When KFLAG=0, SLAP........963800
C                      ICOUNT will not be altered.                       SLAP........963900
C                                                                        SLAP........964000
C *Description:                                                          SLAP........964100
C                                                                        SLAP........964200
C   Record that this error occurred and possibly dump and clear the      SLAP........964300
C   tables.                                                              SLAP........964400
C                                                                        SLAP........964500
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........964600
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........964700
C                 Laboratories, 1982.                                    SLAP........964800
C***ROUTINES CALLED  I1MACH, XGETUA                                      SLAP........964900
C***REVISION HISTORY  (YYMMDD)                                           SLAP........965000
C   800319  DATE WRITTEN                                                 SLAP........965100
C   861211  REVISION DATE from Version 3.2                               SLAP........965200
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........965300
C   900413  Routine modified to remove reference to KFLAG.  (WRB)        SLAP........965400
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling        SLAP........965500
C           sequence, use IF-THEN-ELSE, make number of saved entries     SLAP........965600
C           easily changeable, changed routine name from XERSAV to       SLAP........965700
C           XERSVE.  (RWC)                                               SLAP........965800
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)            SLAP........965900
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........966000
C***END PROLOGUE  XERSVE                                                 SLAP........966100
      PARAMETER (LENTAB=10)                                              SLAP........966200
      INTEGER LUN(5)                                                     SLAP........966300
      CHARACTER*(*) LIBRAR, SUBROU, MESSG                                SLAP........966400
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB              SLAP........966500
      CHARACTER*20 MESTAB(LENTAB), MES                                   SLAP........966600
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)            SLAP........966700
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG   SLAP........966800
      DATA KOUNTX/0/, NMSG/0/                                            SLAP........966900
C***FIRST EXECUTABLE STATEMENT  XERSVE                                   SLAP........967000
C                                                                        SLAP........967100
      IF (KFLAG.LE.0) THEN                                               SLAP........967200
C                                                                        SLAP........967300
C        Dump the table.                                                 SLAP........967400
C                                                                        SLAP........967500
         IF (NMSG.EQ.0) RETURN                                           SLAP........967600
C                                                                        SLAP........967700
C        Print to each unit.                                             SLAP........967800
C                                                                        SLAP........967900
         CALL XGETUA (LUN, NUNIT)                                        SLAP........968000
         DO 20 KUNIT = 1,NUNIT                                           SLAP........968100
            IUNIT = LUN(KUNIT)                                           SLAP........968200
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)                            SLAP........968300
C                                                                        SLAP........968400
C           Print the table header.                                      SLAP........968500
C                                                                        SLAP........968600
            WRITE (IUNIT,9000)                                           SLAP........968700
C                                                                        SLAP........968800
C           Print body of table.                                         SLAP........968900
C                                                                        SLAP........969000
            DO 10 I = 1,NMSG                                             SLAP........969100
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),       SLAP........969200
     *            NERTAB(I),LEVTAB(I),KOUNT(I)                           SLAP........969300
   10       CONTINUE                                                     SLAP........969400
C                                                                        SLAP........969500
C           Print number of other errors.                                SLAP........969600
C                                                                        SLAP........969700
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX                   SLAP........969800
            WRITE (IUNIT,9030)                                           SLAP........969900
   20    CONTINUE                                                        SLAP........970000
C                                                                        SLAP........970100
C        Clear the error tables.                                         SLAP........970200
C                                                                        SLAP........970300
         IF (KFLAG.EQ.0) THEN                                            SLAP........970400
            NMSG = 0                                                     SLAP........970500
            KOUNTX = 0                                                   SLAP........970600
         ENDIF                                                           SLAP........970700
      ELSE                                                               SLAP........970800
C                                                                        SLAP........970900
C        PROCESS A MESSAGE...                                            SLAP........971000
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,    SLAP........971100
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.                 SLAP........971200
C                                                                        SLAP........971300
         LIB = LIBRAR                                                    SLAP........971400
         SUB = SUBROU                                                    SLAP........971500
         MES = MESSG                                                     SLAP........971600
         DO 30 I = 1,NMSG                                                SLAP........971700
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.            SLAP........971800
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.            SLAP........971900
     *         LEVEL.EQ.LEVTAB(I)) THEN                                  SLAP........972000
                  KOUNT(I) = KOUNT(I) + 1                                SLAP........972100
                  ICOUNT = KOUNT(I)                                      SLAP........972200
                  RETURN                                                 SLAP........972300
            ENDIF                                                        SLAP........972400
   30    CONTINUE                                                        SLAP........972500
C                                                                        SLAP........972600
         IF (NMSG.LT.LENTAB) THEN                                        SLAP........972700
C                                                                        SLAP........972800
C           Empty slot found for new message.                            SLAP........972900
C                                                                        SLAP........973000
            NMSG = NMSG + 1                                              SLAP........973100
            LIBTAB(I) = LIB                                              SLAP........973200
            SUBTAB(I) = SUB                                              SLAP........973300
            MESTAB(I) = MES                                              SLAP........973400
            NERTAB(I) = NERR                                             SLAP........973500
            LEVTAB(I) = LEVEL                                            SLAP........973600
            KOUNT (I) = 1                                                SLAP........973700
            ICOUNT    = 1                                                SLAP........973800
         ELSE                                                            SLAP........973900
C                                                                        SLAP........974000
C           Table is full.                                               SLAP........974100
C                                                                        SLAP........974200
            KOUNTX = KOUNTX+1                                            SLAP........974300
            ICOUNT = 0                                                   SLAP........974400
         ENDIF                                                           SLAP........974500
      ENDIF                                                              SLAP........974600
      RETURN                                                             SLAP........974700
C                                                                        SLAP........974800
C     Formats.                                                           SLAP........974900
C                                                                        SLAP........975000
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /                       SLAP........975100
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',        SLAP........975200
     +   '     LEVEL     COUNT')                                         SLAP........975300
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)                                       SLAP........975400
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)        SLAP........975500
 9030 FORMAT (1X)                                                        SLAP........975600
      END                                                                SLAP........975700
*DECK XGETUA                                                             SLAP........975800
      SUBROUTINE XGETUA (IUNITA, N)                                      SLAP........975900
C***BEGIN PROLOGUE  XGETUA                                               SLAP........976000
C***PURPOSE  Return unit number(s) to which error messages are being     SLAP........976100
C            sent.                                                       SLAP........976200
C***LIBRARY   SLATEC (XERROR)                                            SLAP........976300
C***CATEGORY  R3C                                                        SLAP........976400
C***TYPE      ALL (XGETUA-A)                                             SLAP........976500
C***KEYWORDS  ERROR, XERROR                                              SLAP........976600
C***AUTHOR  Jones, R. E., (SNLA)                                         SLAP........976700
C***DESCRIPTION                                                          SLAP........976800
C                                                                        SLAP........976900
C     Abstract                                                           SLAP........977000
C        XGETUA may be called to determine the unit number or numbers    SLAP........977100
C        to which error messages are being sent.                         SLAP........977200
C        These unit numbers may have been set by a call to XSETUN,       SLAP........977300
C        or a call to XSETUA, or may be a default value.                 SLAP........977400
C                                                                        SLAP........977500
C     Description of Parameters                                          SLAP........977600
C      --Output--                                                        SLAP........977700
C        IUNIT - an array of one to five unit numbers, depending         SLAP........977800
C                on the value of N.  A value of zero refers to the       SLAP........977900
C                default unit, as defined by the I1MACH machine          SLAP........978000
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are       SLAP........978100
C                defined by XGETUA.  The values of IUNIT(N+1),...,       SLAP........978200
C                IUNIT(5) are not defined (for N .LT. 5) or altered      SLAP........978300
C                in any way by XGETUA.                                   SLAP........978400
C        N     - the number of units to which copies of the              SLAP........978500
C                error messages are being sent.  N will be in the        SLAP........978600
C                range from 1 to 5.                                      SLAP........978700
C                                                                        SLAP........978800
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC        SLAP........978900
C                 Error-handling Package, SAND82-0800, Sandia            SLAP........979000
C                 Laboratories, 1982.                                    SLAP........979100
C***ROUTINES CALLED  J4SAVE                                              SLAP........979200
C***REVISION HISTORY  (YYMMDD)                                           SLAP........979300
C   790801  DATE WRITTEN                                                 SLAP........979400
C   861211  REVISION DATE from Version 3.2                               SLAP........979500
C   891214  Prologue converted to Version 4.0 format.  (BAB)             SLAP........979600
C   920501  Reformatted the REFERENCES section.  (WRB)                   SLAP........979700
C***END PROLOGUE  XGETUA                                                 SLAP........979800
      DIMENSION IUNITA(5)                                                SLAP........979900
C***FIRST EXECUTABLE STATEMENT  XGETUA                                   SLAP........980000
      N = J4SAVE(5,0,.FALSE.)                                            SLAP........980100
      DO 30 I=1,N                                                        SLAP........980200
         INDEX = I+4                                                     SLAP........980300
         IF (I.EQ.1) INDEX = 3                                           SLAP........980400
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)                             SLAP........980500
   30 CONTINUE                                                           SLAP........980600
      RETURN                                                             SLAP........980700
      END                                                                SLAP........980800
