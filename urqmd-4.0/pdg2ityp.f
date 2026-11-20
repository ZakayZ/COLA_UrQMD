      subroutine pdg_to_ityp(pdg_in, ityp_out, iso3_out)

          implicit none
          include 'comres.f'

          integer pdg_in
          integer ityp_out
          integer iso3_out

          integer tab_size
          parameter (TAB_SIZE = 181)

          integer idtab(3, TAB_SIZE)
          integer sort_pdg(TAB_SIZE), sort_orig(TAB_SIZE)
          logical initialized
          save sort_pdg, sort_orig, initialized
          data initialized /.false./

          integer abs_pdg, lo, hi, mid, i, j
          integer tmp_p, tmp_o, ityp_try, iso3_try
          integer, external :: isoit

          data idtab/
c Neutron
     .       1, -1,  2112,
c Proton
     .       1,  1,  2212,
c N*
     .       2, -1, 12112,        2,  1, 12212,
     .       3, -1,  1214,        3,  1,  2124,
     .       4, -1, 22112,        4,  1, 22212,
     .       5, -1, 32112,        5,  1, 32212,
     .       6, -1,  2116,        6,  1,  2216,
     .       7, -1, 12116,        7,  1, 12216,
     .       8, -1, 21214,        8,  1, 22124,
     .       9, -1, 42112,        9,  1, 42212,
     .      10, -1, 31214,       10,  1, 32124,
     .      14, -1,  1218,       14,  1,  2128,
     .      23, -1,  1218,       23,  1,  2128,
c Delta
     .      24, -3,  1114,  24, -1,  2114,  24, 1,  2214,  24, 3,  2224,
     .      25, -3, 31114,  25, -1, 32114,  25, 1, 32214,  25, 3, 32224,
     .      26, -3,  1112,  26, -1,  1212,  26, 1,  2122,  26, 3,  2222,
     .      27, -3, 11114,  27, -1, 12114,  27, 1, 12214,  27, 3, 12224,
     .      28, -3, 11112,  28, -1, 11212,  28, 1, 12122,  28, 3, 12222,
     .      29, -3,  1116,  29, -1,  1216,  29, 1,  2126,  29, 3,  2226,
     .      30, -3, 21112,  30, -1, 21212,  30, 1, 22122,  30, 3, 22222,
     .      31, -3, 21114,  31, -1, 22114,  31, 1, 22214,  31, 3, 22224,
     .      32, -3, 11116,  32, -1, 11216,  32, 1, 12126,  32, 3, 12226,
     .      33, -3,  1118,  33, -1,  2118,  33, 1,  2218,  33, 3,  2228,
     .      40, -3,  1118,  40, -1,  2118,  40, 1,  2218,  40, 3,  2228,
c Lambda
     .      41,  0,  3122,
     .      42,  0, 13122,
     .      43,  0,  3124,
     .      44,  0, 23122,
     .      45,  0, 33122,
     .      46,  0, 13124,
     .      47,  0, 43122,
     .      48,  0, 53122,
     .      49,  0,  3126,
     .      50,  0, 13126,
     .      51,  0, 23124,
     .      52,  0,  3128,
     .      53,  0, 23126,
c Sigma
     .      54, -2,  3112,    54,  0,  3212,    54,  2,  3222,
     .      55, -2,  3114,    55,  0,  3214,    55,  2,  3224,
     .      56, -2, 13112,    56,  0, 13212,    56,  2, 13222,
     .      57, -2, 13114,    57,  0, 13214,    57,  2, 13224,
     .      58, -2, 23112,    58,  0, 23212,    58,  2, 23222,
     .      59, -2,  3116,    59,  0,  3216,    59,  2,  3226,
     .      60, -2, 13116,    60,  0, 13216,    60,  2, 13226,
     .      61, -2, 23114,    61,  0, 23214,    61,  2, 23224,
     .      62, -2,  3118,    62,  0,  3218,    62,  2,  3228,
c Xi
     .      63, -1,  3312,     63,  1,  3322,
     .      64, -1,  3314,     64,  1,  3324,
     .      66, -1, 13314,     66,  1, 13324,
c Omega
     .      69,  0,  3334,
c  Lambda_C
     .      70,  0,  4122,
c gamma
     .     100,  0,    22,
c pion
     .     101, -2,  -211,    101,  0,   111,    101,  2,   211,
c eta
     .     102,  0,   221,
c omega
     .     103,  0,   223,
c rho
     .     104, -2,  -213,    104,  0,   113,    104,  2,   213,
c f0(980)
     .     105,  0, 10221,
c kaon
     .     106, -1,   311,    106,  1,   321,
c eta'
     .     107,  0,   331,
c k*(892)
     .     108, -1,   313,    108,  1,   323,
c phi
     .     109,  0,   333,
c k0*(1430)
     .     110, -1, 10311,    110,  1, 10321,
c a0(980)
     .     111, -2,-10211,    111,  0, 10111,    111,  2, 10211,
c f0(1370)
     .     112,  0, 20221,
c k1(1270)
     .     113, -1, 10313,    113,  1, 10323,
c a1(1260)
     .     114, -2,-20213,    114,  0, 20113,    114,  2, 20213,
c f1(1285)
     .     115,  0, 20223,
c f1'(1510)
     .     116,  0, 40223,
c k2*(1430)
     .     117, -1,   315,    117,  1,   325,
c a2(1329)
     .     118, -2,  -215,    118,  0,   115,    118,  2,   215,
c f2(1270)
     .     119,  0,   225,
c f2'(1525)
     .     120,  0,   335,
c k1(1400)
     .     121, -1, 20313,    121,  1, 20323,
c b1
     .     122, -2,-10213,    122,  0, 10113,    122,  2, 10213,
c h1
     .     123,  0, 10223,
c K* (1410)
     .     125, -1, 30313,    125,  1, 30323,
c rho (1450)
     .     126, -2,-40213,    126,  0, 40113,    126,  2, 40213,
c omega (1420)
     .     127,  0, 50223,
c phi(1680)
     .     128,  0, 10333,
c k*(1680)
     .     129, -1, 40313,    129,  1, 40323,
c rho(1700)
     .     130, -2,-30213,    130,  0, 30113,    130,  2, 30213,
c omega(1600)
     .     131,  0, 60223,
c phi(1850)
     .     132,  0,   337,
c D
     .     133,  -1,   421,    133,   1,   411,
c D*
     .     134,  -1, 10421,    134,   1, 10411,
c J/Psi
     .     135,  0, 443,
c Chi_c
     .     136,  0, 10441,
c Psi'
     .     137,  0, 100443,
c Ds
     .     138,   0,   431,
c Ds*
     .     139,   0,   433 /

          ityp_out = 0
          iso3_out = 0
          if (pdg_in .eq. 0) then
               ityp_out = 11
               iso3_out = 1
               return
          endif
          if (.not. initialized) then
             do i = 1, TAB_SIZE
                sort_pdg(i)  = abs(idtab(3, i))
                sort_orig(i) = i
             enddo
             do i = 2, TAB_SIZE
                tmp_p = sort_pdg(i)
                tmp_o = sort_orig(i)
                j = i - 1
                do while (j .ge. 1)
                   if (sort_pdg(j) .gt. tmp_p) then
                      sort_pdg(j+1)  = sort_pdg(j)
                      sort_orig(j+1) = sort_orig(j)
                      j = j - 1
                   elseif ((sort_pdg(j) .eq. tmp_p) .and.
     &                     (idtab(1, sort_orig(j)) .gt.
     &                      idtab(1, tmp_o))) then
                      sort_pdg(j+1)  = sort_pdg(j)
                      sort_orig(j+1) = sort_orig(j)
                      j = j - 1
                   else
                      exit
                   endif
                enddo
                sort_pdg(j+1)  = tmp_p
                sort_orig(j+1) = tmp_o
             enddo
             initialized = .true.
          endif

c-------------------------------------------------------------------
c  Binary search for |pdg|
c-------------------------------------------------------------------
          abs_pdg = abs(pdg_in)
          lo = 1
          hi = TAB_SIZE
          mid = 0
          do while (lo .le. hi)
             mid = (lo + hi) / 2
             if (sort_pdg(mid) .lt. abs_pdg) then
                lo = mid + 1
             elseif (sort_pdg(mid) .gt. abs_pdg) then
                hi = mid - 1
             else
                exit
             endif
          enddo
          if (sort_pdg(mid) .ne. abs_pdg) goto 200

c  walk back to the first entry with this |pdg|
          do while (mid .gt. 1)
             if (sort_pdg(mid-1) .eq. abs_pdg) then
                mid = mid - 1
             else
                exit
             endif
          enddo


          do i = mid, TAB_SIZE
             if (sort_pdg(i) .ne. abs_pdg) exit
             if (idtab(3, sort_orig(i)) .eq. pdg_in) then
                ityp_out = idtab(1, sort_orig(i))
                iso3_out = idtab(2, sort_orig(i))
                return
             endif
          enddo


          if (pdg_in .lt. 0) then
             do i = mid, TAB_SIZE
                if (sort_pdg(i) .ne. abs_pdg) exit
                if (idtab(3, sort_orig(i)) .eq. abs_pdg) then
                   ityp_try = idtab(1, sort_orig(i))
                   iso3_try = idtab(2, sort_orig(i))
                   if ((ityp_try .ge. minmes) .and.
     &                 (mod(isoit(ityp_try), 2) .eq. 0)) cycle
                   ityp_out = -ityp_try
                   iso3_out = -iso3_try
                   return
                endif
             enddo
          endif

c-------------------------------------------------------------------
c  unknown PDG, encode as |pdg|+1000
c-------------------------------------------------------------------
 200      continue
          ityp_out = sign(1, pdg_in) * (abs_pdg + 1000)
          iso3_out = 0
          return
          end