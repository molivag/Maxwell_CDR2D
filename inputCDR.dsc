! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! 2d Convection-Diffusion-Reaction simulation                 !
! Input data file                                             !
!                                   MAOG       Bcn, Dic. 2021 !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

# > > > > > > > Geometry
ElemType = Quadrilateral !Triangle, Quadrilateral
DimPr    = 2             !Dimension del problema 
nelem    = 100           !Number of lnods
nnodes   = 121           !Total number of nodal points
nne      = 4             !Number of nodes in the element
ndofn    = 1             !Degrees of freedom
totGp    = 4             !1,4,9 for Q, 1,3,4,6 for P    
maxband  = 26             !Maximo ancho de banda   

# > > > > > > > Stabilization
kstab    = 3             !Type of stabilization 
ktaum    = 1             !Type of tau matrix
patau    = 1.0           !Parameter to obtain tau
hnatu    = 2.0           !Reference element length

# > > > > > > > Physical Properties
#DIFMA_11                        !Diffusion tensor
1.0E-4  , 0.0E-0  , 0.0E-0
0.0E-0  , 1.0E-0  , 0.0E-0
0.0E-0  , 0.0E-0  , 1.0E-0
#DIFMA_12    
0.0e-0 ,  0.0e-0 ,  0.0e-0
0.0e-0  , 0.0e-0  , 0.0e-0
0.0e-0  , 0.0e-0  , 0.0e-0
#DIFMA_22       
1.0E-4 , 0.0E-0 ,  0.0E-0
0.0E-0 , 1.0E-0 ,  0.0E-0
0.0E-0 , 0.0E-0 ,  1.0E-0
#COMAT_1                        !Convection tensor
1.0e+1 , 0.0e-0 , 0.0e-0
0.0e-0 , 4.0e+2 , 0.0e-0
0.0e-0 , 0.0e-0 , 4.0e+2
#COMAT_2      
0.0e+0 , 0.0e-0 , 0.0e-0
0.0e-0 , 0.0e+1 , 0.0e-0
0.0e-0 , 0.0e-0 , 0.0e+3
#REAMA                          !Reaction tensor
0.0e+4 , 0.0e-0 , 0.0e-0
0.0e-0 , 8.0e+4 , 0.0e-0
0.0e-0 , 0.0e-0 , 8.0e+4
#FORCE                          !Force tensor
1.0e+0 , 0.0e+0 , 0.0e+0
1 91 75 73 88	# > > > > > > > Mesh -> Element's Nodes -> Coordinate's Nodes
2 94 78 75 91
3 98 83 78 94
4 102 84 83 98
5 105 96 84 102
6 110 103 96 105
7 115 108 103 110
8 117 112 108 115
9 119 116 112 117
10 121 120 116 119
11 75 60 58 73
12 78 63 60 75
13 83 66 63 78
14 84 71 66 83
15 96 80 71 84
16 103 87 80 96
17 108 100 87 103
18 112 107 100 108
19 116 113 107 112
20 120 118 113 116
21 60 47 45 58
22 63 51 47 60
23 66 54 51 63
24 71 59 54 66
25 80 68 59 71
26 87 79 68 80
27 100 86 79 87
28 107 99 86 100
29 113 109 99 107
30 118 114 109 113
31 47 36 34 45
32 51 38 36 47
33 54 42 38 51
34 59 50 42 54
35 68 56 50 59
36 79 65 56 68
37 86 76 65 79
38 99 90 76 86
39 109 104 90 99
40 114 111 104 109
41 36 28 25 34
42 38 30 28 36
43 42 33 30 38
44 50 41 33 42
45 56 48 41 50
46 65 55 48 56
47 76 69 55 65
48 90 81 69 76
49 104 95 81 90
50 111 106 95 104
51 28 18 17 25
52 30 21 18 28
53 33 23 21 30
54 41 31 23 33
55 48 40 31 41
56 55 49 40 48
57 69 62 49 55
58 81 70 62 69
59 95 85 70 81
60 106 101 85 95
61 18 12 11 17
62 21 14 12 18
63 23 20 14 21
64 31 26 20 23
65 40 32 26 31
66 49 43 32 40
67 62 53 43 49
68 70 67 53 62
69 85 82 67 70
70 101 97 82 85
71 12 8 5 11
72 14 9 8 12
73 20 15 9 14
74 26 22 15 20
75 32 29 22 26
76 43 39 29 32
77 53 52 39 43
78 67 64 52 53
79 82 77 64 67
80 97 93 77 82
81 8 4 3 5
82 9 7 4 8
83 15 13 7 9
84 22 19 13 15
85 29 27 19 22
86 39 37 27 29
87 52 46 37 39
88 64 61 46 52
89 77 74 61 64
90 93 92 74 77
91 4 2 1 3
92 7 6 2 4
93 13 10 6 7
94 19 16 10 13
95 27 24 16 19
96 37 35 24 27
97 46 44 35 37
98 61 57 44 46
99 74 72 57 61
100 92 89 72 74
    1               0              10               0
    2               1              10               0
    3               0               9               0
    4               1               9               0
    5               0               8               0
    6               2              10               0
    7               2               9               0
    8               1               8               0
    9               2               8               0
   10               3              10               0
   11               0               7               0
   12               1               7               0
   13               3               9               0
   14               2               7               0
   15               3               8               0
   16               4              10               0
   17               0               6               0
   18               1               6               0
   19               4               9               0
   20               3               7               0
   21               2               6               0
   22               4               8               0
   23               3               6               0
   24               5              10               0
   25               0               5               0
   26               4               7               0
   27               5               9               0
   28               1               5               0
   29               5               8               0
   30               2               5               0
   31               4               6               0
   32               5               7               0
   33               3               5               0
   34               0               4               0
   35               6              10               0
   36               1               4               0
   37               6               9               0
   38               2               4               0
   39               6               8               0
   40               5               6               0
   41               4               5               0
   42               3               4               0
   43               6               7               0
   44               7              10               0
   45               0               3               0
   46               7               9               0
   47               1               3               0
   48               5               5               0
   49               6               6               0
   50               4               4               0
   51               2               3               0
   52               7               8               0
   53               7               7               0
   54               3               3               0
   55               6               5               0
   56               5               4               0
   57               8              10               0
   58               0               2               0
   59               4               3               0
   60               1               2               0
   61               8               9               0
   62               7               6               0
   63               2               2               0
   64               8               8               0
   65               6               4               0
   66               3               2               0
   67               8               7               0
   68               5               3               0
   69               7               5               0
   70               8               6               0
   71               4               2               0
   72               9              10               0
   73               0               1               0
   74               9               9               0
   75               1               1               0
   76               7               4               0
   77               9               8               0
   78               2               1               0
   79               6               3               0
   80               5               2               0
   81               8               5               0
   82               9               7               0
   83               3               1               0
   84               4               1               0
   85               9               6               0
   86               7               3               0
   87               6               2               0
   88               0               0               0
   89              10              10               0
   90               8               4               0
   91               1               0               0
   92              10               9               0
   93              10               8               0
   94               2               0               0
   95               9               5               0
   96               5               1               0
   97              10               7               0
   98               3               0               0
   99               8               3               0
  100               7               2               0
  101              10               6               0
  102               4               0               0
  103               6               1               0
  104               9               4               0
  105               5               0               0
  106              10               5               0
  107               8               2               0
  108               7               1               0
  109               9               3               0
  110               6               0               0
  111              10               4               0
  112               8               1               0
  113               9               2               0
  114              10               3               0
  115               7               0               0
  116               9               1               0
  117               8               0               0
  118              10               2               0
  119               9               0               0
  120              10               1               0
  121              10               0               0