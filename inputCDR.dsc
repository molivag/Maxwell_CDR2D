! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! 2d Convection-Diffusion-Reaction simulation                 !
! Input data file                                             !
!                                   MAOG       Bcn, Dic. 2021 !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

# > > > > > > > Geometry
ElemType = Quadrilateral !Triangle, Quadrilateral
DimPr    = 2             !Dimension del problema 
nelem    = 4           !Number of lnods
nnodes   = 9           !Total number of nodal points
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
0.0e+1 , 0.0e-0 , 0.0e-0
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
1 2 5 4 1 # > > > > > > > Mesh -> Element's Nodes -> Coordinate's Nodes
2 3 6 5 2
3 5 8 7 4
4 6 9 8 5
1	0.0	0.0
2	0.5	0.0
3	1.0	0.0
4	0.0	0.5
5	0.5	0.5
6	1.0	0.5
7	0.0	1.0
8	0.5	1.0
9	1.0	1.0
