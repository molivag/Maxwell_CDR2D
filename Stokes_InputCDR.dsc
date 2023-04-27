! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! λ∇x∇xE + ∇P + β∇∇⋅E = f   in  Ω  
! 2d Convection-Diffusion-Reaction simulation                 !          ∇⋅E - ɣΔP  = 0   in  Ω  
! Input data file                                             !                 nxE = 0   on  ∂Ω 
!                                   MAOG       Bcn, Dic. 2021 ! with:                                                             
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       λ, β, ɣ coefficients.    
 
# > > > > > > > Model Parameters
ElemType = QUAD             !Could be QUAD or TRIA
ProbType = TIME             !Problem type TIME=transient, other=static
DimPr    = 2                !Dimension del problema
ndofn    = 3                !Degrees of freedom
totGp    = 4                !1,4,9 for Q, 1,3,7 for P
simul    = 5                !1=LdomT2; 2=LdomT1; 3=SimpleCuad; 4=PolyMaxwell; 5=PolyStokes; 6=SourceLoc
elemSour = 1                !Number of elements will contain the source
skipline = 84               !Lines must be skipped until read the mesh in Geometry module

# > > > > > > >Geometry
nelem    = 400              !Number of nodes
nnodes   = 441              !Total number of nodal points
nne      = 4                !Nodes per element Q:4-9; P:3-6
i_exp    = 0                !Exponent of characteristic mesh size 3,4,5 or 6 2^(-i)
hnatu    = 2.0              !Reference element length
refiType = NO               !NONE; PS=Powell-Sabin; CB=Crossed-Box

# > > > > > > > Time Discretization
theta    = 1                !BDF1=1 ;BDF2=2 ;CN=3
time_ini = 0.0              !Starting time simulation
time_fin = 1.0              !Ending time simulation
max_time = 20               !Max time of simulation
u0cond   = 0.0              !Value of initial condition

# > > > > > > > Stabilization
kstab    = 3                !Stabilization: 0(NONE), 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG), 6(MVAF)
ktaum    = 2                !Tau matrix: 0, 1, 2 
patau    = 1.0              !Parameter to obtain tau
n_val    = 0.0              !n parameter in exact solution, for simul=1
Cu       = 0.0              !Algorithmic constant
ell      = 0.0              !Constante de longitud  
1/mu=λ   = 1.0              !Reluctivity of the medium	µ0=4πE-7 = 795774,71545 [T•m•A^-1]

# > > > > > > > Name outPut Files
pospro   = 2                !Execution of post-processing routine 1=yes, 2=no
testID   = Stokes_testx     !data file with input parameters in each iteration Res/results
Postpro  = Stokes_testx
Error    = xxxxxxxxxxxx
Cordina  = xxxxxxxxxxxx
Conecti  = xxxxxxxxxxxx

# > > > > > > > Physical Properties
#DIFMA_11                !Diffusion tensor
1.0  , 0.0  , 0.0
0.0  , 1.0 , 0.0
0.0  , 0.0  , 0.0
#DIFMA_12         
0.0  , 0.0  , 0.0 
0.0  , 0.0  , 0.0 
0.0  , 0.0  , 0.0 
#DIFMA_21         
0.0  , 0.0  , 0.0 
0.0  , 0.0  , 0.0 
0.0  , 0.0  , 0.0 
#DIFMA_22
1.0 , 0.0  , 0.0
0.0  , 1.0  , 0.0
0.0  , 0.0  , 0.0
#COMAT_1                 !Convection tensor
0.0  , 0.0 , 1.0
0.0  , 0.0 , 0.0
1.0 , 0.0 , 0.0
#COMAT_2
0.0 , 0.0  , 0.0
0.0 , 0.0  , 1.0
0.0 , 1.0 , 0.0
#REAMA                   !Reaction tensor
0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0
#FORCE                   !Force tensor
1.0 , 1.0 , 0.0

# > > > > > > > Element Source Location
3121

# > > > > > > > Mesh -> Element's Nodes -> Coordinate's Nodes 
  1     1    2   23   22
  2     2    3   24   23
  3     3    4   25   24
  4     4    5   26   25
  5     5    6   27   26
  6     6    7   28   27
  7     7    8   29   28
  8     8    9   30   29
  9     9   10   31   30
 10    10   11   32   31
 11    11   12   33   32
 12    12   13   34   33
 13    13   14   35   34
 14    14   15   36   35
 15    15   16   37   36
 16    16   17   38   37
 17    17   18   39   38
 18    18   19   40   39
 19    19   20   41   40
 20    20   21   42   41
 21    22   23   44   43
 22    23   24   45   44
 23    24   25   46   45
 24    25   26   47   46
 25    26   27   48   47
 26    27   28   49   48
 27    28   29   50   49
 28    29   30   51   50
 29    30   31   52   51
 30    31   32   53   52
 31    32   33   54   53
 32    33   34   55   54
 33    34   35   56   55
 34    35   36   57   56
 35    36   37   58   57
 36    37   38   59   58
 37    38   39   60   59
 38    39   40   61   60
 39    40   41   62   61
 40    41   42   63   62
 41    43   44   65   64
 42    44   45   66   65
 43    45   46   67   66
 44    46   47   68   67
 45    47   48   69   68
 46    48   49   70   69
 47    49   50   71   70
 48    50   51   72   71
 49    51   52   73   72
 50    52   53   74   73
 51    53   54   75   74
 52    54   55   76   75
 53    55   56   77   76
 54    56   57   78   77
 55    57   58   79   78
 56    58   59   80   79
 57    59   60   81   80
 58    60   61   82   81
 59    61   62   83   82
 60    62   63   84   83
 61    64   65   86   85
 62    65   66   87   86
 63    66   67   88   87
 64    67   68   89   88
 65    68   69   90   89
 66    69   70   91   90
 67    70   71   92   91
 68    71   72   93   92
 69    72   73   94   93
 70    73   74   95   94
 71    74   75   96   95
 72    75   76   97   96
 73    76   77   98   97
 74    77   78   99   98
 75    78   79  100   99
 76    79   80  101  100
 77    80   81  102  101
 78    81   82  103  102
 79    82   83  104  103
 80    83   84  105  104
 81    85   86  107  106
 82    86   87  108  107
 83    87   88  109  108
 84    88   89  110  109
 85    89   90  111  110
 86    90   91  112  111
 87    91   92  113  112
 88    92   93  114  113
 89    93   94  115  114
 90    94   95  116  115
 91    95   96  117  116
 92    96   97  118  117
 93    97   98  119  118
 94    98   99  120  119
 95    99  100  121  120
 96   100  101  122  121
 97   101  102  123  122
 98   102  103  124  123
 99   103  104  125  124
100   104  105  126  125
101   106  107  128  127
102   107  108  129  128
103   108  109  130  129
104   109  110  131  130
105   110  111  132  131
106   111  112  133  132
107   112  113  134  133
108   113  114  135  134
109   114  115  136  135
110   115  116  137  136
111   116  117  138  137
112   117  118  139  138
113   118  119  140  139
114   119  120  141  140
115   120  121  142  141
116   121  122  143  142
117   122  123  144  143
118   123  124  145  144
119   124  125  146  145
120   125  126  147  146
121   127  128  149  148
122   128  129  150  149
123   129  130  151  150
124   130  131  152  151
125   131  132  153  152
126   132  133  154  153
127   133  134  155  154
128   134  135  156  155
129   135  136  157  156
130   136  137  158  157
131   137  138  159  158
132   138  139  160  159
133   139  140  161  160
134   140  141  162  161
135   141  142  163  162
136   142  143  164  163
137   143  144  165  164
138   144  145  166  165
139   145  146  167  166
140   146  147  168  167
141   148  149  170  169
142   149  150  171  170
143   150  151  172  171
144   151  152  173  172
145   152  153  174  173
146   153  154  175  174
147   154  155  176  175
148   155  156  177  176
149   156  157  178  177
150   157  158  179  178
151   158  159  180  179
152   159  160  181  180
153   160  161  182  181
154   161  162  183  182
155   162  163  184  183
156   163  164  185  184
157   164  165  186  185
158   165  166  187  186
159   166  167  188  187
160   167  168  189  188
161   169  170  191  190
162   170  171  192  191
163   171  172  193  192
164   172  173  194  193
165   173  174  195  194
166   174  175  196  195
167   175  176  197  196
168   176  177  198  197
169   177  178  199  198
170   178  179  200  199
171   179  180  201  200
172   180  181  202  201
173   181  182  203  202
174   182  183  204  203
175   183  184  205  204
176   184  185  206  205
177   185  186  207  206
178   186  187  208  207
179   187  188  209  208
180   188  189  210  209
181   190  191  212  211
182   191  192  213  212
183   192  193  214  213
184   193  194  215  214
185   194  195  216  215
186   195  196  217  216
187   196  197  218  217
188   197  198  219  218
189   198  199  220  219
190   199  200  221  220
191   200  201  222  221
192   201  202  223  222
193   202  203  224  223
194   203  204  225  224
195   204  205  226  225
196   205  206  227  226
197   206  207  228  227
198   207  208  229  228
199   208  209  230  229
200   209  210  231  230
201   211  212  233  232
202   212  213  234  233
203   213  214  235  234
204   214  215  236  235
205   215  216  237  236
206   216  217  238  237
207   217  218  239  238
208   218  219  240  239
209   219  220  241  240
210   220  221  242  241
211   221  222  243  242
212   222  223  244  243
213   223  224  245  244
214   224  225  246  245
215   225  226  247  246
216   226  227  248  247
217   227  228  249  248
218   228  229  250  249
219   229  230  251  250
220   230  231  252  251
221   232  233  254  253
222   233  234  255  254
223   234  235  256  255
224   235  236  257  256
225   236  237  258  257
226   237  238  259  258
227   238  239  260  259
228   239  240  261  260
229   240  241  262  261
230   241  242  263  262
231   242  243  264  263
232   243  244  265  264
233   244  245  266  265
234   245  246  267  266
235   246  247  268  267
236   247  248  269  268
237   248  249  270  269
238   249  250  271  270
239   250  251  272  271
240   251  252  273  272
241   253  254  275  274
242   254  255  276  275
243   255  256  277  276
244   256  257  278  277
245   257  258  279  278
246   258  259  280  279
247   259  260  281  280
248   260  261  282  281
249   261  262  283  282
250   262  263  284  283
251   263  264  285  284
252   264  265  286  285
253   265  266  287  286
254   266  267  288  287
255   267  268  289  288
256   268  269  290  289
257   269  270  291  290
258   270  271  292  291
259   271  272  293  292
260   272  273  294  293
261   274  275  296  295
262   275  276  297  296
263   276  277  298  297
264   277  278  299  298
265   278  279  300  299
266   279  280  301  300
267   280  281  302  301
268   281  282  303  302
269   282  283  304  303
270   283  284  305  304
271   284  285  306  305
272   285  286  307  306
273   286  287  308  307
274   287  288  309  308
275   288  289  310  309
276   289  290  311  310
277   290  291  312  311
278   291  292  313  312
279   292  293  314  313
280   293  294  315  314
281   295  296  317  316
282   296  297  318  317
283   297  298  319  318
284   298  299  320  319
285   299  300  321  320
286   300  301  322  321
287   301  302  323  322
288   302  303  324  323
289   303  304  325  324
290   304  305  326  325
291   305  306  327  326
292   306  307  328  327
293   307  308  329  328
294   308  309  330  329
295   309  310  331  330
296   310  311  332  331
297   311  312  333  332
298   312  313  334  333
299   313  314  335  334
300   314  315  336  335
301   316  317  338  337
302   317  318  339  338
303   318  319  340  339
304   319  320  341  340
305   320  321  342  341
306   321  322  343  342
307   322  323  344  343
308   323  324  345  344
309   324  325  346  345
310   325  326  347  346
311   326  327  348  347
312   327  328  349  348
313   328  329  350  349
314   329  330  351  350
315   330  331  352  351
316   331  332  353  352
317   332  333  354  353
318   333  334  355  354
319   334  335  356  355
320   335  336  357  356
321   337  338  359  358
322   338  339  360  359
323   339  340  361  360
324   340  341  362  361
325   341  342  363  362
326   342  343  364  363
327   343  344  365  364
328   344  345  366  365
329   345  346  367  366
330   346  347  368  367
331   347  348  369  368
332   348  349  370  369
333   349  350  371  370
334   350  351  372  371
335   351  352  373  372
336   352  353  374  373
337   353  354  375  374
338   354  355  376  375
339   355  356  377  376
340   356  357  378  377
341   358  359  380  379
342   359  360  381  380
343   360  361  382  381
344   361  362  383  382
345   362  363  384  383
346   363  364  385  384
347   364  365  386  385
348   365  366  387  386
349   366  367  388  387
350   367  368  389  388
351   368  369  390  389
352   369  370  391  390
353   370  371  392  391
354   371  372  393  392
355   372  373  394  393
356   373  374  395  394
357   374  375  396  395
358   375  376  397  396
359   376  377  398  397
360   377  378  399  398
361   379  380  401  400
362   380  381  402  401
363   381  382  403  402
364   382  383  404  403
365   383  384  405  404
366   384  385  406  405
367   385  386  407  406
368   386  387  408  407
369   387  388  409  408
370   388  389  410  409
371   389  390  411  410
372   390  391  412  411
373   391  392  413  412
374   392  393  414  413
375   393  394  415  414
376   394  395  416  415
377   395  396  417  416
378   396  397  418  417
379   397  398  419  418
380   398  399  420  419
381   400  401  422  421
382   401  402  423  422
383   402  403  424  423
384   403  404  425  424
385   404  405  426  425
386   405  406  427  426
387   406  407  428  427
388   407  408  429  428
389   408  409  430  429
390   409  410  431  430
391   410  411  432  431
392   411  412  433  432
393   412  413  434  433
394   413  414  435  434
395   414  415  436  435
396   415  416  437  436
397   416  417  438  437
398   417  418  439  438
399   418  419  440  439
400   419  420  441  440
  1  	 0.00000   0.00000
  2  	 0.05000   0.00000
  3  	 0.10000   0.00000
  4  	 0.15000   0.00000
  5  	 0.20000   0.00000
  6  	 0.25000   0.00000
  7  	 0.30000   0.00000
  8  	 0.35000   0.00000
  9  	 0.40000   0.00000
 10  	 0.45000   0.00000
 11  	 0.50000   0.00000
 12  	 0.55000   0.00000
 13  	 0.60000   0.00000
 14  	 0.65000   0.00000
 15  	 0.70000   0.00000
 16  	 0.75000   0.00000
 17  	 0.80000   0.00000
 18  	 0.85000   0.00000
 19  	 0.90000   0.00000
 20  	 0.95000   0.00000
 21  	 1.00000   0.00000
 22  	 0.00000   0.05000
 23  	 0.05000   0.05000
 24  	 0.10000   0.05000
 25  	 0.15000   0.05000
 26  	 0.20000   0.05000
 27  	 0.25000   0.05000
 28  	 0.30000   0.05000
 29  	 0.35000   0.05000
 30  	 0.40000   0.05000
 31  	 0.45000   0.05000
 32  	 0.50000   0.05000
 33  	 0.55000   0.05000
 34  	 0.60000   0.05000
 35  	 0.65000   0.05000
 36  	 0.70000   0.05000
 37  	 0.75000   0.05000
 38  	 0.80000   0.05000
 39  	 0.85000   0.05000
 40  	 0.90000   0.05000
 41  	 0.95000   0.05000
 42  	 1.00000   0.05000
 43  	 0.00000   0.10000
 44  	 0.05000   0.10000
 45  	 0.10000   0.10000
 46  	 0.15000   0.10000
 47  	 0.20000   0.10000
 48  	 0.25000   0.10000
 49  	 0.30000   0.10000
 50  	 0.35000   0.10000
 51  	 0.40000   0.10000
 52  	 0.45000   0.10000
 53  	 0.50000   0.10000
 54  	 0.55000   0.10000
 55  	 0.60000   0.10000
 56  	 0.65000   0.10000
 57  	 0.70000   0.10000
 58  	 0.75000   0.10000
 59  	 0.80000   0.10000
 60  	 0.85000   0.10000
 61  	 0.90000   0.10000
 62  	 0.95000   0.10000
 63  	 1.00000   0.10000
 64  	 0.00000   0.15000
 65  	 0.05000   0.15000
 66  	 0.10000   0.15000
 67  	 0.15000   0.15000
 68  	 0.20000   0.15000
 69  	 0.25000   0.15000
 70  	 0.30000   0.15000
 71  	 0.35000   0.15000
 72  	 0.40000   0.15000
 73  	 0.45000   0.15000
 74  	 0.50000   0.15000
 75  	 0.55000   0.15000
 76  	 0.60000   0.15000
 77  	 0.65000   0.15000
 78  	 0.70000   0.15000
 79  	 0.75000   0.15000
 80  	 0.80000   0.15000
 81  	 0.85000   0.15000
 82  	 0.90000   0.15000
 83  	 0.95000   0.15000
 84  	 1.00000   0.15000
 85  	 0.00000   0.20000
 86  	 0.05000   0.20000
 87  	 0.10000   0.20000
 88  	 0.15000   0.20000
 89  	 0.20000   0.20000
 90  	 0.25000   0.20000
 91  	 0.30000   0.20000
 92  	 0.35000   0.20000
 93  	 0.40000   0.20000
 94  	 0.45000   0.20000
 95  	 0.50000   0.20000
 96  	 0.55000   0.20000
 97  	 0.60000   0.20000
 98  	 0.65000   0.20000
 99  	 0.70000   0.20000
100  	 0.75000   0.20000
101  	 0.80000   0.20000
102  	 0.85000   0.20000
103  	 0.90000   0.20000
104  	 0.95000   0.20000
105  	 1.00000   0.20000
106  	 0.00000   0.25000
107  	 0.05000   0.25000
108  	 0.10000   0.25000
109  	 0.15000   0.25000
110  	 0.20000   0.25000
111  	 0.25000   0.25000
112  	 0.30000   0.25000
113  	 0.35000   0.25000
114  	 0.40000   0.25000
115  	 0.45000   0.25000
116  	 0.50000   0.25000
117  	 0.55000   0.25000
118  	 0.60000   0.25000
119  	 0.65000   0.25000
120  	 0.70000   0.25000
121  	 0.75000   0.25000
122  	 0.80000   0.25000
123  	 0.85000   0.25000
124  	 0.90000   0.25000
125  	 0.95000   0.25000
126  	 1.00000   0.25000
127  	 0.00000   0.30000
128  	 0.05000   0.30000
129  	 0.10000   0.30000
130  	 0.15000   0.30000
131  	 0.20000   0.30000
132  	 0.25000   0.30000
133  	 0.30000   0.30000
134  	 0.35000   0.30000
135  	 0.40000   0.30000
136  	 0.45000   0.30000
137  	 0.50000   0.30000
138  	 0.55000   0.30000
139  	 0.60000   0.30000
140  	 0.65000   0.30000
141  	 0.70000   0.30000
142  	 0.75000   0.30000
143  	 0.80000   0.30000
144  	 0.85000   0.30000
145  	 0.90000   0.30000
146  	 0.95000   0.30000
147  	 1.00000   0.30000
148  	 0.00000   0.35000
149  	 0.05000   0.35000
150  	 0.10000   0.35000
151  	 0.15000   0.35000
152  	 0.20000   0.35000
153  	 0.25000   0.35000
154  	 0.30000   0.35000
155  	 0.35000   0.35000
156  	 0.40000   0.35000
157  	 0.45000   0.35000
158  	 0.50000   0.35000
159  	 0.55000   0.35000
160  	 0.60000   0.35000
161  	 0.65000   0.35000
162  	 0.70000   0.35000
163  	 0.75000   0.35000
164  	 0.80000   0.35000
165  	 0.85000   0.35000
166  	 0.90000   0.35000
167  	 0.95000   0.35000
168  	 1.00000   0.35000
169  	 0.00000   0.40000
170  	 0.05000   0.40000
171  	 0.10000   0.40000
172  	 0.15000   0.40000
173  	 0.20000   0.40000
174  	 0.25000   0.40000
175  	 0.30000   0.40000
176  	 0.35000   0.40000
177  	 0.40000   0.40000
178  	 0.45000   0.40000
179  	 0.50000   0.40000
180  	 0.55000   0.40000
181  	 0.60000   0.40000
182  	 0.65000   0.40000
183  	 0.70000   0.40000
184  	 0.75000   0.40000
185  	 0.80000   0.40000
186  	 0.85000   0.40000
187  	 0.90000   0.40000
188  	 0.95000   0.40000
189  	 1.00000   0.40000
190  	 0.00000   0.45000
191  	 0.05000   0.45000
192  	 0.10000   0.45000
193  	 0.15000   0.45000
194  	 0.20000   0.45000
195  	 0.25000   0.45000
196  	 0.30000   0.45000
197  	 0.35000   0.45000
198  	 0.40000   0.45000
199  	 0.45000   0.45000
200  	 0.50000   0.45000
201  	 0.55000   0.45000
202  	 0.60000   0.45000
203  	 0.65000   0.45000
204  	 0.70000   0.45000
205  	 0.75000   0.45000
206  	 0.80000   0.45000
207  	 0.85000   0.45000
208  	 0.90000   0.45000
209  	 0.95000   0.45000
210  	 1.00000   0.45000
211  	 0.00000   0.50000
212  	 0.05000   0.50000
213  	 0.10000   0.50000
214  	 0.15000   0.50000
215  	 0.20000   0.50000
216  	 0.25000   0.50000
217  	 0.30000   0.50000
218  	 0.35000   0.50000
219  	 0.40000   0.50000
220  	 0.45000   0.50000
221  	 0.50000   0.50000
222  	 0.55000   0.50000
223  	 0.60000   0.50000
224  	 0.65000   0.50000
225  	 0.70000   0.50000
226  	 0.75000   0.50000
227  	 0.80000   0.50000
228  	 0.85000   0.50000
229  	 0.90000   0.50000
230  	 0.95000   0.50000
231  	 1.00000   0.50000
232  	 0.00000   0.55000
233  	 0.05000   0.55000
234  	 0.10000   0.55000
235  	 0.15000   0.55000
236  	 0.20000   0.55000
237  	 0.25000   0.55000
238  	 0.30000   0.55000
239  	 0.35000   0.55000
240  	 0.40000   0.55000
241  	 0.45000   0.55000
242  	 0.50000   0.55000
243  	 0.55000   0.55000
244  	 0.60000   0.55000
245  	 0.65000   0.55000
246  	 0.70000   0.55000
247  	 0.75000   0.55000
248  	 0.80000   0.55000
249  	 0.85000   0.55000
250  	 0.90000   0.55000
251  	 0.95000   0.55000
252  	 1.00000   0.55000
253  	 0.00000   0.60000
254  	 0.05000   0.60000
255  	 0.10000   0.60000
256  	 0.15000   0.60000
257  	 0.20000   0.60000
258  	 0.25000   0.60000
259  	 0.30000   0.60000
260  	 0.35000   0.60000
261  	 0.40000   0.60000
262  	 0.45000   0.60000
263  	 0.50000   0.60000
264  	 0.55000   0.60000
265  	 0.60000   0.60000
266  	 0.65000   0.60000
267  	 0.70000   0.60000
268  	 0.75000   0.60000
269  	 0.80000   0.60000
270  	 0.85000   0.60000
271  	 0.90000   0.60000
272  	 0.95000   0.60000
273  	 1.00000   0.60000
274  	 0.00000   0.65000
275  	 0.05000   0.65000
276  	 0.10000   0.65000
277  	 0.15000   0.65000
278  	 0.20000   0.65000
279  	 0.25000   0.65000
280  	 0.30000   0.65000
281  	 0.35000   0.65000
282  	 0.40000   0.65000
283  	 0.45000   0.65000
284  	 0.50000   0.65000
285  	 0.55000   0.65000
286  	 0.60000   0.65000
287  	 0.65000   0.65000
288  	 0.70000   0.65000
289  	 0.75000   0.65000
290  	 0.80000   0.65000
291  	 0.85000   0.65000
292  	 0.90000   0.65000
293  	 0.95000   0.65000
294  	 1.00000   0.65000
295  	 0.00000   0.70000
296  	 0.05000   0.70000
297  	 0.10000   0.70000
298  	 0.15000   0.70000
299  	 0.20000   0.70000
300  	 0.25000   0.70000
301  	 0.30000   0.70000
302  	 0.35000   0.70000
303  	 0.40000   0.70000
304  	 0.45000   0.70000
305  	 0.50000   0.70000
306  	 0.55000   0.70000
307  	 0.60000   0.70000
308  	 0.65000   0.70000
309  	 0.70000   0.70000
310  	 0.75000   0.70000
311  	 0.80000   0.70000
312  	 0.85000   0.70000
313  	 0.90000   0.70000
314  	 0.95000   0.70000
315  	 1.00000   0.70000
316  	 0.00000   0.75000
317  	 0.05000   0.75000
318  	 0.10000   0.75000
319  	 0.15000   0.75000
320  	 0.20000   0.75000
321  	 0.25000   0.75000
322  	 0.30000   0.75000
323  	 0.35000   0.75000
324  	 0.40000   0.75000
325  	 0.45000   0.75000
326  	 0.50000   0.75000
327  	 0.55000   0.75000
328  	 0.60000   0.75000
329  	 0.65000   0.75000
330  	 0.70000   0.75000
331  	 0.75000   0.75000
332  	 0.80000   0.75000
333  	 0.85000   0.75000
334  	 0.90000   0.75000
335  	 0.95000   0.75000
336  	 1.00000   0.75000
337  	 0.00000   0.80000
338  	 0.05000   0.80000
339  	 0.10000   0.80000
340  	 0.15000   0.80000
341  	 0.20000   0.80000
342  	 0.25000   0.80000
343  	 0.30000   0.80000
344  	 0.35000   0.80000
345  	 0.40000   0.80000
346  	 0.45000   0.80000
347  	 0.50000   0.80000
348  	 0.55000   0.80000
349  	 0.60000   0.80000
350  	 0.65000   0.80000
351  	 0.70000   0.80000
352  	 0.75000   0.80000
353  	 0.80000   0.80000
354  	 0.85000   0.80000
355  	 0.90000   0.80000
356  	 0.95000   0.80000
357  	 1.00000   0.80000
358  	 0.00000   0.85000
359  	 0.05000   0.85000
360  	 0.10000   0.85000
361  	 0.15000   0.85000
362  	 0.20000   0.85000
363  	 0.25000   0.85000
364  	 0.30000   0.85000
365  	 0.35000   0.85000
366  	 0.40000   0.85000
367  	 0.45000   0.85000
368  	 0.50000   0.85000
369  	 0.55000   0.85000
370  	 0.60000   0.85000
371  	 0.65000   0.85000
372  	 0.70000   0.85000
373  	 0.75000   0.85000
374  	 0.80000   0.85000
375  	 0.85000   0.85000
376  	 0.90000   0.85000
377  	 0.95000   0.85000
378  	 1.00000   0.85000
379  	 0.00000   0.90000
380  	 0.05000   0.90000
381  	 0.10000   0.90000
382  	 0.15000   0.90000
383  	 0.20000   0.90000
384  	 0.25000   0.90000
385  	 0.30000   0.90000
386  	 0.35000   0.90000
387  	 0.40000   0.90000
388  	 0.45000   0.90000
389  	 0.50000   0.90000
390  	 0.55000   0.90000
391  	 0.60000   0.90000
392  	 0.65000   0.90000
393  	 0.70000   0.90000
394  	 0.75000   0.90000
395  	 0.80000   0.90000
396  	 0.85000   0.90000
397  	 0.90000   0.90000
398  	 0.95000   0.90000
399  	 1.00000   0.90000
400  	 0.00000   0.95000
401  	 0.05000   0.95000
402  	 0.10000   0.95000
403  	 0.15000   0.95000
404  	 0.20000   0.95000
405  	 0.25000   0.95000
406  	 0.30000   0.95000
407  	 0.35000   0.95000
408  	 0.40000   0.95000
409  	 0.45000   0.95000
410  	 0.50000   0.95000
411  	 0.55000   0.95000
412  	 0.60000   0.95000
413  	 0.65000   0.95000
414  	 0.70000   0.95000
415  	 0.75000   0.95000
416  	 0.80000   0.95000
417  	 0.85000   0.95000
418  	 0.90000   0.95000
419  	 0.95000   0.95000
420  	 1.00000   0.95000
421  	 0.00000   1.00000
422  	 0.05000   1.00000
423  	 0.10000   1.00000
424  	 0.15000   1.00000
425  	 0.20000   1.00000
426  	 0.25000   1.00000
427  	 0.30000   1.00000
428  	 0.35000   1.00000
429  	 0.40000   1.00000
430  	 0.45000   1.00000
431  	 0.50000   1.00000
432  	 0.55000   1.00000
433  	 0.60000   1.00000
434  	 0.65000   1.00000
435  	 0.70000   1.00000
436  	 0.75000   1.00000
437  	 0.80000   1.00000
438  	 0.85000   1.00000
439  	 0.90000   1.00000
440  	 0.95000   1.00000
441  	 1.00000   1.00000
