  ëD  Œ   k820309              2021.6.0    Êrc                                                                                                          
       src/differential_operators/central_differential_operator_mod.f90 CENTRAL_DIFFERENTIAL_OPERATOR_MOD                                                     
       DIFFERENTIAL_OPERATOR_T          @                                         
       FIELD_T          @                                         
       DOMAIN_T                   @                               '                    #NAME    #APPLY                  $                                                        1         À                                                 #APPLY_I    #         @                                      	               #THIS    #OUT 	   #IN    #DOMAIN    #DIRECTION              
                                                   #DIFFERENTIAL_OPERATOR_T              
                                	     p               #FIELD_T 
             
                                      p              #FIELD_T 
             
                                      Ø              #DOMAIN_T              
                                                                             @               D           
     'p                    #F    #IS    #IE    #JS    #JE    #INIT    #INIT_ON_DOMAIN    #INIT_REAL     #COPY %   #CREATE_SIMILAR )   #UPDATE_S1 -   #UPDATE_S1V1 2   #UPDATE_S1V1S2V2 8   #UPDATE_S1V1V2 @   #UPDATE G   #ASSIGN_S1 H   #ASSIGN_V1 M   #ASSIGN_S1V1 R   #ASSIGN_S1V1V2 X   #ASSIGN_S1V1S2V2 _   #ASSIGN g                                                                          
            &                   &                                                                                           `                                                             d                                                             h                                                             l             1         À    $                                              #INIT    #         @                                                      #THIS    #IS    #IE    #JS    #JE                                                   p               #FIELD_T              
                                                      
                                                      
                                                      
                                            1         À    $                                              #INIT_ON_DOMAIN    #         @                                                      #THIS    #DOMAIN                                                   p               #FIELD_T              
                                       Ø              #DOMAIN_T    1         À    $                                               #INIT_REAL !   #         @                                  !                    #THIS "   #R #   #DOMAIN $             
                                "     p               #FIELD_T              
                                 #     
                
                                  $     Ø              #DOMAIN_T    1         À    $                            %             	     #COPY &   #         @                                  &                    #THIS '   #FIN (             
                                '     p               #FIELD_T              
                                 (     p              #FIELD_T    1         À    $                            )             
     #CREATE_SIMILAR *   #         @                                  *                    #THIS +   #DESTINATION ,             
                                 +     p              #FIELD_T              
                                ,     p               #FIELD_T    1         À    $                            -                  #UPDATE_FIELD_S1 .   #         @                                  .                    #THIS /   #SCALAR1 0   #DOMAIN 1             
                                /     p               #FIELD_T 
             
                                 0     
                
                                  1     Ø              #DOMAIN_T    1         À    $                            2                  #UPDATE_FIELD_S1V1 3   #         @                                  3                    #THIS 4   #SCALAR1 5   #V1 6   #DOMAIN 7             
                                4     p               #FIELD_T 
             
                                 5     
                
                                  6     p              #FIELD_T 
             
                                  7     Ø              #DOMAIN_T    1         À    $                            8                  #UPDATE_FIELD_S1V1S2V2 9   #         @                                  9                    #THIS :   #SCALAR1 ;   #V1 <   #SCALAR2 =   #V2 >   #DOMAIN ?             
                                :     p               #FIELD_T 
             
                                 ;     
                
                                  <     p              #FIELD_T 
             
                                 =     
                
                                  >     p              #FIELD_T 
             
                                  ?     Ø              #DOMAIN_T    1         À    $                            @              	    #UPDATE_FIELD_S1V1V2 A   #         @                                  A                    #THIS B   #SCALAR1 C   #F1 D   #F2 E   #DOMAIN F             
                                B     p               #FIELD_T 
             
                                 C     
                
                                  D     p              #FIELD_T 
             
                                  E     p              #FIELD_T 
             
                                  F     Ø              #DOMAIN_T    4             $                         @    G                    3             $                         @             u #FIELD_T 
   #UPDATE_S1 -   #UPDATE_S1V1V2 @   #UPDATE_S1V1 2   #UPDATE_S1V1S2V2 8   1         À    $                            H              
    #ASSIGN_FIELD_S1 I   #         @                                  I                    #THIS J   #SCALAR1 K   #DOMAIN L             
                                J     p               #FIELD_T 
             
                                 K     
                
                                  L     Ø              #DOMAIN_T    1         À    $                            M                  #ASSIGN_FIELD_V1 N   #         @                                  N                    #THIS O   #V1 P   #DOMAIN Q             
                                O     p               #FIELD_T 
             
                                  P     p              #FIELD_T 
             
                                  Q     Ø              #DOMAIN_T    1         À    $                            R                  #ASSIGN_FIELD_S1V1 S   #         @                                  S                    #THIS T   #SCALAR1 U   #V1 V   #DOMAIN W             
                                T     p               #FIELD_T 
             
                                 U     
                
                                  V     p              #FIELD_T 
             
                                  W     Ø              #DOMAIN_T    1         À    $                            X                  #ASSIGN_FIELD_S1V1V2 Y   #         @                                  Y                    #THIS Z   #SCALAR1 [   #F1 \   #F2 ]   #DOMAIN ^             
                                Z     p               #FIELD_T 
             
                                 [     
                
                                  \     p              #FIELD_T 
             
                                  ]     p              #FIELD_T 
             
                                  ^     Ø              #DOMAIN_T    1         À    $                            _                  #ASSIGN_FIELD_S1V1S2V2 `   #         @                                  `                    #THIS a   #SCALAR1 b   #V1 c   #SCALAR2 d   #V2 e   #DOMAIN f             
                                a     p               #FIELD_T 
             
                                 b     
                
                                  c     p              #FIELD_T 
             
                                 d     
                
                                  e     p              #FIELD_T 
             
                                  f     Ø              #DOMAIN_T    4             $                         @    g                    3             $                         @             u #FIELD_T 
   #ASSIGN_S1V1 R   #ASSIGN_S1V1V2 X   #ASSIGN_S1 H   #ASSIGN_V1 M   #ASSIGN_S1V1S2V2 _                     @               @                'p                    #F h   #IS i   #IE j   #JS k   #JE l   #INIT m   #INIT_ON_DOMAIN n   #INIT_REAL o   #COPY p   #CREATE_SIMILAR q   #UPDATE_S1 r   #UPDATE_S1V1 s   #UPDATE_S1V1S2V2 t   #UPDATE_S1V1V2 u   #UPDATE v   #ASSIGN_S1 w   #ASSIGN_V1 x   #ASSIGN_S1V1 y   #ASSIGN_S1V1V2 z   #ASSIGN_S1V1S2V2 {   #ASSIGN |                                            h                              
            &                   &                                                                                      i     `                                                        j     d                                                        k     h                                                        l     l             1         À    $                            m                  #INIT    1         À    $                            n                  #INIT_ON_DOMAIN    1         À    $                            o                  #INIT_REAL !   1         À    $                            p             	     #COPY &   1         À    $                            q             
     #CREATE_SIMILAR *   1         À    $                            r                  #UPDATE_FIELD_S1 .   1         À    $                            s                  #UPDATE_FIELD_S1V1 3   1         À    $                            t                  #UPDATE_FIELD_S1V1S2V2 9   1         À    $                            u              	    #UPDATE_FIELD_S1V1V2 A   4             $                         @    v                    3             $                         @             u #FIELD_T    #UPDATE_S1 r   #UPDATE_S1V1V2 u   #UPDATE_S1V1 s   #UPDATE_S1V1S2V2 t   1         À    $                            w              
    #ASSIGN_FIELD_S1 I   1         À    $                            x                  #ASSIGN_FIELD_V1 N   1         À    $                            y                  #ASSIGN_FIELD_S1V1 S   1         À    $                            z                  #ASSIGN_FIELD_S1V1V2 Y   1         À    $                            {                  #ASSIGN_FIELD_S1V1S2V2 `   4             $                         @    |                    3             $                         @             u #FIELD_T    #ASSIGN_S1V1 y   #ASSIGN_S1V1V2 z   #ASSIGN_S1 w   #ASSIGN_V1 x   #ASSIGN_S1V1S2V2 {                     @               D                'Ø                    #XS }   #XE ~   #YS    #YE    #DX    #DY    #IS    #IE    #JS    #JE    #NX    #NY    #X    #Y    #INIT                                               }                
                                              ~               
                                                             
                                                             
                                                              
                                                   (          
                                                   0                                                             4                                                             8       	                                                      <       
                                                      @                                                             D                                                                  H                 
            &                                                                                                                 
            &                                           1         À                                                  #INIT    #         @                                                   	   #THIS    #XS    #XE    #IS    #IE    #YS    #YE    #JS    #JE                                                   Ø               #DOMAIN_T              
                                      
                
                                      
                
                                                      
                                                      
                                      
                
                                      
                
                                                      
                                                              @               @                'Ø                    #XS    #XE    #YS    #YE    #DX    #DY    #IS    #IE    #JS    #JE     #NX ¡   #NY ¢   #X £   #Y €   #INIT ¥                                                              
                                                             
                                                             
                                                             
                                                              
                                                   (          
                                                   0                                                             4                                                             8       	                                                       <       
                                                 ¡     @                                                        ¢     D                                                      £            H                 
            &                                                                                    €                             
            &                                           1         À                                ¥                  #INIT                      @                          Š     '                    #DIFFERENTIAL_OPERATOR_T §   #APPLY š                                               §                           #DIFFERENTIAL_OPERATOR_T    1         À                               š                  #APPLY_CENTRAL2 ©   #         @                                   ©                    #THIS ª   #OUT «   #IN ¬   #DOMAIN ­   #DIRECTION ®             
                                 ª                   #CENTRAL2_T Š             
D                                 «     p               #FIELD_T              
                                  ¬     p              #FIELD_T              
                                  ­     Ø              #DOMAIN_T              
                                 ®                                             @                          ¯     '                    #DIFFERENTIAL_OPERATOR_T °   #APPLY ±                                               °                           #DIFFERENTIAL_OPERATOR_T    1         À                               ±                  #APPLY_CENTRAL4 ²   #         @                                   ²                    #THIS ³   #OUT Ž   #IN µ   #DOMAIN ¶   #DIRECTION ·             
                                 ³                   #CENTRAL4_T ¯             
D                                 Ž     p               #FIELD_T              
                                  µ     p              #FIELD_T              
                                  ¶     Ø              #DOMAIN_T              
                                 ·                                  k      fn#fn *     X   J  DIFFERENTIAL_OPERATOR_MOD    c  H   J  FIELD_MOD    «  I   J  DOMAIN_MOD B   ô  e       DIFFERENTIAL_OPERATOR_T+DIFFERENTIAL_OPERATOR_MOD G   Y  P   a   DIFFERENTIAL_OPERATOR_T%NAME+DIFFERENTIAL_OPERATOR_MOD H   ©  U   a   DIFFERENTIAL_OPERATOR_T%APPLY+DIFFERENTIAL_OPERATOR_MOD 2   þ  ~       APPLY_I+DIFFERENTIAL_OPERATOR_MOD 7   |  e   a   APPLY_I%THIS+DIFFERENTIAL_OPERATOR_MOD 6   á  U   a   APPLY_I%OUT+DIFFERENTIAL_OPERATOR_MOD 5   6  U   a   APPLY_I%IN+DIFFERENTIAL_OPERATOR_MOD 9     V   a   APPLY_I%DOMAIN+DIFFERENTIAL_OPERATOR_MOD <   á  P   a   APPLY_I%DIRECTION+DIFFERENTIAL_OPERATOR_MOD "   1  y      FIELD_T+FIELD_MOD $   ª  ¬   a   FIELD_T%F+FIELD_MOD %   V  H   a   FIELD_T%IS+FIELD_MOD %     H   a   FIELD_T%IE+FIELD_MOD %   æ  H   a   FIELD_T%JS+FIELD_MOD %   .  H   a   FIELD_T%JE+FIELD_MOD '   v  R   a   FIELD_T%INIT+FIELD_MOD    È  r       INIT+FIELD_MOD $   :	  U   a   INIT%THIS+FIELD_MOD "   	  @   a   INIT%IS+FIELD_MOD "   Ï	  @   a   INIT%IE+FIELD_MOD "   
  @   a   INIT%JS+FIELD_MOD "   O
  @   a   INIT%JE+FIELD_MOD 1   
  \   a   FIELD_T%INIT_ON_DOMAIN+FIELD_MOD )   ë
  ^       INIT_ON_DOMAIN+FIELD_MOD .   I  U   a   INIT_ON_DOMAIN%THIS+FIELD_MOD 0     V   a   INIT_ON_DOMAIN%DOMAIN+FIELD_MOD ,   ô  W   a   FIELD_T%INIT_REAL+FIELD_MOD $   K  e       INIT_REAL+FIELD_MOD )   °  U   a   INIT_REAL%THIS+FIELD_MOD &     @   a   INIT_REAL%R+FIELD_MOD +   E  V   a   INIT_REAL%DOMAIN+FIELD_MOD '     R   a   FIELD_T%COPY+FIELD_MOD    í  [       COPY+FIELD_MOD $   H  U   a   COPY%THIS+FIELD_MOD #     U   a   COPY%FIN+FIELD_MOD 1   ò  \   a   FIELD_T%CREATE_SIMILAR+FIELD_MOD )   N  c       CREATE_SIMILAR+FIELD_MOD .   ±  U   a   CREATE_SIMILAR%THIS+FIELD_MOD 5     U   a   CREATE_SIMILAR%DESTINATION+FIELD_MOD ,   [  ]   a   FIELD_T%UPDATE_S1+FIELD_MOD *   ž  k       UPDATE_FIELD_S1+FIELD_MOD /   #  U   a   UPDATE_FIELD_S1%THIS+FIELD_MOD 2   x  @   a   UPDATE_FIELD_S1%SCALAR1+FIELD_MOD 1   ž  V   a   UPDATE_FIELD_S1%DOMAIN+FIELD_MOD .     _   a   FIELD_T%UPDATE_S1V1+FIELD_MOD ,   m  s       UPDATE_FIELD_S1V1+FIELD_MOD 1   à  U   a   UPDATE_FIELD_S1V1%THIS+FIELD_MOD 4   5  @   a   UPDATE_FIELD_S1V1%SCALAR1+FIELD_MOD /   u  U   a   UPDATE_FIELD_S1V1%V1+FIELD_MOD 3   Ê  V   a   UPDATE_FIELD_S1V1%DOMAIN+FIELD_MOD 2      c   a   FIELD_T%UPDATE_S1V1S2V2+FIELD_MOD 0            UPDATE_FIELD_S1V1S2V2+FIELD_MOD 5     U   a   UPDATE_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   `  @   a   UPDATE_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3      U   a   UPDATE_FIELD_S1V1S2V2%V1+FIELD_MOD 8   õ  @   a   UPDATE_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   5  U   a   UPDATE_FIELD_S1V1S2V2%V2+FIELD_MOD 7     V   a   UPDATE_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD 0   à  a   a   FIELD_T%UPDATE_S1V1V2+FIELD_MOD .   A  {       UPDATE_FIELD_S1V1V2+FIELD_MOD 3   Œ  U   a   UPDATE_FIELD_S1V1V2%THIS+FIELD_MOD 6     @   a   UPDATE_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   Q  U   a   UPDATE_FIELD_S1V1V2%F1+FIELD_MOD 1   Š  U   a   UPDATE_FIELD_S1V1V2%F2+FIELD_MOD 5   û  V   a   UPDATE_FIELD_S1V1V2%DOMAIN+FIELD_MOD )   Q  H   a   FIELD_T%UPDATE+FIELD_MOD 5        `   gen@UPDATE+DIFFERENTIAL_OPERATOR_MOD ,   .  ]   a   FIELD_T%ASSIGN_S1+FIELD_MOD *     k       ASSIGN_FIELD_S1+FIELD_MOD /   ö  U   a   ASSIGN_FIELD_S1%THIS+FIELD_MOD 2   K  @   a   ASSIGN_FIELD_S1%SCALAR1+FIELD_MOD 1     V   a   ASSIGN_FIELD_S1%DOMAIN+FIELD_MOD ,   á  ]   a   FIELD_T%ASSIGN_V1+FIELD_MOD *   >  f       ASSIGN_FIELD_V1+FIELD_MOD /   €  U   a   ASSIGN_FIELD_V1%THIS+FIELD_MOD -   ù  U   a   ASSIGN_FIELD_V1%V1+FIELD_MOD 1   N  V   a   ASSIGN_FIELD_V1%DOMAIN+FIELD_MOD .   €  _   a   FIELD_T%ASSIGN_S1V1+FIELD_MOD ,     s       ASSIGN_FIELD_S1V1+FIELD_MOD 1   v  U   a   ASSIGN_FIELD_S1V1%THIS+FIELD_MOD 4   Ë  @   a   ASSIGN_FIELD_S1V1%SCALAR1+FIELD_MOD /     U   a   ASSIGN_FIELD_S1V1%V1+FIELD_MOD 3   `  V   a   ASSIGN_FIELD_S1V1%DOMAIN+FIELD_MOD 0   ¶  a   a   FIELD_T%ASSIGN_S1V1V2+FIELD_MOD .      {       ASSIGN_FIELD_S1V1V2+FIELD_MOD 3      U   a   ASSIGN_FIELD_S1V1V2%THIS+FIELD_MOD 6   ç   @   a   ASSIGN_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   '!  U   a   ASSIGN_FIELD_S1V1V2%F1+FIELD_MOD 1   |!  U   a   ASSIGN_FIELD_S1V1V2%F2+FIELD_MOD 5   Ñ!  V   a   ASSIGN_FIELD_S1V1V2%DOMAIN+FIELD_MOD 2   '"  c   a   FIELD_T%ASSIGN_S1V1S2V2+FIELD_MOD 0   "         ASSIGN_FIELD_S1V1S2V2+FIELD_MOD 5   #  U   a   ASSIGN_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   g#  @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3   §#  U   a   ASSIGN_FIELD_S1V1S2V2%V1+FIELD_MOD 8   ü#  @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   <$  U   a   ASSIGN_FIELD_S1V1S2V2%V2+FIELD_MOD 7   $  V   a   ASSIGN_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD )   ç$  H   a   FIELD_T%ASSIGN+FIELD_MOD 5   /%  €   `   gen@ASSIGN+DIFFERENTIAL_OPERATOR_MOD "   Ó%  y      FIELD_T+FIELD_MOD $   L'  ¬   a   FIELD_T%F+FIELD_MOD %   ø'  H   a   FIELD_T%IS+FIELD_MOD %   @(  H   a   FIELD_T%IE+FIELD_MOD %   (  H   a   FIELD_T%JS+FIELD_MOD %   Ð(  H   a   FIELD_T%JE+FIELD_MOD '   )  R   a   FIELD_T%INIT+FIELD_MOD 1   j)  \   a   FIELD_T%INIT_ON_DOMAIN+FIELD_MOD ,   Æ)  W   a   FIELD_T%INIT_REAL+FIELD_MOD '   *  R   a   FIELD_T%COPY+FIELD_MOD 1   o*  \   a   FIELD_T%CREATE_SIMILAR+FIELD_MOD ,   Ë*  ]   a   FIELD_T%UPDATE_S1+FIELD_MOD .   (+  _   a   FIELD_T%UPDATE_S1V1+FIELD_MOD 2   +  c   a   FIELD_T%UPDATE_S1V1S2V2+FIELD_MOD 0   ê+  a   a   FIELD_T%UPDATE_S1V1V2+FIELD_MOD )   K,  H   a   FIELD_T%UPDATE+FIELD_MOD %   ,     `   gen@UPDATE+FIELD_MOD ,   (-  ]   a   FIELD_T%ASSIGN_S1+FIELD_MOD ,   -  ]   a   FIELD_T%ASSIGN_V1+FIELD_MOD .   â-  _   a   FIELD_T%ASSIGN_S1V1+FIELD_MOD 0   A.  a   a   FIELD_T%ASSIGN_S1V1V2+FIELD_MOD 2   ¢.  c   a   FIELD_T%ASSIGN_S1V1S2V2+FIELD_MOD )   /  H   a   FIELD_T%ASSIGN+FIELD_MOD %   M/  €   `   gen@ASSIGN+FIELD_MOD $   ñ/  È       DOMAIN_T+DOMAIN_MOD '   ¹0  H   a   DOMAIN_T%XS+DOMAIN_MOD '   1  H   a   DOMAIN_T%XE+DOMAIN_MOD '   I1  H   a   DOMAIN_T%YS+DOMAIN_MOD '   1  H   a   DOMAIN_T%YE+DOMAIN_MOD '   Ù1  H   a   DOMAIN_T%DX+DOMAIN_MOD '   !2  H   a   DOMAIN_T%DY+DOMAIN_MOD '   i2  H   a   DOMAIN_T%IS+DOMAIN_MOD '   ±2  H   a   DOMAIN_T%IE+DOMAIN_MOD '   ù2  H   a   DOMAIN_T%JS+DOMAIN_MOD '   A3  H   a   DOMAIN_T%JE+DOMAIN_MOD '   3  H   a   DOMAIN_T%NX+DOMAIN_MOD '   Ñ3  H   a   DOMAIN_T%NY+DOMAIN_MOD &   4     a   DOMAIN_T%X+DOMAIN_MOD &   ­4     a   DOMAIN_T%Y+DOMAIN_MOD )   A5  R   a   DOMAIN_T%INIT+DOMAIN_MOD     5         INIT+DOMAIN_MOD %   %6  V   a   INIT%THIS+DOMAIN_MOD #   {6  @   a   INIT%XS+DOMAIN_MOD #   »6  @   a   INIT%XE+DOMAIN_MOD #   û6  @   a   INIT%IS+DOMAIN_MOD #   ;7  @   a   INIT%IE+DOMAIN_MOD #   {7  @   a   INIT%YS+DOMAIN_MOD #   »7  @   a   INIT%YE+DOMAIN_MOD #   û7  @   a   INIT%JS+DOMAIN_MOD #   ;8  @   a   INIT%JE+DOMAIN_MOD $   {8  È       DOMAIN_T+DOMAIN_MOD '   C9  H   a   DOMAIN_T%XS+DOMAIN_MOD '   9  H   a   DOMAIN_T%XE+DOMAIN_MOD '   Ó9  H   a   DOMAIN_T%YS+DOMAIN_MOD '   :  H   a   DOMAIN_T%YE+DOMAIN_MOD '   c:  H   a   DOMAIN_T%DX+DOMAIN_MOD '   «:  H   a   DOMAIN_T%DY+DOMAIN_MOD '   ó:  H   a   DOMAIN_T%IS+DOMAIN_MOD '   ;;  H   a   DOMAIN_T%IE+DOMAIN_MOD '   ;  H   a   DOMAIN_T%JS+DOMAIN_MOD '   Ë;  H   a   DOMAIN_T%JE+DOMAIN_MOD '   <  H   a   DOMAIN_T%NX+DOMAIN_MOD '   [<  H   a   DOMAIN_T%NY+DOMAIN_MOD &   £<     a   DOMAIN_T%X+DOMAIN_MOD &   7=     a   DOMAIN_T%Y+DOMAIN_MOD )   Ë=  R   a   DOMAIN_T%INIT+DOMAIN_MOD    >  x       CENTRAL2_T 3   >  m   a   CENTRAL2_T%DIFFERENTIAL_OPERATOR_T !   ?  \   a   CENTRAL2_T%APPLY    ^?  ~       APPLY_CENTRAL2 $   Ü?  X   a   APPLY_CENTRAL2%THIS #   4@  U   a   APPLY_CENTRAL2%OUT "   @  U   a   APPLY_CENTRAL2%IN &   Þ@  V   a   APPLY_CENTRAL2%DOMAIN )   4A  P   a   APPLY_CENTRAL2%DIRECTION    A  x       CENTRAL4_T 3   üA  m   a   CENTRAL4_T%DIFFERENTIAL_OPERATOR_T !   iB  \   a   CENTRAL4_T%APPLY    ÅB  ~       APPLY_CENTRAL4 $   CC  X   a   APPLY_CENTRAL4%THIS #   C  U   a   APPLY_CENTRAL4%OUT "   ðC  U   a   APPLY_CENTRAL4%IN &   ED  V   a   APPLY_CENTRAL4%DOMAIN )   D  P   a   APPLY_CENTRAL4%DIRECTION 