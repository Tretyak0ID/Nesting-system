  <H  Å   k820309              2021.6.0    'p½c                                                                                                          
       src/differential_operators/sbp_differential_operator_mod.f90 SBP_DIFFERENTIAL_OPERATOR_MOD                                                     
       DIFFERENTIAL_OPERATOR_T          @                                         
       FIELD_T          @                                         
       DOMAIN_T                   @                               '                    #NAME    #APPLY                  $                                                        1         Ą                                                 #APPLY_I    #         @                                      	               #THIS    #OUT 	   #IN    #DOMAIN    #DIRECTION              
                                                   #DIFFERENTIAL_OPERATOR_T              
                                	     p               #FIELD_T 
             
                                      p              #FIELD_T 
             
                                      Ų              #DOMAIN_T              
                                                                             @               D           
     'p                    #F    #IS    #IE    #JS    #JE    #INIT    #INIT_ON_DOMAIN    #INIT_REAL     #COPY %   #CREATE_SIMILAR )   #UPDATE_S1 -   #UPDATE_S1V1 2   #UPDATE_S1V1S2V2 8   #UPDATE_S1V1V2 @   #UPDATE G   #ASSIGN_S1 H   #ASSIGN_V1 M   #ASSIGN_S1V1 R   #ASSIGN_S1V1V2 X   #ASSIGN_S1V1S2V2 _   #ASSIGN g                                                                          
            &                   &                                                                                           `                                                             d                                                             h                                                             l             1         Ą    $                                              #INIT    #         @                                                      #THIS    #IS    #IE    #JS    #JE                                                   p               #FIELD_T              
                                                      
                                                      
                                                      
                                            1         Ą    $                                              #INIT_ON_DOMAIN    #         @                                                      #THIS    #DOMAIN                                                   p               #FIELD_T              
                                       Ų              #DOMAIN_T    1         Ą    $                                               #INIT_REAL !   #         @                                  !                    #THIS "   #R #   #DOMAIN $             
                                "     p               #FIELD_T              
                                 #     
                
                                  $     Ų              #DOMAIN_T    1         Ą    $                            %             	     #COPY &   #         @                                  &                    #THIS '   #FIN (             
                                '     p               #FIELD_T              
                                 (     p              #FIELD_T    1         Ą    $                            )             
     #CREATE_SIMILAR *   #         @                                  *                    #THIS +   #DESTINATION ,             
                                 +     p              #FIELD_T              
                                ,     p               #FIELD_T    1         Ą    $                            -                  #UPDATE_FIELD_S1 .   #         @                                  .                    #THIS /   #SCALAR1 0   #DOMAIN 1             
                                /     p               #FIELD_T 
             
                                 0     
                
                                  1     Ų              #DOMAIN_T    1         Ą    $                            2                  #UPDATE_FIELD_S1V1 3   #         @                                  3                    #THIS 4   #SCALAR1 5   #V1 6   #DOMAIN 7             
                                4     p               #FIELD_T 
             
                                 5     
                
                                  6     p              #FIELD_T 
             
                                  7     Ų              #DOMAIN_T    1         Ą    $                            8                  #UPDATE_FIELD_S1V1S2V2 9   #         @                                  9                    #THIS :   #SCALAR1 ;   #V1 <   #SCALAR2 =   #V2 >   #DOMAIN ?             
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
                                  ?     Ų              #DOMAIN_T    1         Ą    $                            @              	    #UPDATE_FIELD_S1V1V2 A   #         @                                  A                    #THIS B   #SCALAR1 C   #F1 D   #F2 E   #DOMAIN F             
                                B     p               #FIELD_T 
             
                                 C     
                
                                  D     p              #FIELD_T 
             
                                  E     p              #FIELD_T 
             
                                  F     Ų              #DOMAIN_T    4             $                         @    G                    3             $                         @             u #FIELD_T 
   #UPDATE_S1 -   #UPDATE_S1V1V2 @   #UPDATE_S1V1 2   #UPDATE_S1V1S2V2 8   1         Ą    $                            H              
    #ASSIGN_FIELD_S1 I   #         @                                  I                    #THIS J   #SCALAR1 K   #DOMAIN L             
                                J     p               #FIELD_T 
             
                                 K     
                
                                  L     Ų              #DOMAIN_T    1         Ą    $                            M                  #ASSIGN_FIELD_V1 N   #         @                                  N                    #THIS O   #V1 P   #DOMAIN Q             
                                O     p               #FIELD_T 
             
                                  P     p              #FIELD_T 
             
                                  Q     Ų              #DOMAIN_T    1         Ą    $                            R                  #ASSIGN_FIELD_S1V1 S   #         @                                  S                    #THIS T   #SCALAR1 U   #V1 V   #DOMAIN W             
                                T     p               #FIELD_T 
             
                                 U     
                
                                  V     p              #FIELD_T 
             
                                  W     Ų              #DOMAIN_T    1         Ą    $                            X                  #ASSIGN_FIELD_S1V1V2 Y   #         @                                  Y                    #THIS Z   #SCALAR1 [   #F1 \   #F2 ]   #DOMAIN ^             
                                Z     p               #FIELD_T 
             
                                 [     
                
                                  \     p              #FIELD_T 
             
                                  ]     p              #FIELD_T 
             
                                  ^     Ų              #DOMAIN_T    1         Ą    $                            _                  #ASSIGN_FIELD_S1V1S2V2 `   #         @                                  `                    #THIS a   #SCALAR1 b   #V1 c   #SCALAR2 d   #V2 e   #DOMAIN f             
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
                                  f     Ų              #DOMAIN_T    4             $                         @    g                    3             $                         @             u #FIELD_T 
   #ASSIGN_S1V1 R   #ASSIGN_S1V1V2 X   #ASSIGN_S1 H   #ASSIGN_V1 M   #ASSIGN_S1V1S2V2 _                     @               @                'p                    #F h   #IS i   #IE j   #JS k   #JE l   #INIT m   #INIT_ON_DOMAIN n   #INIT_REAL o   #COPY p   #CREATE_SIMILAR q   #UPDATE_S1 r   #UPDATE_S1V1 s   #UPDATE_S1V1S2V2 t   #UPDATE_S1V1V2 u   #UPDATE v   #ASSIGN_S1 w   #ASSIGN_V1 x   #ASSIGN_S1V1 y   #ASSIGN_S1V1V2 z   #ASSIGN_S1V1S2V2 {   #ASSIGN |                                            h                              
            &                   &                                                                                      i     `                                                        j     d                                                        k     h                                                        l     l             1         Ą    $                            m                  #INIT    1         Ą    $                            n                  #INIT_ON_DOMAIN    1         Ą    $                            o                  #INIT_REAL !   1         Ą    $                            p             	     #COPY &   1         Ą    $                            q             
     #CREATE_SIMILAR *   1         Ą    $                            r                  #UPDATE_FIELD_S1 .   1         Ą    $                            s                  #UPDATE_FIELD_S1V1 3   1         Ą    $                            t                  #UPDATE_FIELD_S1V1S2V2 9   1         Ą    $                            u              	    #UPDATE_FIELD_S1V1V2 A   4             $                         @    v                    3             $                         @             u #FIELD_T    #UPDATE_S1 r   #UPDATE_S1V1V2 u   #UPDATE_S1V1 s   #UPDATE_S1V1S2V2 t   1         Ą    $                            w              
    #ASSIGN_FIELD_S1 I   1         Ą    $                            x                  #ASSIGN_FIELD_V1 N   1         Ą    $                            y                  #ASSIGN_FIELD_S1V1 S   1         Ą    $                            z                  #ASSIGN_FIELD_S1V1V2 Y   1         Ą    $                            {                  #ASSIGN_FIELD_S1V1S2V2 `   4             $                         @    |                    3             $                         @             u #FIELD_T    #ASSIGN_S1V1 y   #ASSIGN_S1V1V2 z   #ASSIGN_S1 w   #ASSIGN_V1 x   #ASSIGN_S1V1S2V2 {                     @               D                'Ų                    #XS }   #XE ~   #YS    #YE    #DX    #DY    #IS    #IE    #JS    #JE    #NX    #NY    #X    #Y    #INIT                                               }                
                                              ~               
                                                             
                                                             
                                                              
                                                   (          
                                                   0                                                             4                                                             8       	                                                      <       
                                                      @                                                             D                                                                  H                 
            &                                                                                                                 
            &                                           1         Ą                                                  #INIT    #         @                                                   	   #THIS    #XS    #XE    #IS    #IE    #YS    #YE    #JS    #JE                                                   Ų               #DOMAIN_T              
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
                                                              @               @                'Ų                    #XS    #XE    #YS    #YE    #DX    #DY    #IS    #IE    #JS    #JE     #NX ”   #NY ¢   #X £   #Y ¤   #INIT „                                                              
                                                             
                                                             
                                                             
                                                              
                                                   (          
                                                   0                                                             4                                                             8       	                                                       <       
                                                 ”     @                                                        ¢     D                                                      £            H                 
            &                                                                                    ¤                             
            &                                           1         Ą                                „                  #INIT                      @                          ¦     '                    #DIFFERENTIAL_OPERATOR_T §   #APPLY Ø                                               §                           #DIFFERENTIAL_OPERATOR_T    1         Ą                               Ø                  #APPLY_SBP21 ©   #         @                                   ©                    #THIS Ŗ   #OUT «   #IN ¬   #DOMAIN ­   #DIRECTION ®             
                                 Ŗ                   #SBP21_T ¦             
D                                 «     p               #FIELD_T              
                                  ¬     p              #FIELD_T              
                                  ­     Ų              #DOMAIN_T              
                                 ®                                             @                          Æ     '                    #DIFFERENTIAL_OPERATOR_T °   #APPLY ±                                               °                           #DIFFERENTIAL_OPERATOR_T    1         Ą                               ±                  #APPLY_SBP21_2 ²   #         @                                   ²                    #THIS ³   #OUT “   #IN µ   #DOMAIN ¶   #DIRECTION ·             
                                 ³                   #SBP21_2_T Æ             
D                                 “     p               #FIELD_T              
                                  µ     p              #FIELD_T              
                                  ¶     Ų              #DOMAIN_T              
                                 ·                                             @                          ø     '                    #DIFFERENTIAL_OPERATOR_T ¹   #APPLY ŗ                                               ¹                           #DIFFERENTIAL_OPERATOR_T    1         Ą                               ŗ                  #APPLY_SBP42 »   #         @                                   »                    #THIS ¼   #OUT ½   #IN ¾   #DOMAIN æ   #DIRECTION Ą             
                                 ¼                   #SBP42_T ø             
D                                 ½     p               #FIELD_T              
                                  ¾     p              #FIELD_T              
                                  æ     Ų              #DOMAIN_T              
                                 Ą                                  c      fn#fn *     X   J  DIFFERENTIAL_OPERATOR_MOD    [  H   J  FIELD_MOD    £  I   J  DOMAIN_MOD B   ģ  e       DIFFERENTIAL_OPERATOR_T+DIFFERENTIAL_OPERATOR_MOD G   Q  P   a   DIFFERENTIAL_OPERATOR_T%NAME+DIFFERENTIAL_OPERATOR_MOD H   ”  U   a   DIFFERENTIAL_OPERATOR_T%APPLY+DIFFERENTIAL_OPERATOR_MOD 2   ö  ~       APPLY_I+DIFFERENTIAL_OPERATOR_MOD 7   t  e   a   APPLY_I%THIS+DIFFERENTIAL_OPERATOR_MOD 6   Ł  U   a   APPLY_I%OUT+DIFFERENTIAL_OPERATOR_MOD 5   .  U   a   APPLY_I%IN+DIFFERENTIAL_OPERATOR_MOD 9     V   a   APPLY_I%DOMAIN+DIFFERENTIAL_OPERATOR_MOD <   Ł  P   a   APPLY_I%DIRECTION+DIFFERENTIAL_OPERATOR_MOD "   )  y      FIELD_T+FIELD_MOD $   ¢  ¬   a   FIELD_T%F+FIELD_MOD %   N  H   a   FIELD_T%IS+FIELD_MOD %     H   a   FIELD_T%IE+FIELD_MOD %   Ž  H   a   FIELD_T%JS+FIELD_MOD %   &  H   a   FIELD_T%JE+FIELD_MOD '   n  R   a   FIELD_T%INIT+FIELD_MOD    Ą  r       INIT+FIELD_MOD $   2	  U   a   INIT%THIS+FIELD_MOD "   	  @   a   INIT%IS+FIELD_MOD "   Ē	  @   a   INIT%IE+FIELD_MOD "   
  @   a   INIT%JS+FIELD_MOD "   G
  @   a   INIT%JE+FIELD_MOD 1   
  \   a   FIELD_T%INIT_ON_DOMAIN+FIELD_MOD )   ć
  ^       INIT_ON_DOMAIN+FIELD_MOD .   A  U   a   INIT_ON_DOMAIN%THIS+FIELD_MOD 0     V   a   INIT_ON_DOMAIN%DOMAIN+FIELD_MOD ,   ģ  W   a   FIELD_T%INIT_REAL+FIELD_MOD $   C  e       INIT_REAL+FIELD_MOD )   Ø  U   a   INIT_REAL%THIS+FIELD_MOD &   ż  @   a   INIT_REAL%R+FIELD_MOD +   =  V   a   INIT_REAL%DOMAIN+FIELD_MOD '     R   a   FIELD_T%COPY+FIELD_MOD    å  [       COPY+FIELD_MOD $   @  U   a   COPY%THIS+FIELD_MOD #     U   a   COPY%FIN+FIELD_MOD 1   ź  \   a   FIELD_T%CREATE_SIMILAR+FIELD_MOD )   F  c       CREATE_SIMILAR+FIELD_MOD .   ©  U   a   CREATE_SIMILAR%THIS+FIELD_MOD 5   ž  U   a   CREATE_SIMILAR%DESTINATION+FIELD_MOD ,   S  ]   a   FIELD_T%UPDATE_S1+FIELD_MOD *   °  k       UPDATE_FIELD_S1+FIELD_MOD /     U   a   UPDATE_FIELD_S1%THIS+FIELD_MOD 2   p  @   a   UPDATE_FIELD_S1%SCALAR1+FIELD_MOD 1   °  V   a   UPDATE_FIELD_S1%DOMAIN+FIELD_MOD .     _   a   FIELD_T%UPDATE_S1V1+FIELD_MOD ,   e  s       UPDATE_FIELD_S1V1+FIELD_MOD 1   Ų  U   a   UPDATE_FIELD_S1V1%THIS+FIELD_MOD 4   -  @   a   UPDATE_FIELD_S1V1%SCALAR1+FIELD_MOD /   m  U   a   UPDATE_FIELD_S1V1%V1+FIELD_MOD 3   Ā  V   a   UPDATE_FIELD_S1V1%DOMAIN+FIELD_MOD 2     c   a   FIELD_T%UPDATE_S1V1S2V2+FIELD_MOD 0   {         UPDATE_FIELD_S1V1S2V2+FIELD_MOD 5     U   a   UPDATE_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   X  @   a   UPDATE_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3     U   a   UPDATE_FIELD_S1V1S2V2%V1+FIELD_MOD 8   ķ  @   a   UPDATE_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   -  U   a   UPDATE_FIELD_S1V1S2V2%V2+FIELD_MOD 7     V   a   UPDATE_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD 0   Ų  a   a   FIELD_T%UPDATE_S1V1V2+FIELD_MOD .   9  {       UPDATE_FIELD_S1V1V2+FIELD_MOD 3   “  U   a   UPDATE_FIELD_S1V1V2%THIS+FIELD_MOD 6   	  @   a   UPDATE_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   I  U   a   UPDATE_FIELD_S1V1V2%F1+FIELD_MOD 1     U   a   UPDATE_FIELD_S1V1V2%F2+FIELD_MOD 5   ó  V   a   UPDATE_FIELD_S1V1V2%DOMAIN+FIELD_MOD )   I  H   a   FIELD_T%UPDATE+FIELD_MOD 5        `   gen@UPDATE+DIFFERENTIAL_OPERATOR_MOD ,   &  ]   a   FIELD_T%ASSIGN_S1+FIELD_MOD *     k       ASSIGN_FIELD_S1+FIELD_MOD /   ī  U   a   ASSIGN_FIELD_S1%THIS+FIELD_MOD 2   C  @   a   ASSIGN_FIELD_S1%SCALAR1+FIELD_MOD 1     V   a   ASSIGN_FIELD_S1%DOMAIN+FIELD_MOD ,   Ł  ]   a   FIELD_T%ASSIGN_V1+FIELD_MOD *   6  f       ASSIGN_FIELD_V1+FIELD_MOD /     U   a   ASSIGN_FIELD_V1%THIS+FIELD_MOD -   ń  U   a   ASSIGN_FIELD_V1%V1+FIELD_MOD 1   F  V   a   ASSIGN_FIELD_V1%DOMAIN+FIELD_MOD .     _   a   FIELD_T%ASSIGN_S1V1+FIELD_MOD ,   ū  s       ASSIGN_FIELD_S1V1+FIELD_MOD 1   n  U   a   ASSIGN_FIELD_S1V1%THIS+FIELD_MOD 4   Ć  @   a   ASSIGN_FIELD_S1V1%SCALAR1+FIELD_MOD /     U   a   ASSIGN_FIELD_S1V1%V1+FIELD_MOD 3   X  V   a   ASSIGN_FIELD_S1V1%DOMAIN+FIELD_MOD 0   ®  a   a   FIELD_T%ASSIGN_S1V1V2+FIELD_MOD .      {       ASSIGN_FIELD_S1V1V2+FIELD_MOD 3      U   a   ASSIGN_FIELD_S1V1V2%THIS+FIELD_MOD 6   ß   @   a   ASSIGN_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   !  U   a   ASSIGN_FIELD_S1V1V2%F1+FIELD_MOD 1   t!  U   a   ASSIGN_FIELD_S1V1V2%F2+FIELD_MOD 5   É!  V   a   ASSIGN_FIELD_S1V1V2%DOMAIN+FIELD_MOD 2   "  c   a   FIELD_T%ASSIGN_S1V1S2V2+FIELD_MOD 0   "         ASSIGN_FIELD_S1V1S2V2+FIELD_MOD 5   
#  U   a   ASSIGN_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   _#  @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3   #  U   a   ASSIGN_FIELD_S1V1S2V2%V1+FIELD_MOD 8   ō#  @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   4$  U   a   ASSIGN_FIELD_S1V1S2V2%V2+FIELD_MOD 7   $  V   a   ASSIGN_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD )   ß$  H   a   FIELD_T%ASSIGN+FIELD_MOD 5   '%  ¤   `   gen@ASSIGN+DIFFERENTIAL_OPERATOR_MOD "   Ė%  y      FIELD_T+FIELD_MOD $   D'  ¬   a   FIELD_T%F+FIELD_MOD %   š'  H   a   FIELD_T%IS+FIELD_MOD %   8(  H   a   FIELD_T%IE+FIELD_MOD %   (  H   a   FIELD_T%JS+FIELD_MOD %   Č(  H   a   FIELD_T%JE+FIELD_MOD '   )  R   a   FIELD_T%INIT+FIELD_MOD 1   b)  \   a   FIELD_T%INIT_ON_DOMAIN+FIELD_MOD ,   ¾)  W   a   FIELD_T%INIT_REAL+FIELD_MOD '   *  R   a   FIELD_T%COPY+FIELD_MOD 1   g*  \   a   FIELD_T%CREATE_SIMILAR+FIELD_MOD ,   Ć*  ]   a   FIELD_T%UPDATE_S1+FIELD_MOD .    +  _   a   FIELD_T%UPDATE_S1V1+FIELD_MOD 2   +  c   a   FIELD_T%UPDATE_S1V1S2V2+FIELD_MOD 0   ā+  a   a   FIELD_T%UPDATE_S1V1V2+FIELD_MOD )   C,  H   a   FIELD_T%UPDATE+FIELD_MOD %   ,     `   gen@UPDATE+FIELD_MOD ,    -  ]   a   FIELD_T%ASSIGN_S1+FIELD_MOD ,   }-  ]   a   FIELD_T%ASSIGN_V1+FIELD_MOD .   Ś-  _   a   FIELD_T%ASSIGN_S1V1+FIELD_MOD 0   9.  a   a   FIELD_T%ASSIGN_S1V1V2+FIELD_MOD 2   .  c   a   FIELD_T%ASSIGN_S1V1S2V2+FIELD_MOD )   ż.  H   a   FIELD_T%ASSIGN+FIELD_MOD %   E/  ¤   `   gen@ASSIGN+FIELD_MOD $   é/  Č       DOMAIN_T+DOMAIN_MOD '   ±0  H   a   DOMAIN_T%XS+DOMAIN_MOD '   ł0  H   a   DOMAIN_T%XE+DOMAIN_MOD '   A1  H   a   DOMAIN_T%YS+DOMAIN_MOD '   1  H   a   DOMAIN_T%YE+DOMAIN_MOD '   Ń1  H   a   DOMAIN_T%DX+DOMAIN_MOD '   2  H   a   DOMAIN_T%DY+DOMAIN_MOD '   a2  H   a   DOMAIN_T%IS+DOMAIN_MOD '   ©2  H   a   DOMAIN_T%IE+DOMAIN_MOD '   ń2  H   a   DOMAIN_T%JS+DOMAIN_MOD '   93  H   a   DOMAIN_T%JE+DOMAIN_MOD '   3  H   a   DOMAIN_T%NX+DOMAIN_MOD '   É3  H   a   DOMAIN_T%NY+DOMAIN_MOD &   4     a   DOMAIN_T%X+DOMAIN_MOD &   „4     a   DOMAIN_T%Y+DOMAIN_MOD )   95  R   a   DOMAIN_T%INIT+DOMAIN_MOD     5         INIT+DOMAIN_MOD %   6  V   a   INIT%THIS+DOMAIN_MOD #   s6  @   a   INIT%XS+DOMAIN_MOD #   ³6  @   a   INIT%XE+DOMAIN_MOD #   ó6  @   a   INIT%IS+DOMAIN_MOD #   37  @   a   INIT%IE+DOMAIN_MOD #   s7  @   a   INIT%YS+DOMAIN_MOD #   ³7  @   a   INIT%YE+DOMAIN_MOD #   ó7  @   a   INIT%JS+DOMAIN_MOD #   38  @   a   INIT%JE+DOMAIN_MOD $   s8  Č       DOMAIN_T+DOMAIN_MOD '   ;9  H   a   DOMAIN_T%XS+DOMAIN_MOD '   9  H   a   DOMAIN_T%XE+DOMAIN_MOD '   Ė9  H   a   DOMAIN_T%YS+DOMAIN_MOD '   :  H   a   DOMAIN_T%YE+DOMAIN_MOD '   [:  H   a   DOMAIN_T%DX+DOMAIN_MOD '   £:  H   a   DOMAIN_T%DY+DOMAIN_MOD '   ė:  H   a   DOMAIN_T%IS+DOMAIN_MOD '   3;  H   a   DOMAIN_T%IE+DOMAIN_MOD '   {;  H   a   DOMAIN_T%JS+DOMAIN_MOD '   Ć;  H   a   DOMAIN_T%JE+DOMAIN_MOD '   <  H   a   DOMAIN_T%NX+DOMAIN_MOD '   S<  H   a   DOMAIN_T%NY+DOMAIN_MOD &   <     a   DOMAIN_T%X+DOMAIN_MOD &   /=     a   DOMAIN_T%Y+DOMAIN_MOD )   Ć=  R   a   DOMAIN_T%INIT+DOMAIN_MOD    >  x       SBP21_T 0   >  m   a   SBP21_T%DIFFERENTIAL_OPERATOR_T    ś>  Y   a   SBP21_T%APPLY    S?  ~       APPLY_SBP21 !   Ń?  U   a   APPLY_SBP21%THIS     &@  U   a   APPLY_SBP21%OUT    {@  U   a   APPLY_SBP21%IN #   Š@  V   a   APPLY_SBP21%DOMAIN &   &A  P   a   APPLY_SBP21%DIRECTION    vA  x       SBP21_2_T 2   īA  m   a   SBP21_2_T%DIFFERENTIAL_OPERATOR_T     [B  [   a   SBP21_2_T%APPLY    ¶B  ~       APPLY_SBP21_2 #   4C  W   a   APPLY_SBP21_2%THIS "   C  U   a   APPLY_SBP21_2%OUT !   ąC  U   a   APPLY_SBP21_2%IN %   5D  V   a   APPLY_SBP21_2%DOMAIN (   D  P   a   APPLY_SBP21_2%DIRECTION    ŪD  x       SBP42_T 0   SE  m   a   SBP42_T%DIFFERENTIAL_OPERATOR_T    ĄE  Y   a   SBP42_T%APPLY    F  ~       APPLY_SBP42 !   F  U   a   APPLY_SBP42%THIS     ģF  U   a   APPLY_SBP42%OUT    AG  U   a   APPLY_SBP42%IN #   G  V   a   APPLY_SBP42%DOMAIN &   ģG  P   a   APPLY_SBP42%DIRECTION 