  _3  �   k820309              2021.6.0    \jc                                                                                                          
       src/vec_math_mod.f90 VEC_MATH_MOD                                                     
       FIELD_T          @       �                                  
       DOMAIN_T                   @               @                'p                    #F    #IS    #IE    #JS    #JE    #INIT 	   #INIT_ON_DOMAIN    #INIT_REAL    #COPY    #CREATE_SIMILAR    #UPDATE_S1 "   #UPDATE_S1V1 '   #UPDATE_S1V1S2V2 -   #UPDATE_S1V1V2 5   #UPDATE <   #ASSIGN_S1 =   #ASSIGN_V1 B   #ASSIGN_S1V1 G   #ASSIGN_S1V1V2 M   #ASSIGN_S1V1S2V2 T   #ASSIGN \              �                                                            
            &                   &                                                        �                                   `                          �                                   d                          �                                   h                          �                                   l             1         �   � $                      �      	                  #INIT 
   #         @                                  
                    #THIS    #IS    #IE    #JS    #JE                                                   p               #FIELD_T              
                                                      
                                                      
                                                      
                                            1         �   � $                      �                        #INIT_ON_DOMAIN    #         @                                                      #THIS    #DOMAIN                                                   p               #FIELD_T              
                                       �              #DOMAIN_T    1         �   � $                      �                        #INIT_REAL    #         @                                                      #THIS    #R    #DOMAIN              
                                     p               #FIELD_T              
                                      
                
                                       �              #DOMAIN_T    1         �   � $                      �                   	     #COPY    #         @                                                      #THIS    #FIN              
                                     p               #FIELD_T              
                                      p              #FIELD_T    1         �   � $                      �                   
     #CREATE_SIMILAR    #         @                                                      #THIS     #DESTINATION !             
                                       p              #FIELD_T              
                                !     p               #FIELD_T    1         �   � $                      �      "                  #UPDATE_FIELD_S1 #   #         @                                  #                    #THIS $   #SCALAR1 %   #DOMAIN &             
                                $     p               #FIELD_T              
                                 %     
                
                                  &     �              #DOMAIN_T    1         �   � $                      �      '                  #UPDATE_FIELD_S1V1 (   #         @                                  (                    #THIS )   #SCALAR1 *   #V1 +   #DOMAIN ,             
                                )     p               #FIELD_T              
                                 *     
                
                                  +     p              #FIELD_T              
                                  ,     �              #DOMAIN_T    1         �   � $                      �      -                  #UPDATE_FIELD_S1V1S2V2 .   #         @                                  .                    #THIS /   #SCALAR1 0   #V1 1   #SCALAR2 2   #V2 3   #DOMAIN 4             
                                /     p               #FIELD_T              
                                 0     
                
                                  1     p              #FIELD_T              
                                 2     
                
                                  3     p              #FIELD_T              
                                  4     �              #DOMAIN_T    1         �   � $                      �      5              	    #UPDATE_FIELD_S1V1V2 6   #         @                                  6                    #THIS 7   #SCALAR1 8   #F1 9   #F2 :   #DOMAIN ;             
                                7     p               #FIELD_T              
                                 8     
                
                                  9     p              #FIELD_T              
                                  :     p              #FIELD_T              
                                  ;     �              #DOMAIN_T    4         �   � $                         @    <                    3         �   � $                         @             u #FIELD_T    #UPDATE_S1 "   #UPDATE_S1V1V2 5   #UPDATE_S1V1 '   #UPDATE_S1V1S2V2 -   1         �   � $                      �      =              
    #ASSIGN_FIELD_S1 >   #         @                                  >                    #THIS ?   #SCALAR1 @   #DOMAIN A             
                                ?     p               #FIELD_T              
                                 @     
                
                                  A     �              #DOMAIN_T    1         �   � $                      �      B                  #ASSIGN_FIELD_V1 C   #         @                                  C                    #THIS D   #V1 E   #DOMAIN F             
                                D     p               #FIELD_T              
                                  E     p              #FIELD_T              
                                  F     �              #DOMAIN_T    1         �   � $                      �      G                  #ASSIGN_FIELD_S1V1 H   #         @                                  H                    #THIS I   #SCALAR1 J   #V1 K   #DOMAIN L             
                                I     p               #FIELD_T              
                                 J     
                
                                  K     p              #FIELD_T              
                                  L     �              #DOMAIN_T    1         �   � $                      �      M                  #ASSIGN_FIELD_S1V1V2 N   #         @                                  N                    #THIS O   #SCALAR1 P   #F1 Q   #F2 R   #DOMAIN S             
                                O     p               #FIELD_T              
                                 P     
                
                                  Q     p              #FIELD_T              
                                  R     p              #FIELD_T              
                                  S     �              #DOMAIN_T    1         �   � $                      �      T                  #ASSIGN_FIELD_S1V1S2V2 U   #         @                                  U                    #THIS V   #SCALAR1 W   #V1 X   #SCALAR2 Y   #V2 Z   #DOMAIN [             
                                V     p               #FIELD_T              
                                 W     
                
                                  X     p              #FIELD_T              
                                 Y     
                
                                  Z     p              #FIELD_T              
                                  [     �              #DOMAIN_T    4         �   � $                         @    \                    3         �   � $                         @             u #FIELD_T    #ASSIGN_S1V1 G   #ASSIGN_S1V1V2 M   #ASSIGN_S1 =   #ASSIGN_V1 B   #ASSIGN_S1V1S2V2 T                     @               D                '�                    #XS ]   #XE ^   #YS _   #YE `   #DX a   #DY b   #IS c   #IE d   #JS e   #JE f   #NX g   #NY h   #X i   #Y j   #INIT k                �                              ]                
                �                              ^               
                �                              _               
                �                              `               
                �                              a                
                �                              b     (          
                �                              c     0                          �                              d     4                          �                              e     8       	                   �                              f     <       
                   �                              g     @                          �                              h     D                        �                              i            H                 
            &                                                      �                              j            �                 
            &                                           1         �   �                       �      k                  #INIT l   #         @                                  l                 	   #THIS m   #XS o   #XE p   #IS q   #IE r   #YS s   #YE t   #JS u   #JE v                                             m     �               #DOMAIN_T n             
                                 o     
                
                                 p     
                
                                 q                     
                                 r                     
                                 s     
                
                                 t     
                
                                 u                     
                                 v                             @               @           n     '�                    #XS w   #XE x   #YS y   #YE z   #DX {   #DY |   #IS }   #IE ~   #JS    #JE �   #NX �   #NY �   #X �   #Y �   #INIT �                �                              w                
                �                              x               
                �                              y               
                �                              z               
                �                              {                
                �                              |     (          
                �                              }     0                          �                              ~     4                          �                                   8       	                   �                              �     <       
                   �                              �     @                          �                              �     D                        �                              �            H                 
            &                                                      �                              �            �                 
            &                                           1         �   �                       �      �                  #INIT l   %         @                                 �                    
       #FIELD �   #DOMAIN �             
                                  �     p              #FIELD_T              
                                  �     �              #DOMAIN_T n   %         @                                 �                    
       #FIELD �   #DOMAIN �             
                                  �     p              #FIELD_T              
                                  �     �              #DOMAIN_T n   %         @                                 �                    
       #FIELD �   #DOMAIN �             
                                  �     p              #FIELD_T              
                                  �     �              #DOMAIN_T n      �   *      fn#fn    �   H   J  FIELD_MOD      I   J  DOMAIN_MOD "   [  y      FIELD_T+FIELD_MOD $   �  �   a   FIELD_T%F+FIELD_MOD %   �  H   a   FIELD_T%IS+FIELD_MOD %   �  H   a   FIELD_T%IE+FIELD_MOD %     H   a   FIELD_T%JS+FIELD_MOD %   X  H   a   FIELD_T%JE+FIELD_MOD '   �  R   a   FIELD_T%INIT+FIELD_MOD    �  r       INIT+FIELD_MOD $   d  U   a   INIT%THIS+FIELD_MOD "   �  @   a   INIT%IS+FIELD_MOD "   �  @   a   INIT%IE+FIELD_MOD "   9  @   a   INIT%JS+FIELD_MOD "   y  @   a   INIT%JE+FIELD_MOD 1   �  \   a   FIELD_T%INIT_ON_DOMAIN+FIELD_MOD )     ^       INIT_ON_DOMAIN+FIELD_MOD .   s  U   a   INIT_ON_DOMAIN%THIS+FIELD_MOD 0   �  V   a   INIT_ON_DOMAIN%DOMAIN+FIELD_MOD ,     W   a   FIELD_T%INIT_REAL+FIELD_MOD $   u  e       INIT_REAL+FIELD_MOD )   �  U   a   INIT_REAL%THIS+FIELD_MOD &   /	  @   a   INIT_REAL%R+FIELD_MOD +   o	  V   a   INIT_REAL%DOMAIN+FIELD_MOD '   �	  R   a   FIELD_T%COPY+FIELD_MOD    
  [       COPY+FIELD_MOD $   r
  U   a   COPY%THIS+FIELD_MOD #   �
  U   a   COPY%FIN+FIELD_MOD 1     \   a   FIELD_T%CREATE_SIMILAR+FIELD_MOD )   x  c       CREATE_SIMILAR+FIELD_MOD .   �  U   a   CREATE_SIMILAR%THIS+FIELD_MOD 5   0  U   a   CREATE_SIMILAR%DESTINATION+FIELD_MOD ,   �  ]   a   FIELD_T%UPDATE_S1+FIELD_MOD *   �  k       UPDATE_FIELD_S1+FIELD_MOD /   M  U   a   UPDATE_FIELD_S1%THIS+FIELD_MOD 2   �  @   a   UPDATE_FIELD_S1%SCALAR1+FIELD_MOD 1   �  V   a   UPDATE_FIELD_S1%DOMAIN+FIELD_MOD .   8  _   a   FIELD_T%UPDATE_S1V1+FIELD_MOD ,   �  s       UPDATE_FIELD_S1V1+FIELD_MOD 1   
  U   a   UPDATE_FIELD_S1V1%THIS+FIELD_MOD 4   _  @   a   UPDATE_FIELD_S1V1%SCALAR1+FIELD_MOD /   �  U   a   UPDATE_FIELD_S1V1%V1+FIELD_MOD 3   �  V   a   UPDATE_FIELD_S1V1%DOMAIN+FIELD_MOD 2   J  c   a   FIELD_T%UPDATE_S1V1S2V2+FIELD_MOD 0   �  �       UPDATE_FIELD_S1V1S2V2+FIELD_MOD 5   5  U   a   UPDATE_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   �  @   a   UPDATE_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3   �  U   a   UPDATE_FIELD_S1V1S2V2%V1+FIELD_MOD 8     @   a   UPDATE_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   _  U   a   UPDATE_FIELD_S1V1S2V2%V2+FIELD_MOD 7   �  V   a   UPDATE_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD 0   
  a   a   FIELD_T%UPDATE_S1V1V2+FIELD_MOD .   k  {       UPDATE_FIELD_S1V1V2+FIELD_MOD 3   �  U   a   UPDATE_FIELD_S1V1V2%THIS+FIELD_MOD 6   ;  @   a   UPDATE_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   {  U   a   UPDATE_FIELD_S1V1V2%F1+FIELD_MOD 1   �  U   a   UPDATE_FIELD_S1V1V2%F2+FIELD_MOD 5   %  V   a   UPDATE_FIELD_S1V1V2%DOMAIN+FIELD_MOD )   {  H   a   FIELD_T%UPDATE+FIELD_MOD %   �  �   `   gen@UPDATE+FIELD_MOD ,   X  ]   a   FIELD_T%ASSIGN_S1+FIELD_MOD *   �  k       ASSIGN_FIELD_S1+FIELD_MOD /      U   a   ASSIGN_FIELD_S1%THIS+FIELD_MOD 2   u  @   a   ASSIGN_FIELD_S1%SCALAR1+FIELD_MOD 1   �  V   a   ASSIGN_FIELD_S1%DOMAIN+FIELD_MOD ,     ]   a   FIELD_T%ASSIGN_V1+FIELD_MOD *   h  f       ASSIGN_FIELD_V1+FIELD_MOD /   �  U   a   ASSIGN_FIELD_V1%THIS+FIELD_MOD -   #  U   a   ASSIGN_FIELD_V1%V1+FIELD_MOD 1   x  V   a   ASSIGN_FIELD_V1%DOMAIN+FIELD_MOD .   �  _   a   FIELD_T%ASSIGN_S1V1+FIELD_MOD ,   -  s       ASSIGN_FIELD_S1V1+FIELD_MOD 1   �  U   a   ASSIGN_FIELD_S1V1%THIS+FIELD_MOD 4   �  @   a   ASSIGN_FIELD_S1V1%SCALAR1+FIELD_MOD /   5  U   a   ASSIGN_FIELD_S1V1%V1+FIELD_MOD 3   �  V   a   ASSIGN_FIELD_S1V1%DOMAIN+FIELD_MOD 0   �  a   a   FIELD_T%ASSIGN_S1V1V2+FIELD_MOD .   A  {       ASSIGN_FIELD_S1V1V2+FIELD_MOD 3   �  U   a   ASSIGN_FIELD_S1V1V2%THIS+FIELD_MOD 6     @   a   ASSIGN_FIELD_S1V1V2%SCALAR1+FIELD_MOD 1   Q  U   a   ASSIGN_FIELD_S1V1V2%F1+FIELD_MOD 1   �  U   a   ASSIGN_FIELD_S1V1V2%F2+FIELD_MOD 5   �  V   a   ASSIGN_FIELD_S1V1V2%DOMAIN+FIELD_MOD 2   Q  c   a   FIELD_T%ASSIGN_S1V1S2V2+FIELD_MOD 0   �  �       ASSIGN_FIELD_S1V1S2V2+FIELD_MOD 5   <  U   a   ASSIGN_FIELD_S1V1S2V2%THIS+FIELD_MOD 8   �  @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR1+FIELD_MOD 3   �  U   a   ASSIGN_FIELD_S1V1S2V2%V1+FIELD_MOD 8   &   @   a   ASSIGN_FIELD_S1V1S2V2%SCALAR2+FIELD_MOD 3   f   U   a   ASSIGN_FIELD_S1V1S2V2%V2+FIELD_MOD 7   �   V   a   ASSIGN_FIELD_S1V1S2V2%DOMAIN+FIELD_MOD )   !  H   a   FIELD_T%ASSIGN+FIELD_MOD %   Y!  �   `   gen@ASSIGN+FIELD_MOD $   �!  �       DOMAIN_T+DOMAIN_MOD '   �"  H   a   DOMAIN_T%XS+DOMAIN_MOD '   #  H   a   DOMAIN_T%XE+DOMAIN_MOD '   U#  H   a   DOMAIN_T%YS+DOMAIN_MOD '   �#  H   a   DOMAIN_T%YE+DOMAIN_MOD '   �#  H   a   DOMAIN_T%DX+DOMAIN_MOD '   -$  H   a   DOMAIN_T%DY+DOMAIN_MOD '   u$  H   a   DOMAIN_T%IS+DOMAIN_MOD '   �$  H   a   DOMAIN_T%IE+DOMAIN_MOD '   %  H   a   DOMAIN_T%JS+DOMAIN_MOD '   M%  H   a   DOMAIN_T%JE+DOMAIN_MOD '   �%  H   a   DOMAIN_T%NX+DOMAIN_MOD '   �%  H   a   DOMAIN_T%NY+DOMAIN_MOD &   %&  �   a   DOMAIN_T%X+DOMAIN_MOD &   �&  �   a   DOMAIN_T%Y+DOMAIN_MOD )   M'  R   a   DOMAIN_T%INIT+DOMAIN_MOD     �'  �       INIT+DOMAIN_MOD %   1(  V   a   INIT%THIS+DOMAIN_MOD #   �(  @   a   INIT%XS+DOMAIN_MOD #   �(  @   a   INIT%XE+DOMAIN_MOD #   )  @   a   INIT%IS+DOMAIN_MOD #   G)  @   a   INIT%IE+DOMAIN_MOD #   �)  @   a   INIT%YS+DOMAIN_MOD #   �)  @   a   INIT%YE+DOMAIN_MOD #   *  @   a   INIT%JS+DOMAIN_MOD #   G*  @   a   INIT%JE+DOMAIN_MOD $   �*  �       DOMAIN_T+DOMAIN_MOD '   O+  H   a   DOMAIN_T%XS+DOMAIN_MOD '   �+  H   a   DOMAIN_T%XE+DOMAIN_MOD '   �+  H   a   DOMAIN_T%YS+DOMAIN_MOD '   ',  H   a   DOMAIN_T%YE+DOMAIN_MOD '   o,  H   a   DOMAIN_T%DX+DOMAIN_MOD '   �,  H   a   DOMAIN_T%DY+DOMAIN_MOD '   �,  H   a   DOMAIN_T%IS+DOMAIN_MOD '   G-  H   a   DOMAIN_T%IE+DOMAIN_MOD '   �-  H   a   DOMAIN_T%JS+DOMAIN_MOD '   �-  H   a   DOMAIN_T%JE+DOMAIN_MOD '   .  H   a   DOMAIN_T%NX+DOMAIN_MOD '   g.  H   a   DOMAIN_T%NY+DOMAIN_MOD &   �.  �   a   DOMAIN_T%X+DOMAIN_MOD &   C/  �   a   DOMAIN_T%Y+DOMAIN_MOD )   �/  R   a   DOMAIN_T%INIT+DOMAIN_MOD     )0  g       CALC_MASS_FIELD &   �0  U   a   CALC_MASS_FIELD%FIELD '   �0  V   a   CALC_MASS_FIELD%DOMAIN "   ;1  g       CALC_C_NORM_FIELD (   �1  U   a   CALC_C_NORM_FIELD%FIELD )   �1  V   a   CALC_C_NORM_FIELD%DOMAIN (   M2  g       CALC_SQRT_L2_NORM_FIELD .   �2  U   a   CALC_SQRT_L2_NORM_FIELD%FIELD /   	3  V   a   CALC_SQRT_L2_NORM_FIELD%DOMAIN 