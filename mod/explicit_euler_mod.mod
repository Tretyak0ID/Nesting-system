  4  �   k820309              2021.6.0    N�c                                                                                                          
       src/timeschemes/explicit_Euler_mod.f90 EXPLICIT_EULER_MOD                                                     
       STVEC_T                                                     
       TIMESCHEME_T          @       �                                  
       OPERATOR_T          @       �                                  
       MULTI_DOMAIN_T                   @                               '                      #COPY    #CREATE_SIMILAR    #UPDATE_S1V1    #UPDATE_S1V1S2V2    #UPDATE_S1V1V2    #UPDATE %   #ASSIGN_S1 &   #ASSIGN_S1V1 +   #ASSIGN_S1V1S2V2 1   #ASSIGN_S1V1V2 9   #ASSIGN @   1         �   � $                      �                        #COPY_STVEC    #         @                                                      #THIS    #FIN 	   #MULTI_DOMAIN 
             
                                                     #STVEC_T              
                                 	                    #STVEC_T              
                                 
     �             #MULTI_DOMAIN_T    1         �   � $                      �                        #CREATE_SIMILAR_STVEC    #         @                                                      #THIS    #DESTINATION              
                                                     #STVEC_T              
                                                     #STVEC_T    1         �   � $                     �                        #UPDATE_STVEC_S1V1    #         @                                                     #THIS    #SCALAR1    #V1    #MULTI_DOMAIN              
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      �             #MULTI_DOMAIN_T    1         �   � $                      �                        #UPDATE_STVEC_S1V1S2V2    #         @                                                      #THIS    #SCALAR1    #V1    #SCALAR2    #V2    #MULTI_DOMAIN              
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      �             #MULTI_DOMAIN_T    1         �   � $                      �                        #UPDATE_STVEC_S1V1V2    #         @                                                      #THIS     #SCALAR1 !   #V1 "   #V2 #   #MULTI_DOMAIN $             
                                                      #STVEC_T              
                                 !     
                
                                 "                    #STVEC_T              
                                 #                    #STVEC_T              
                                 $     �             #MULTI_DOMAIN_T    4         �   � $                         @    %                    3         �   � $                         @             u #STVEC_T    #UPDATE_S1V1    #UPDATE_S1V1V2    #UPDATE_S1V1S2V2    1         �   � $                      �      &                  #ASSIGN_STVEC_S1 '   #         @                                  '                    #THIS (   #SCALAR1 )   #MULTI_DOMAIN *             
                                (                     #STVEC_T              
                                 )     
                
                                 *     �             #MULTI_DOMAIN_T    1         �   � $                      �      +                  #ASSIGN_STVEC_S1V1 ,   #         @                                  ,                    #THIS -   #SCALAR1 .   #V1 /   #MULTI_DOMAIN 0             
                                -                     #STVEC_T              
                                 .     
                
                                 /                    #STVEC_T              
                                 0     �             #MULTI_DOMAIN_T    1         �   � $                      �      1             	     #ASSIGN_STVEC_S1V1S2V2 2   #         @                                  2                    #THIS 3   #SCALAR1 4   #V1 5   #SCALAR2 6   #V2 7   #MULTI_DOMAIN 8             
                                3                     #STVEC_T              
                                 4     
                
                                 5                    #STVEC_T              
                                 6     
                
                                 7                    #STVEC_T              
                                 8     �             #MULTI_DOMAIN_T    1         �   � $                      �      9             
 	    #ASSIGN_STVEC_S1V1V2 :   #         @                                  :                    #THIS ;   #SCALAR1 <   #V1 =   #V2 >   #MULTI_DOMAIN ?             
                                ;                     #STVEC_T              
                                 <     
                
                                 =                    #STVEC_T              
                                 >                    #STVEC_T              
                                 ?     �             #MULTI_DOMAIN_T    4         �   � $                         @    @                    3         �   � $                         @             u #STVEC_T    #ASSIGN_S1V1 +   #ASSIGN_S1 &   #ASSIGN_S1V1S2V2 1   #ASSIGN_S1V1V2 9                     @                          A     '                      #STEP B   1         �   �                       �     B                  #STEP C   #         @                                 C     	               #THIS D   #V0 E   #OPERATOR F   #MULTI_DOMAIN H   #DT I             
                               D                     #TIMESCHEME_T A             
                               E                     #STVEC_T              
                               F                     #OPERATOR_T G             
                                 H     �             #MULTI_DOMAIN_T              
                                I     
                        @                         G     '                      #APPLY J   1         �   �                       �     J                  #APPLY_I K   #         @                                K     	               #THIS L   #OUT M   #IN N   #MULTI_DOMAIN O             
                               L                     #OPERATOR_T G             
                               M                     #STVEC_T              
                               N                     #STVEC_T              
                                 O     �             #MULTI_DOMAIN_T                      @                          P     '                      #APPLY Q   1         �   �                      �     Q                  #APPLY_I K                     @               �                '�                   #GLOBAL_DOMAIN R   #SUBDOMAINS m   #NUM_SUB_X n   #NUM_SUB_Y o   #DEGREE_CONDENSATION p   #INIT q                �                               R     �                      #DOMAIN_T S                     @              D           S     '�                    #XS T   #XE U   #YS V   #YE W   #DX X   #DY Y   #IS Z   #IE [   #JS \   #JE ]   #NX ^   #NY _   #X `   #Y a   #INIT b                �                              T                
                �                              U               
                �                              V               
                �                              W               
                �                              X                
                �                              Y     (          
                �                              Z     0                          �                              [     4                          �                              \     8       	                   �                              ]     <       
                   �                              ^     @                          �                              _     D                        �                              `            H                 
            &                                                      �                              a            �                 
            &                                           1         �   �                       �      b                  #INIT c   #         @                                  c                 	   #THIS d   #XS e   #XE f   #IS g   #IE h   #YS i   #YE j   #JS k   #JE l                                             d     �               #DOMAIN_T S             
                                 e     
                
                                 f     
                
                                 g                     
                                 h                     
                                 i     
                
                                 j     
                
                                 k                     
                                 l                      �                               m            �       �             #DOMAIN_T S             &                   &                                                        �                              n     8                         �                              o     <                       �                              p            @                            &                   &                                           1         �   � $                      �      q                  #INIT_MULTI_DOMAIN r   #         @                                  r                    #THIS s   #GLOBAL_DOMAIN t   #NUM_SUB_X u   #NUM_SUB_Y v   #DEGREE_CONDENSATION w             
                                s     �              #MULTI_DOMAIN_T              
                                 t     �              #DOMAIN_T S             
                                 u                     
                                 v                   
                                 w                                 &                   &                                                             @               �           x     '�                   #GLOBAL_DOMAIN y   #SUBDOMAINS z   #NUM_SUB_X {   #NUM_SUB_Y |   #DEGREE_CONDENSATION }   #INIT ~                �                               y     �                      #DOMAIN_T S              �                               z            �       �             #DOMAIN_T S             &                   &                                                        �                              {     8                         �                              |     <                       �                              }            @                            &                   &                                           1         �   � $                      �      ~                  #INIT_MULTI_DOMAIN r                     @               @              '�                    #TIMESCHEME_T �   #TENDENCY �   #STEP �                �                               �                            #TIMESCHEME_T A               �                              �                            #STVEC_T    1         �   � $                      �     �                  #STEP_EXPLICIT_EULER �   #         @                                   �                    #THIS �   #V0 �   #OPERATOR �   #MULTI_DOMAIN �   #DT �             
D @                              �     �               #EXPLICIT_EULER_T              
D @                              �                     #STVEC_T              
D @                              �                     #OPERATOR_T P             
  @                               �     �             #MULTI_DOMAIN_T x             
  @                              �     
         �   B      fn#fn    �   H   J  STVEC_MOD    *  M   J  TIMESCHEME_MOD    w  K   J  OPERATOR_MOD !   �  O   J  MULTI_DOMAIN_MOD "           STVEC_T+STVEC_MOD '     X   a   STVEC_T%COPY+STVEC_MOD %   p  m       COPY_STVEC+STVEC_MOD *   �  U   a   COPY_STVEC%THIS+STVEC_MOD )   2  U   a   COPY_STVEC%FIN+STVEC_MOD 2   �  \   a   COPY_STVEC%MULTI_DOMAIN+STVEC_MOD 1   �  b   a   STVEC_T%CREATE_SIMILAR+STVEC_MOD /   E  c       CREATE_SIMILAR_STVEC+STVEC_MOD 4   �  U   a   CREATE_SIMILAR_STVEC%THIS+STVEC_MOD ;   �  U   a   CREATE_SIMILAR_STVEC%DESTINATION+STVEC_MOD .   R  _   a   STVEC_T%UPDATE_S1V1+STVEC_MOD ,   �  y       UPDATE_STVEC_S1V1+STVEC_MOD 1   *  U   a   UPDATE_STVEC_S1V1%THIS+STVEC_MOD 4     @   a   UPDATE_STVEC_S1V1%SCALAR1+STVEC_MOD /   �  U   a   UPDATE_STVEC_S1V1%V1+STVEC_MOD 9     \   a   UPDATE_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   p  c   a   STVEC_T%UPDATE_S1V1S2V2+STVEC_MOD 0   �  �       UPDATE_STVEC_S1V1S2V2+STVEC_MOD 5   a	  U   a   UPDATE_STVEC_S1V1S2V2%THIS+STVEC_MOD 8   �	  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   �	  U   a   UPDATE_STVEC_S1V1S2V2%V1+STVEC_MOD 8   K
  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3   �
  U   a   UPDATE_STVEC_S1V1S2V2%V2+STVEC_MOD =   �
  \   a   UPDATE_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   <  a   a   STVEC_T%UPDATE_S1V1V2+STVEC_MOD .   �  �       UPDATE_STVEC_S1V1V2+STVEC_MOD 3     U   a   UPDATE_STVEC_S1V1V2%THIS+STVEC_MOD 6   s  @   a   UPDATE_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1   �  U   a   UPDATE_STVEC_S1V1V2%V1+STVEC_MOD 1     U   a   UPDATE_STVEC_S1V1V2%V2+STVEC_MOD ;   ]  \   a   UPDATE_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )   �  H   a   STVEC_T%UPDATE+STVEC_MOD %     �   `   gen@UPDATE+STVEC_MOD ,   �  ]   a   STVEC_T%ASSIGN_S1+STVEC_MOD *   �  q       ASSIGN_STVEC_S1+STVEC_MOD /   U  U   a   ASSIGN_STVEC_S1%THIS+STVEC_MOD 2   �  @   a   ASSIGN_STVEC_S1%SCALAR1+STVEC_MOD 7   �  \   a   ASSIGN_STVEC_S1%MULTI_DOMAIN+STVEC_MOD .   F  _   a   STVEC_T%ASSIGN_S1V1+STVEC_MOD ,   �  y       ASSIGN_STVEC_S1V1+STVEC_MOD 1     U   a   ASSIGN_STVEC_S1V1%THIS+STVEC_MOD 4   s  @   a   ASSIGN_STVEC_S1V1%SCALAR1+STVEC_MOD /   �  U   a   ASSIGN_STVEC_S1V1%V1+STVEC_MOD 9     \   a   ASSIGN_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   d  c   a   STVEC_T%ASSIGN_S1V1S2V2+STVEC_MOD 0   �  �       ASSIGN_STVEC_S1V1S2V2+STVEC_MOD 5   U  U   a   ASSIGN_STVEC_S1V1S2V2%THIS+STVEC_MOD 8   �  @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   �  U   a   ASSIGN_STVEC_S1V1S2V2%V1+STVEC_MOD 8   ?  @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3     U   a   ASSIGN_STVEC_S1V1S2V2%V2+STVEC_MOD =   �  \   a   ASSIGN_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   0  a   a   STVEC_T%ASSIGN_S1V1V2+STVEC_MOD .   �  �       ASSIGN_STVEC_S1V1V2+STVEC_MOD 3     U   a   ASSIGN_STVEC_S1V1V2%THIS+STVEC_MOD 6   g  @   a   ASSIGN_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1   �  U   a   ASSIGN_STVEC_S1V1V2%V1+STVEC_MOD 1   �  U   a   ASSIGN_STVEC_S1V1V2%V2+STVEC_MOD ;   Q  \   a   ASSIGN_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )   �  H   a   STVEC_T%ASSIGN+STVEC_MOD %   �  �   `   gen@ASSIGN+STVEC_MOD ,   �  Z       TIMESCHEME_T+TIMESCHEME_MOD 1   �  R   a   TIMESCHEME_T%STEP+TIMESCHEME_MOD $   6  �       STEP+TIMESCHEME_MOD )   �  Z   a   STEP%THIS+TIMESCHEME_MOD '     U   a   STEP%V0+TIMESCHEME_MOD -   g  X   a   STEP%OPERATOR+TIMESCHEME_MOD 1   �  \   a   STEP%MULTI_DOMAIN+TIMESCHEME_MOD '     @   a   STEP%DT+TIMESCHEME_MOD (   [  [       OPERATOR_T+OPERATOR_MOD .   �  U   a   OPERATOR_T%APPLY+OPERATOR_MOD %     u       APPLY_I+OPERATOR_MOD *   �  X   a   APPLY_I%THIS+OPERATOR_MOD )   �  U   a   APPLY_I%OUT+OPERATOR_MOD (   -  U   a   APPLY_I%IN+OPERATOR_MOD 2   �  \   a   APPLY_I%MULTI_DOMAIN+OPERATOR_MOD (   �  [       OPERATOR_T+OPERATOR_MOD .   9  U   a   OPERATOR_T%APPLY+OPERATOR_MOD 0   �  �       MULTI_DOMAIN_T+MULTI_DOMAIN_MOD >   B  ^   a   MULTI_DOMAIN_T%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD $   �  �       DOMAIN_T+DOMAIN_MOD '   h   H   a   DOMAIN_T%XS+DOMAIN_MOD '   �   H   a   DOMAIN_T%XE+DOMAIN_MOD '   �   H   a   DOMAIN_T%YS+DOMAIN_MOD '   @!  H   a   DOMAIN_T%YE+DOMAIN_MOD '   �!  H   a   DOMAIN_T%DX+DOMAIN_MOD '   �!  H   a   DOMAIN_T%DY+DOMAIN_MOD '   "  H   a   DOMAIN_T%IS+DOMAIN_MOD '   `"  H   a   DOMAIN_T%IE+DOMAIN_MOD '   �"  H   a   DOMAIN_T%JS+DOMAIN_MOD '   �"  H   a   DOMAIN_T%JE+DOMAIN_MOD '   8#  H   a   DOMAIN_T%NX+DOMAIN_MOD '   �#  H   a   DOMAIN_T%NY+DOMAIN_MOD &   �#  �   a   DOMAIN_T%X+DOMAIN_MOD &   \$  �   a   DOMAIN_T%Y+DOMAIN_MOD )   �$  R   a   DOMAIN_T%INIT+DOMAIN_MOD     B%  �       INIT+DOMAIN_MOD %   �%  V   a   INIT%THIS+DOMAIN_MOD #   *&  @   a   INIT%XS+DOMAIN_MOD #   j&  @   a   INIT%XE+DOMAIN_MOD #   �&  @   a   INIT%IS+DOMAIN_MOD #   �&  @   a   INIT%IE+DOMAIN_MOD #   *'  @   a   INIT%YS+DOMAIN_MOD #   j'  @   a   INIT%YE+DOMAIN_MOD #   �'  @   a   INIT%JS+DOMAIN_MOD #   �'  @   a   INIT%JE+DOMAIN_MOD ;   *(  �   a   MULTI_DOMAIN_T%SUBDOMAINS+MULTI_DOMAIN_MOD :   �(  H   a   MULTI_DOMAIN_T%NUM_SUB_X+MULTI_DOMAIN_MOD :   ,)  H   a   MULTI_DOMAIN_T%NUM_SUB_Y+MULTI_DOMAIN_MOD D   t)  �   a   MULTI_DOMAIN_T%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 5    *  _   a   MULTI_DOMAIN_T%INIT+MULTI_DOMAIN_MOD 3   *  �       INIT_MULTI_DOMAIN+MULTI_DOMAIN_MOD 8   +  \   a   INIT_MULTI_DOMAIN%THIS+MULTI_DOMAIN_MOD A   w+  V   a   INIT_MULTI_DOMAIN%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD =   �+  @   a   INIT_MULTI_DOMAIN%NUM_SUB_X+MULTI_DOMAIN_MOD =   ,  @   a   INIT_MULTI_DOMAIN%NUM_SUB_Y+MULTI_DOMAIN_MOD G   M,  �   a   INIT_MULTI_DOMAIN%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 0   �,  �       MULTI_DOMAIN_T+MULTI_DOMAIN_MOD >   �-  ^   a   MULTI_DOMAIN_T%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD ;   .  �   a   MULTI_DOMAIN_T%SUBDOMAINS+MULTI_DOMAIN_MOD :   �.  H   a   MULTI_DOMAIN_T%NUM_SUB_X+MULTI_DOMAIN_MOD :   /  H   a   MULTI_DOMAIN_T%NUM_SUB_Y+MULTI_DOMAIN_MOD D   M/  �   a   MULTI_DOMAIN_T%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 5   �/  _   a   MULTI_DOMAIN_T%INIT+MULTI_DOMAIN_MOD !   X0  z       EXPLICIT_EULER_T .   �0  b   a   EXPLICIT_EULER_T%TIMESCHEME_T *   41  ]   a   EXPLICIT_EULER_T%TENDENCY &   �1  a   a   EXPLICIT_EULER_T%STEP $   �1  �       STEP_EXPLICIT_EULER )   t2  ^   a   STEP_EXPLICIT_EULER%THIS '   �2  U   a   STEP_EXPLICIT_EULER%V0 -   '3  X   a   STEP_EXPLICIT_EULER%OPERATOR 1   3  \   a   STEP_EXPLICIT_EULER%MULTI_DOMAIN '   �3  @   a   STEP_EXPLICIT_EULER%DT 