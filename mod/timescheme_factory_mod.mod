  �8  �   k820309              2021.6.0    �lc                                                                                                          
       src/timeschemes/timesheme_factory_mod.f90 TIMESCHEME_FACTORY_MOD                                                     
       STVEC_T                                                     
       TIMESCHEME_T                                                     
       EXPLICIT_EULER_T                      @                              
       RK4_T                   @                               '                      #COPY    #CREATE_SIMILAR    #UPDATE_S1V1    #UPDATE_S1V1S2V2    #UPDATE_S1V1V2    #UPDATE %   #ASSIGN_S1 &   #ASSIGN_S1V1 +   #ASSIGN_S1V1S2V2 1   #ASSIGN_S1V1V2 9   #ASSIGN @   1         �   � $                      �                        #COPY_STVEC    #         @                                                      #THIS    #FIN 	   #MULTI_DOMAIN 
             
                                                     #STVEC_T              
                                 	                    #STVEC_T              
                                 
     �             #MULTI_DOMAIN_T    1         �   � $                     �                        #CREATE_SIMILAR_STVEC    #         @                                                     #THIS    #DESTINATION              
                                                     #STVEC_T              
                                                     #STVEC_T    1         �   � $                      �                        #UPDATE_STVEC_S1V1    #         @                                                      #THIS    #SCALAR1    #V1    #MULTI_DOMAIN              
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
                        @               @         J     '�                    #TIMESCHEME_T K   #TENDENCY L   #STEP M                �                               K                            #TIMESCHEME_T A               �                              L                            #STVEC_T    1         �   � $                      �     M                  #STEP_EXPLICIT_EULER N   #         @                                  N                    #THIS O   #V0 P   #OPERATOR Q   #MULTI_DOMAIN R   #DT S             
                                O     �               #EXPLICIT_EULER_T J             
                                P                     #STVEC_T              
                                Q                     #OPERATOR_T G             
                                  R     �             #MULTI_DOMAIN_T              
                                 S     
                        @               @         T     '�                   #TIMESCHEME_T U   #K1 V   #K2 W   #K3 X   #K4 Y   #Y Z   #STEP [                � $                              U                            #TIMESCHEME_T A               � $                             V                            #STVEC_T                � $                             W             �              #STVEC_T                � $                             X                           #STVEC_T                � $                             Y             �             #STVEC_T                � $                             Z                           #STVEC_T    1         �   � $                      �     [                  #STEP_RK4 \   #         @     @                            \                    #THIS ]   #V0 ^   #OPERATOR _   #MULTI_DOMAIN `   #DT a             
                                ]     �              #RK4_T T             
                                ^                     #STVEC_T              
                                _                     #OPERATOR_T G             
                                  `     �             #MULTI_DOMAIN_T              
                                 a     
                        @               �                '�                   #GLOBAL_DOMAIN b   #SUBDOMAINS }   #NUM_SUB_X ~   #NUM_SUB_Y    #DEGREE_CONDENSATION �   #INIT �                �                               b     �                      #DOMAIN_T c                     @              D           c     '�                    #XS d   #XE e   #YS f   #YE g   #DX h   #DY i   #IS j   #IE k   #JS l   #JE m   #NX n   #NY o   #X p   #Y q   #INIT r                �                              d                
                �                              e               
                �                              f               
                �                              g               
                �                              h                
                �                              i     (          
                �                              j     0                          �                              k     4                          �                              l     8       	                   �                              m     <       
                   �                              n     @                          �                              o     D                        �                              p            H                 
            &                                                      �                              q            �                 
            &                                           1         �   �                       �      r                  #INIT s   #         @                                  s                 	   #THIS t   #XS u   #XE v   #IS w   #IE x   #YS y   #YE z   #JS {   #JE |                                             t     �               #DOMAIN_T c             
                                 u     
                
                                 v     
                
                                 w                     
                                 x                     
                                 y     
                
                                 z     
                
                                 {                     
                                 |                      �                               }            �       �             #DOMAIN_T c             &                   &                                                        �                              ~     8                         �                                   <                       �                              �            @                            &                   &                                           1         �   � $                      �      �                  #INIT_MULTI_DOMAIN �   #         @                                  �                    #THIS �   #GLOBAL_DOMAIN �   #NUM_SUB_X �   #NUM_SUB_Y �   #DEGREE_CONDENSATION �             
                                �     �              #MULTI_DOMAIN_T              
                                 �     �              #DOMAIN_T c             
                                 �                     
                                 �                   
                                 �                                 &                   &                                                             @                         G     '                      #APPLY �   1         �   �                       �     �                  #APPLY_I �   #         @                                 �     	               #THIS �   #OUT �   #IN �   #MULTI_DOMAIN �             
                               �                     #OPERATOR_T G             
                               �                     #STVEC_T              
                               �                     #STVEC_T              
                                 �     �             #MULTI_DOMAIN_T    #         @                                   �                    #TIMESCHEME �   #V �   #TIMESCHEME_NAME �            D @                              �                     #TIMESCHEME_T A             
  @                              �                    #STVEC_T              
                                �                    1 #         @                                  �                    #TIMESCHEME �   #V �            D @                              �                     #TIMESCHEME_T A             
                                 �                    #STVEC_T    #         @                                  �                    #TIMESCHEME �   #V �            D @                              �                     #TIMESCHEME_T A             
                                 �                    #STVEC_T       �   I      fn#fn    �   H   J  STVEC_MOD    1  M   J  TIMESCHEME_MOD #   ~  Q   J  EXPLICIT_EULER_MOD    �  F   J  RK4_MOD "           STVEC_T+STVEC_MOD '     X   a   STVEC_T%COPY+STVEC_MOD %   t  m       COPY_STVEC+STVEC_MOD *   �  U   a   COPY_STVEC%THIS+STVEC_MOD )   6  U   a   COPY_STVEC%FIN+STVEC_MOD 2   �  \   a   COPY_STVEC%MULTI_DOMAIN+STVEC_MOD 1   �  b   a   STVEC_T%CREATE_SIMILAR+STVEC_MOD /   I  c       CREATE_SIMILAR_STVEC+STVEC_MOD 4   �  U   a   CREATE_SIMILAR_STVEC%THIS+STVEC_MOD ;     U   a   CREATE_SIMILAR_STVEC%DESTINATION+STVEC_MOD .   V  _   a   STVEC_T%UPDATE_S1V1+STVEC_MOD ,   �  y       UPDATE_STVEC_S1V1+STVEC_MOD 1   .  U   a   UPDATE_STVEC_S1V1%THIS+STVEC_MOD 4   �  @   a   UPDATE_STVEC_S1V1%SCALAR1+STVEC_MOD /   �  U   a   UPDATE_STVEC_S1V1%V1+STVEC_MOD 9     \   a   UPDATE_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   t  c   a   STVEC_T%UPDATE_S1V1S2V2+STVEC_MOD 0   �  �       UPDATE_STVEC_S1V1S2V2+STVEC_MOD 5   e	  U   a   UPDATE_STVEC_S1V1S2V2%THIS+STVEC_MOD 8   �	  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   �	  U   a   UPDATE_STVEC_S1V1S2V2%V1+STVEC_MOD 8   O
  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3   �
  U   a   UPDATE_STVEC_S1V1S2V2%V2+STVEC_MOD =   �
  \   a   UPDATE_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   @  a   a   STVEC_T%UPDATE_S1V1V2+STVEC_MOD .   �  �       UPDATE_STVEC_S1V1V2+STVEC_MOD 3   "  U   a   UPDATE_STVEC_S1V1V2%THIS+STVEC_MOD 6   w  @   a   UPDATE_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1   �  U   a   UPDATE_STVEC_S1V1V2%V1+STVEC_MOD 1     U   a   UPDATE_STVEC_S1V1V2%V2+STVEC_MOD ;   a  \   a   UPDATE_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )   �  H   a   STVEC_T%UPDATE+STVEC_MOD %     �   `   gen@UPDATE+STVEC_MOD ,   �  ]   a   STVEC_T%ASSIGN_S1+STVEC_MOD *   �  q       ASSIGN_STVEC_S1+STVEC_MOD /   Y  U   a   ASSIGN_STVEC_S1%THIS+STVEC_MOD 2   �  @   a   ASSIGN_STVEC_S1%SCALAR1+STVEC_MOD 7   �  \   a   ASSIGN_STVEC_S1%MULTI_DOMAIN+STVEC_MOD .   J  _   a   STVEC_T%ASSIGN_S1V1+STVEC_MOD ,   �  y       ASSIGN_STVEC_S1V1+STVEC_MOD 1   "  U   a   ASSIGN_STVEC_S1V1%THIS+STVEC_MOD 4   w  @   a   ASSIGN_STVEC_S1V1%SCALAR1+STVEC_MOD /   �  U   a   ASSIGN_STVEC_S1V1%V1+STVEC_MOD 9     \   a   ASSIGN_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   h  c   a   STVEC_T%ASSIGN_S1V1S2V2+STVEC_MOD 0   �  �       ASSIGN_STVEC_S1V1S2V2+STVEC_MOD 5   Y  U   a   ASSIGN_STVEC_S1V1S2V2%THIS+STVEC_MOD 8   �  @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   �  U   a   ASSIGN_STVEC_S1V1S2V2%V1+STVEC_MOD 8   C  @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3   �  U   a   ASSIGN_STVEC_S1V1S2V2%V2+STVEC_MOD =   �  \   a   ASSIGN_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   4  a   a   STVEC_T%ASSIGN_S1V1V2+STVEC_MOD .   �  �       ASSIGN_STVEC_S1V1V2+STVEC_MOD 3     U   a   ASSIGN_STVEC_S1V1V2%THIS+STVEC_MOD 6   k  @   a   ASSIGN_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1   �  U   a   ASSIGN_STVEC_S1V1V2%V1+STVEC_MOD 1      U   a   ASSIGN_STVEC_S1V1V2%V2+STVEC_MOD ;   U  \   a   ASSIGN_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )   �  H   a   STVEC_T%ASSIGN+STVEC_MOD %   �  �   `   gen@ASSIGN+STVEC_MOD ,   �  Z       TIMESCHEME_T+TIMESCHEME_MOD 1   �  R   a   TIMESCHEME_T%STEP+TIMESCHEME_MOD $   :  �       STEP+TIMESCHEME_MOD )   �  Z   a   STEP%THIS+TIMESCHEME_MOD '     U   a   STEP%V0+TIMESCHEME_MOD -   k  X   a   STEP%OPERATOR+TIMESCHEME_MOD 1   �  \   a   STEP%MULTI_DOMAIN+TIMESCHEME_MOD '     @   a   STEP%DT+TIMESCHEME_MOD 4   _  z       EXPLICIT_EULER_T+EXPLICIT_EULER_MOD A   �  b   a   EXPLICIT_EULER_T%TIMESCHEME_T+EXPLICIT_EULER_MOD =   ;  ]   a   EXPLICIT_EULER_T%TENDENCY+EXPLICIT_EULER_MOD 9   �  a   a   EXPLICIT_EULER_T%STEP+EXPLICIT_EULER_MOD 7   �  �       STEP_EXPLICIT_EULER+EXPLICIT_EULER_MOD <   {  ^   a   STEP_EXPLICIT_EULER%THIS+EXPLICIT_EULER_MOD :   �  U   a   STEP_EXPLICIT_EULER%V0+EXPLICIT_EULER_MOD @   .  X   a   STEP_EXPLICIT_EULER%OPERATOR+EXPLICIT_EULER_MOD D   �  \   a   STEP_EXPLICIT_EULER%MULTI_DOMAIN+EXPLICIT_EULER_MOD :   �  @   a   STEP_EXPLICIT_EULER%DT+EXPLICIT_EULER_MOD    "  �       RK4_T+RK4_MOD +   �  b   a   RK4_T%TIMESCHEME_T+RK4_MOD !      ]   a   RK4_T%K1+RK4_MOD !   t   ]   a   RK4_T%K2+RK4_MOD !   �   ]   a   RK4_T%K3+RK4_MOD !   .!  ]   a   RK4_T%K4+RK4_MOD     �!  ]   a   RK4_T%Y+RK4_MOD #   �!  V   a   RK4_T%STEP+RK4_MOD !   >"  �      STEP_RK4+RK4_MOD &   �"  S   a   STEP_RK4%THIS+RK4_MOD $   #  U   a   STEP_RK4%V0+RK4_MOD *   h#  X   a   STEP_RK4%OPERATOR+RK4_MOD .   �#  \   a   STEP_RK4%MULTI_DOMAIN+RK4_MOD $   $  @   a   STEP_RK4%DT+RK4_MOD 0   \$  �       MULTI_DOMAIN_T+MULTI_DOMAIN_MOD >   %  ^   a   MULTI_DOMAIN_T%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD $   n%  �       DOMAIN_T+DOMAIN_MOD '   6&  H   a   DOMAIN_T%XS+DOMAIN_MOD '   ~&  H   a   DOMAIN_T%XE+DOMAIN_MOD '   �&  H   a   DOMAIN_T%YS+DOMAIN_MOD '   '  H   a   DOMAIN_T%YE+DOMAIN_MOD '   V'  H   a   DOMAIN_T%DX+DOMAIN_MOD '   �'  H   a   DOMAIN_T%DY+DOMAIN_MOD '   �'  H   a   DOMAIN_T%IS+DOMAIN_MOD '   .(  H   a   DOMAIN_T%IE+DOMAIN_MOD '   v(  H   a   DOMAIN_T%JS+DOMAIN_MOD '   �(  H   a   DOMAIN_T%JE+DOMAIN_MOD '   )  H   a   DOMAIN_T%NX+DOMAIN_MOD '   N)  H   a   DOMAIN_T%NY+DOMAIN_MOD &   �)  �   a   DOMAIN_T%X+DOMAIN_MOD &   **  �   a   DOMAIN_T%Y+DOMAIN_MOD )   �*  R   a   DOMAIN_T%INIT+DOMAIN_MOD     +  �       INIT+DOMAIN_MOD %   �+  V   a   INIT%THIS+DOMAIN_MOD #   �+  @   a   INIT%XS+DOMAIN_MOD #   8,  @   a   INIT%XE+DOMAIN_MOD #   x,  @   a   INIT%IS+DOMAIN_MOD #   �,  @   a   INIT%IE+DOMAIN_MOD #   �,  @   a   INIT%YS+DOMAIN_MOD #   8-  @   a   INIT%YE+DOMAIN_MOD #   x-  @   a   INIT%JS+DOMAIN_MOD #   �-  @   a   INIT%JE+DOMAIN_MOD ;   �-  �   a   MULTI_DOMAIN_T%SUBDOMAINS+MULTI_DOMAIN_MOD :   �.  H   a   MULTI_DOMAIN_T%NUM_SUB_X+MULTI_DOMAIN_MOD :   �.  H   a   MULTI_DOMAIN_T%NUM_SUB_Y+MULTI_DOMAIN_MOD D   B/  �   a   MULTI_DOMAIN_T%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 5   �/  _   a   MULTI_DOMAIN_T%INIT+MULTI_DOMAIN_MOD 3   M0  �       INIT_MULTI_DOMAIN+MULTI_DOMAIN_MOD 8   �0  \   a   INIT_MULTI_DOMAIN%THIS+MULTI_DOMAIN_MOD A   E1  V   a   INIT_MULTI_DOMAIN%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD =   �1  @   a   INIT_MULTI_DOMAIN%NUM_SUB_X+MULTI_DOMAIN_MOD =   �1  @   a   INIT_MULTI_DOMAIN%NUM_SUB_Y+MULTI_DOMAIN_MOD G   2  �   a   INIT_MULTI_DOMAIN%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD (   �2  [       OPERATOR_T+OPERATOR_MOD .   3  U   a   OPERATOR_T%APPLY+OPERATOR_MOD %   o3  u       APPLY_I+OPERATOR_MOD *   �3  X   a   APPLY_I%THIS+OPERATOR_MOD )   <4  U   a   APPLY_I%OUT+OPERATOR_MOD (   �4  U   a   APPLY_I%IN+OPERATOR_MOD 2   �4  \   a   APPLY_I%MULTI_DOMAIN+OPERATOR_MOD "   B5  t       CREATE_TIMESCHEME -   �5  Z   a   CREATE_TIMESCHEME%TIMESCHEME $   6  U   a   CREATE_TIMESCHEME%V 2   e6  L   a   CREATE_TIMESCHEME%TIMESCHEME_NAME 1   �6  _       CREATE_EXPLICIT_EULER_TIMESCHEME <   7  Z   a   CREATE_EXPLICIT_EULER_TIMESCHEME%TIMESCHEME 3   j7  U   a   CREATE_EXPLICIT_EULER_TIMESCHEME%V &   �7  _       CREATE_RK4_TIMESCHEME 1   8  Z   a   CREATE_RK4_TIMESCHEME%TIMESCHEME (   x8  U   a   CREATE_RK4_TIMESCHEME%V 