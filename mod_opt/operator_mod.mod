  1,  u   k820309              2021.6.0    �~c                                                                                                          
       src/operators/operator_mod.f90 OPERATOR_MOD                                                     
       STVEC_T          @       �                                  
       MULTI_DOMAIN_T                   @                               '                      #COPY    #CREATE_SIMILAR 
   #UPDATE_S1V1    #UPDATE_S1V1S2V2    #UPDATE_S1V1V2    #UPDATE #   #ASSIGN_S1 $   #ASSIGN_S1V1 )   #ASSIGN_S1V1S2V2 /   #ASSIGN_S1V1V2 7   #ASSIGN >   1         �   � $                      �                        #COPY_STVEC    #         @                                                      #THIS    #FIN    #MULTI_DOMAIN              
                                                     #STVEC_T              
                                                     #STVEC_T              
                                      �             #MULTI_DOMAIN_T 	   1         �   � $                      �      
                  #CREATE_SIMILAR_STVEC    #         @                                                      #THIS    #DESTINATION              
                                                     #STVEC_T              
                                                     #STVEC_T    1         �   � $                      �                        #UPDATE_STVEC_S1V1    #         @                                                      #THIS    #SCALAR1    #V1    #MULTI_DOMAIN              
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      �             #MULTI_DOMAIN_T 	   1         �   � $                      �                        #UPDATE_STVEC_S1V1S2V2    #         @                                                      #THIS    #SCALAR1    #V1    #SCALAR2    #V2    #MULTI_DOMAIN              
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      
                
                                                     #STVEC_T              
                                      �             #MULTI_DOMAIN_T 	   1         �   � $                      �                        #UPDATE_STVEC_S1V1V2    #         @                                                      #THIS    #SCALAR1    #V1     #V2 !   #MULTI_DOMAIN "             
                                                     #STVEC_T              
                                      
                
                                                      #STVEC_T              
                                 !                    #STVEC_T              
                                 "     �             #MULTI_DOMAIN_T 	   4         �   � $                         @    #                    3         �   � $                         @             u #STVEC_T    #UPDATE_S1V1    #UPDATE_S1V1V2    #UPDATE_S1V1S2V2    1         �   � $                      �      $                  #ASSIGN_STVEC_S1 %   #         @                                  %                    #THIS &   #SCALAR1 '   #MULTI_DOMAIN (             
                                &                     #STVEC_T              
                                 '     
                
                                 (     �             #MULTI_DOMAIN_T 	   1         �   � $                      �      )                  #ASSIGN_STVEC_S1V1 *   #         @                                  *                    #THIS +   #SCALAR1 ,   #V1 -   #MULTI_DOMAIN .             
                                +                     #STVEC_T              
                                 ,     
                
                                 -                    #STVEC_T              
                                 .     �             #MULTI_DOMAIN_T 	   1         �   � $                      �      /             	     #ASSIGN_STVEC_S1V1S2V2 0   #         @                                  0                    #THIS 1   #SCALAR1 2   #V1 3   #SCALAR2 4   #V2 5   #MULTI_DOMAIN 6             
                                1                     #STVEC_T              
                                 2     
                
                                 3                    #STVEC_T              
                                 4     
                
                                 5                    #STVEC_T              
                                 6     �             #MULTI_DOMAIN_T 	   1         �   � $                      �      7             
 	    #ASSIGN_STVEC_S1V1V2 8   #         @                                  8                    #THIS 9   #SCALAR1 :   #V1 ;   #V2 <   #MULTI_DOMAIN =             
                                9                     #STVEC_T              
                                 :     
                
                                 ;                    #STVEC_T              
                                 <                    #STVEC_T              
                                 =     �             #MULTI_DOMAIN_T 	   4         �   � $                         @    >                    3         �   � $                         @             u #STVEC_T    #ASSIGN_S1V1 )   #ASSIGN_S1 $   #ASSIGN_S1V1S2V2 /   #ASSIGN_S1V1V2 7                     @               �           	     '�                   #GLOBAL_DOMAIN ?   #SUBDOMAINS Z   #NUM_SUB_X [   #NUM_SUB_Y \   #DEGREE_CONDENSATION ]   #INIT ^                �                               ?     �                      #DOMAIN_T @                     @              D           @     '�                    #XS A   #XE B   #YS C   #YE D   #DX E   #DY F   #IS G   #IE H   #JS I   #JE J   #NX K   #NY L   #X M   #Y N   #INIT O                �                              A                
                �                              B               
                �                              C               
                �                              D               
                �                              E                
                �                              F     (          
                �                              G     0                          �                              H     4                          �                              I     8       	                   �                              J     <       
                   �                              K     @                          �                              L     D                        �                              M            H                 
            &                                                      �                              N            �                 
            &                                           1         �   �                       �      O                  #INIT P   #         @                                  P                 	   #THIS Q   #XS R   #XE S   #IS T   #IE U   #YS V   #YE W   #JS X   #JE Y                                             Q     �               #DOMAIN_T @             
                                 R     
                
                                 S     
                
                                 T                     
                                 U                     
                                 V     
                
                                 W     
                
                                 X                     
                                 Y                      �                               Z            �       �             #DOMAIN_T @             &                   &                                                        �                              [     8                         �                              \     <                       �                              ]            @                            &                   &                                           1         �   � $                      �      ^                  #INIT_MULTI_DOMAIN _   #         @                                  _                    #THIS `   #GLOBAL_DOMAIN a   #NUM_SUB_X b   #NUM_SUB_Y c   #DEGREE_CONDENSATION d             
                                `     �              #MULTI_DOMAIN_T 	             
                                 a     �              #DOMAIN_T @             
                                 b                     
                                 c                   
                                 d                                 &                   &                                                             @               �           e     '�                   #GLOBAL_DOMAIN f   #SUBDOMAINS g   #NUM_SUB_X h   #NUM_SUB_Y i   #DEGREE_CONDENSATION j   #INIT k                �                               f     �                      #DOMAIN_T @              �                               g            �       �             #DOMAIN_T @             &                   &                                                        �                              h     8                         �                              i     <                       �                              j            @                            &                   &                                           1         �   � $                      �      k                  #INIT_MULTI_DOMAIN _                     @                          l     '                      #APPLY m   1         �   �                       �     m                  #APPLY_I n   #         @                                  n     	               #THIS o   #OUT p   #IN q   #MULTI_DOMAIN r             
                               o                     #OPERATOR_T l             
                               p                     #STVEC_T              
                               q                     #STVEC_T              
                                 r     �             #MULTI_DOMAIN_T e      �   4      fn#fn    �   H   J  STVEC_MOD !     O   J  MULTI_DOMAIN_MOD "   k        STVEC_T+STVEC_MOD '   r  X   a   STVEC_T%COPY+STVEC_MOD %   �  m       COPY_STVEC+STVEC_MOD *   7  U   a   COPY_STVEC%THIS+STVEC_MOD )   �  U   a   COPY_STVEC%FIN+STVEC_MOD 2   �  \   a   COPY_STVEC%MULTI_DOMAIN+STVEC_MOD 1   =  b   a   STVEC_T%CREATE_SIMILAR+STVEC_MOD /   �  c       CREATE_SIMILAR_STVEC+STVEC_MOD 4     U   a   CREATE_SIMILAR_STVEC%THIS+STVEC_MOD ;   W  U   a   CREATE_SIMILAR_STVEC%DESTINATION+STVEC_MOD .   �  _   a   STVEC_T%UPDATE_S1V1+STVEC_MOD ,     y       UPDATE_STVEC_S1V1+STVEC_MOD 1   �  U   a   UPDATE_STVEC_S1V1%THIS+STVEC_MOD 4   �  @   a   UPDATE_STVEC_S1V1%SCALAR1+STVEC_MOD /     U   a   UPDATE_STVEC_S1V1%V1+STVEC_MOD 9   n  \   a   UPDATE_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   �  c   a   STVEC_T%UPDATE_S1V1S2V2+STVEC_MOD 0   -  �       UPDATE_STVEC_S1V1S2V2+STVEC_MOD 5   �  U   a   UPDATE_STVEC_S1V1S2V2%THIS+STVEC_MOD 8   	  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   P	  U   a   UPDATE_STVEC_S1V1S2V2%V1+STVEC_MOD 8   �	  @   a   UPDATE_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3   �	  U   a   UPDATE_STVEC_S1V1S2V2%V2+STVEC_MOD =   :
  \   a   UPDATE_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   �
  a   a   STVEC_T%UPDATE_S1V1V2+STVEC_MOD .   �
  �       UPDATE_STVEC_S1V1V2+STVEC_MOD 3   x  U   a   UPDATE_STVEC_S1V1V2%THIS+STVEC_MOD 6   �  @   a   UPDATE_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1     U   a   UPDATE_STVEC_S1V1V2%V1+STVEC_MOD 1   b  U   a   UPDATE_STVEC_S1V1V2%V2+STVEC_MOD ;   �  \   a   UPDATE_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )     H   a   STVEC_T%UPDATE+STVEC_MOD %   [  �   `   gen@UPDATE+STVEC_MOD ,   �  ]   a   STVEC_T%ASSIGN_S1+STVEC_MOD *   >  q       ASSIGN_STVEC_S1+STVEC_MOD /   �  U   a   ASSIGN_STVEC_S1%THIS+STVEC_MOD 2     @   a   ASSIGN_STVEC_S1%SCALAR1+STVEC_MOD 7   D  \   a   ASSIGN_STVEC_S1%MULTI_DOMAIN+STVEC_MOD .   �  _   a   STVEC_T%ASSIGN_S1V1+STVEC_MOD ,   �  y       ASSIGN_STVEC_S1V1+STVEC_MOD 1   x  U   a   ASSIGN_STVEC_S1V1%THIS+STVEC_MOD 4   �  @   a   ASSIGN_STVEC_S1V1%SCALAR1+STVEC_MOD /     U   a   ASSIGN_STVEC_S1V1%V1+STVEC_MOD 9   b  \   a   ASSIGN_STVEC_S1V1%MULTI_DOMAIN+STVEC_MOD 2   �  c   a   STVEC_T%ASSIGN_S1V1S2V2+STVEC_MOD 0   !  �       ASSIGN_STVEC_S1V1S2V2+STVEC_MOD 5   �  U   a   ASSIGN_STVEC_S1V1S2V2%THIS+STVEC_MOD 8     @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR1+STVEC_MOD 3   D  U   a   ASSIGN_STVEC_S1V1S2V2%V1+STVEC_MOD 8   �  @   a   ASSIGN_STVEC_S1V1S2V2%SCALAR2+STVEC_MOD 3   �  U   a   ASSIGN_STVEC_S1V1S2V2%V2+STVEC_MOD =   .  \   a   ASSIGN_STVEC_S1V1S2V2%MULTI_DOMAIN+STVEC_MOD 0   �  a   a   STVEC_T%ASSIGN_S1V1V2+STVEC_MOD .   �  �       ASSIGN_STVEC_S1V1V2+STVEC_MOD 3   l  U   a   ASSIGN_STVEC_S1V1V2%THIS+STVEC_MOD 6   �  @   a   ASSIGN_STVEC_S1V1V2%SCALAR1+STVEC_MOD 1     U   a   ASSIGN_STVEC_S1V1V2%V1+STVEC_MOD 1   V  U   a   ASSIGN_STVEC_S1V1V2%V2+STVEC_MOD ;   �  \   a   ASSIGN_STVEC_S1V1V2%MULTI_DOMAIN+STVEC_MOD )     H   a   STVEC_T%ASSIGN+STVEC_MOD %   O  �   `   gen@ASSIGN+STVEC_MOD 0   �  �       MULTI_DOMAIN_T+MULTI_DOMAIN_MOD >   �  ^   a   MULTI_DOMAIN_T%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD $   �  �       DOMAIN_T+DOMAIN_MOD '   �  H   a   DOMAIN_T%XS+DOMAIN_MOD '     H   a   DOMAIN_T%XE+DOMAIN_MOD '   N  H   a   DOMAIN_T%YS+DOMAIN_MOD '   �  H   a   DOMAIN_T%YE+DOMAIN_MOD '   �  H   a   DOMAIN_T%DX+DOMAIN_MOD '   &  H   a   DOMAIN_T%DY+DOMAIN_MOD '   n  H   a   DOMAIN_T%IS+DOMAIN_MOD '   �  H   a   DOMAIN_T%IE+DOMAIN_MOD '   �  H   a   DOMAIN_T%JS+DOMAIN_MOD '   F  H   a   DOMAIN_T%JE+DOMAIN_MOD '   �  H   a   DOMAIN_T%NX+DOMAIN_MOD '   �  H   a   DOMAIN_T%NY+DOMAIN_MOD &     �   a   DOMAIN_T%X+DOMAIN_MOD &   �  �   a   DOMAIN_T%Y+DOMAIN_MOD )   F  R   a   DOMAIN_T%INIT+DOMAIN_MOD     �  �       INIT+DOMAIN_MOD %   *  V   a   INIT%THIS+DOMAIN_MOD #   �  @   a   INIT%XS+DOMAIN_MOD #   �  @   a   INIT%XE+DOMAIN_MOD #       @   a   INIT%IS+DOMAIN_MOD #   @   @   a   INIT%IE+DOMAIN_MOD #   �   @   a   INIT%YS+DOMAIN_MOD #   �   @   a   INIT%YE+DOMAIN_MOD #    !  @   a   INIT%JS+DOMAIN_MOD #   @!  @   a   INIT%JE+DOMAIN_MOD ;   �!  �   a   MULTI_DOMAIN_T%SUBDOMAINS+MULTI_DOMAIN_MOD :   :"  H   a   MULTI_DOMAIN_T%NUM_SUB_X+MULTI_DOMAIN_MOD :   �"  H   a   MULTI_DOMAIN_T%NUM_SUB_Y+MULTI_DOMAIN_MOD D   �"  �   a   MULTI_DOMAIN_T%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 5   v#  _   a   MULTI_DOMAIN_T%INIT+MULTI_DOMAIN_MOD 3   �#  �       INIT_MULTI_DOMAIN+MULTI_DOMAIN_MOD 8   q$  \   a   INIT_MULTI_DOMAIN%THIS+MULTI_DOMAIN_MOD A   �$  V   a   INIT_MULTI_DOMAIN%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD =   #%  @   a   INIT_MULTI_DOMAIN%NUM_SUB_X+MULTI_DOMAIN_MOD =   c%  @   a   INIT_MULTI_DOMAIN%NUM_SUB_Y+MULTI_DOMAIN_MOD G   �%  �   a   INIT_MULTI_DOMAIN%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 0   G&  �       MULTI_DOMAIN_T+MULTI_DOMAIN_MOD >   �&  ^   a   MULTI_DOMAIN_T%GLOBAL_DOMAIN+MULTI_DOMAIN_MOD ;   Y'  �   a   MULTI_DOMAIN_T%SUBDOMAINS+MULTI_DOMAIN_MOD :   (  H   a   MULTI_DOMAIN_T%NUM_SUB_X+MULTI_DOMAIN_MOD :   [(  H   a   MULTI_DOMAIN_T%NUM_SUB_Y+MULTI_DOMAIN_MOD D   �(  �   a   MULTI_DOMAIN_T%DEGREE_CONDENSATION+MULTI_DOMAIN_MOD 5   O)  _   a   MULTI_DOMAIN_T%INIT+MULTI_DOMAIN_MOD    �)  [       OPERATOR_T !   	*  U   a   OPERATOR_T%APPLY    ^*  u       APPLY_I    �*  X   a   APPLY_I%THIS    ++  U   a   APPLY_I%OUT    �+  U   a   APPLY_I%IN %   �+  \   a   APPLY_I%MULTI_DOMAIN 