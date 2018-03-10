*Deck DFTHub
      Subroutine DFTHub(IOut,IPrint,NBasis,NTT,IRwS,DA,DB,FA,FB,UMat,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
      logical DEBUG 
C
C     This routine forms the DFT+U correction to the Fock matrix using a
C     parameterized Hubbard U
C
C     Lee Thompson, 2015.
C
      Logical DEBUG
      Dimension V(MDV),DA(NTT),DB(NTT),FA(NTT),FB(NTT),UMat(NTT)
      Save Zero,Half,One,Two
      Data Zero/0.d0/,Half/0.5d0/,One/1.d0/,Two/2.d0/
C
      DEBUG = .True.
      If(DEBUG.or.IPrint.ge.2) Write(IOut,*) 'Entering DFTHub'

Chph+
      TrFPa1 = SCFTrc(FA,DA,NBasis,1)
      TrFPb1 = SCFTrc(FB,DB,NBasis,1)
      Write(IOut,*)' Hrant - At top of DFTHUB: <FP>alpha=',TrFPa1
      Write(IOut,*)' Hrant - At top of DFTHUB: <FP>beta =',TrFPb1
Chph-
C
      If(IPrint.ge.2) 
     $  Call LTOutS(IOut,'UMat in DFTHub',0,NBasis,UMat,0)
C
C     Compute Mulliken population matrix  
C
      ISMat = IV
      IFA = ISMat + NTT
      IFB = IFA + NTT
      IV = IFB + NTT

      Call FileIO(2,-IRwS,NTT,V(ISMat),0)

      ITmp = IV
      IMulA = ITmp + NBasis*NBasis
      IMulB = IMulA + NTT
      IV = IMulB + NTT  

      Call MultLL(NBasis,DA,V(ISMat),V(ITmp),Zero,One)
      Write(IOut,*)' AAT----- MultLL Alpha Before is :', DA
      Call XSyLT1(.False.,NBasis,NBasis,V(ITmp),V(IMulA))
      Call MultLL(NBasis,DB,V(ISMat),V(ITmp),Zero,One)
      Write(IOut,*)' AAT----- MultLL Beta Before is :', DB
      Call XSyLT1(.False.,NBasis,NBasis,V(ITmp),V(IMulB))

      If(DEBUG.or.IPrint.ge.2) then
        Call LTOutS(IOut,'DFTHub: SMat',0,NBasis,V(ISMat),0)
        Call LTOutS(IOut,'DFTHub: DA',0,NBasis,DA,0)
        Call LTOutS(IOut,'DFTHub: DB',0,NBasis,DB,0)
        Call LTOutS(IOut,'DFTHub: FA',0,NBasis,FA,0)
        Call LTOutS(IOut,'DFTHub: FB',0,NBasis,FB,0)
        Call LTOutS(IOut,'DFTHub: MulA',0,NBasis,V(IMulA),0)
        Call LTOutS(IOut,'DFTHub: MulB',0,NBasis,V(IMulB),0)
        endIf

C     Compute modification of the potential and add to Fock matrix 

      ICorrA  = IV
      ICorrB  = ICorrA  + NTT
      ICorrAA = ICorrB  + NTT
      ICorrBB = ICorrAA + NTT
      IV      = ICorrBB + NTT
C      
      K = 0
      Do I = 1, NBasis
        Do J = 1, I
          If (I.eq.J) then 
            RKD = One
          else
            RKD = Zero
            endIf
          V(ICorrA+K) = Half*UMat(K+1)*(RKD-(Two*V(IMulA+K)))
          V(ICorrB+K) = Half*UMat(K+1)*(RKD-(Two*V(IMulB+K)))
          K = K + 1
        endDo
      endDo
CAAT+
      Call LTOutS(IOut,'DFTHub: ICorrA before',0,NBasis,V(ICorrA),0)
      Call LTOutS(IOut,'DFTHub: ICorrB before',0,NBasis,V(ICorrB),0)
c
      Do I = 1, NTT
         V(ICorrAA+I-1) = V(ICorrA+I-1) * V(ISMat+I-1)
         V(ICorrBB+I-1) = V(ICorrB+I-1) * V(ISMat+I-1)
      enddo  
c      Call LTOutS(IOut,'DFTHub: CorrA',0,NBasis,V(ICorrA),0)
c      Call LTOutS(IOut,' SMat',0,NBasis,V(ISMat),0)
c      Call LTOutS(IOut,'DFTHub: CorrAA',0,NBasis,V(ICorrAA),0)
CAAT-      
c      Call MultLL(NBasis,V(ICorrA),V(ISMat),V(ITmp),Zero,One)
c      Call XSyLT1(.False.,NBasis,NBasis,V(ITmp),V(ICorrA))
c      Call MultLL(NBasis,V(ICorrB),V(ISMat),V(ITmp),Zero,One)
c      Call XSyLT1(.False.,NBasis,NBasis,V(ITmp),V(ICorrB))
CAAT+
c      Call LTOutS(IOut,'DFTHub: CorrA',0,NBasis,V(ICorrAA),0)
c      Call LTOutS(IOut,' SMat',0,NBasis,V(ISMat),0)
C     I edited the lines below where I changed ICorrA to ICorrAA
CAAT-      
      Do I = 1, NTT
          FA(I) = FA(I) + V(ICorrAA+I-1) 
          FB(I) = FB(I) + V(ICorrBB+I-1) 
      endDo

      If(DEBUG.or.IPrint.ge.2) then
        Call LTOutS(IOut,'DFTHub: ICorrA',0,NBasis,V(ICorrA),0)
        Call LTOutS(IOut,'DFTHub: ICorrAA',0,NBasis,V(ICorrAA),0)
        Call LTOutS(IOut,'DFTHub: ICorrB',0,NBasis,V(ICorrB),0)
        Call LTOutS(IOut,'DFTHub: ICorrBB',0,NBasis,V(ICorrBB),0)
        Call LTOutS(IOut,'DFTHub: IFA',0,NBasis,FA,0)
        Call LTOutS(IOut,'DFTHub: IFB',0,NBasis,FB,0)
        endIf

Chph+
      Write(IOut,*)' Hrant - At bot of DFTHUB: <UP>alpha=',
     $  SCFTrc(V(ICorrA),DA,NBasis,1)
      Write(IOut,*)' Hrant - At bot of DFTHUB: <UP>beta =',
     $  SCFTrc(V(ICorrB),DB,NBasis,1)
      TrFPa2 = SCFTrc(FA,DA,NBasis,1)
      TrFPb2 = SCFTrc(FB,DB,NBasis,1)
      Write(IOut,*)' Hrant - At bot of DFTHUB: <FP>alpha=',TrFPa2
      Write(IOut,*)' Hrant - At bot of DFTHUB: <FP>beta =',TrFPb2
      Write(IOut,*)'         diff of alpha = ',TrFPa1-TrFPa2
      Write(IOut,*)'         diff of beta  = ',Trfpb1-Trfpb2
Chph-


      Return
      End
