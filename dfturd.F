*Deck DFTURd
      Subroutine DFTURd(IOut,IPrint,NTT,UMat,UValBF,IV,V,MDV)
      Implicit Real*8(A-H,O-Z)
      Parameter (MaxTab=1000,LenLin=1000,MinEl=-2,MaxEl=118)
      Parameter(IAtomI=1,IAtomW=IAtomI+1,IShlTp=IAtomW+2,
     $  IUVal=IShlTp+1,MaxRlt=IUVal+1)
      Dimension ITable(MaxTab),IReslt(MaxRlt),Line(LenLin),UMat(NTT),
     $  UValBF(*),V(MDV),IEl(MinEl:MaxEl)
      Common /IO/ In
      Equivalence (FURlt,IURlt)
      Character*30 String
      Logical DEBUG,NumEl
#include "commonb.inc"
#include "commonmol.inc"
C
      DEBUG = .True.

      If(DEBUG) then
        Write(IOut,*)
        Write(IOut,*)
        Write(IOut,*)' ****************************** Hrant'
        Write(IOut,*)
        endIf

      If(DEBUG.or.IPrint.ge.2) Write(IOut,*) 'Entering DFTURd'
C
C     Initialize memory and mapping
C
      IMapAt = IV
      IMapL = IMapAt + NBasis 
      IMapBs = IMapL + NBasis
      IV = IMapBs + NAtoms+1
      Call TstCor(IV-1,MDV,'DFTURd: Out of memory at A.')
      Call MakMpC(.False.,.True.,1,0,I5DB1,I7FB1,NShell,ShellT,ShellC,
     $  JAn,MaxAtB,MaxL,V(IMapAt),V(IMapL))
      Call MakMpE(NAtoms,NBasis,V(IMapAt),V(IMapBs))
      If(DEBUG.or.IPrint.ge.2) then
        Write(Iout,*)' NBasis = ',NBasis
        Write(IOut,*)' NAtoms = ',NAtoms
        Call IMtOuS(IOut,'Basis to Atom Map (MapAt)',0,V(IMapAt),NBasis,
     $    1,NBasis,1)
        Call IMtOuS(IOut,'Angular Momentum to Atom Map (MapL)',0,
     $    V(IMapL),NBasis,1,NBasis,1)
        Call IMtOuS(IOut,'Atom to basis map (MapBs)',0,V(IMapBs),
     $    NAtoms+1,1,NAtoms+1,1)
        endIf
C
C     Set up rules for parsing input
C
      Call DFTUTb(IOut,MaxTab,ITable,IAtomI,IAtomW,IShlTp,IUVal)
C
C     Read lines from input file
C
      Call RdWEOF(1,In,Line,LenLin,IEof,Len)
   10 Call QPInit(ITable,0,0,0,2)
      Write(IOut,*) 'IEof:', IEof
      Call IClear(MaxRlt,IReslt)
C
C     Parse new line if successfully read in
C
      Call StrOut(6,Line,Len,1)
      IStat = QParse(IReslt,ITable,Line,Len)
      If(IStat.lt.0.or.IStat.gt.1) then
        Call QPErro(Line,Len,IOut)
      else if(IStat.eq.1) then
        Write(IOut,*) 'return'
        endIf

      If(DEBUG.or.IPrint.ge.2) then
        Write(IOut,*) "IReslts from QParse"
        Write(IOut,*) 'Col 1:'
        If(IReslt(IAtomI).gt.0) then
          Write(IOut,*) IReslt(IAtomI)
        Else
          Call StrOut(IOut,IReslt(IAtomW+1),IReslt(IAtomW),1)
          EndIf
        Write(IOut,*) 'Col 2:'
        Write(IOut,*) IReslt(IShlTp)
        Write(IOut,*) 'Col 3:'
        Call PrtFP(IReslt(IUVal),30,9,String)
        Write(IOut,*) String
        EndIf
C
C     Do conversions to store data into U array
C  
      Write(IOut,*) 'IReslt(IAtomW) is ', IReslt(IAtomW) 
C
      If(IReslt(IAtomW).gt.0) then
        Call FillEl(MinEl,MaxEl,IEl)
C        
      If(.not.NumEl(MaxEl,IEl,IReslt(IAtomW+1),Num)) 
     $  Call GauErr('Bad element label in input')
C
C     Building the Vector after adding U
C      Do I = 1,NAtoms
C        If(IAn(I).eq.Num) then
C          Write(IOut,*) 'I = ', I
C          Write(IOut,*) 
C          Do J = IIGet(V(IMapBs+I-1)),
C     $        IIGet(V(IMapBs+I))-1  
C            If(IIGet(V(IMapL+J-1)).eq.IReslt(IShlTp)) then
C              Do K = IIGet(V(IMapBs+IReslt(IAtomI)-1)), J
C                If(IIGet(V(IMapL+K-1)).eq.IReslt(IShlTp)) then
C                 IURlt = IReslt(IUVal)
C                 UValBF(j) = FURlt
C               endif
C              enddo
C            endif
C          enddo
C        endif 
C      enddo

c      Write(IOut,*) 'NAtoms:',NAtoms,'NBasis:',NNBasis
c      Do I = 1,NAtoms
c        Write(IOut,*) I,'IAn(I):', IAn(I), 'Num:', Num
c        If(IAn(I).eq.Num) then
c          Write(IOut,*) 'Identified Atom Number!' 
c          Write(IOut,*) 'JStart:',IIGet(V(IMapBs),I),'JEnd:',
c     $      IIGet(V(IMapBs),I+1)-1
c          Do J = IIGet(V(IMapBs),I),IIGet(V(IMapBs),I+1)-1
c           Write(IOut,*) J,'V(IMapAt):',IIGet(V(IMapAt),J),
c     $      'IReslt(IAn(I)):',IReslt(IAn(I)) 
c           Write(IOut,*) 'Identified Basis Number!' 
c            If(IIGet(V(IMapL),J).eq.IReslt(IShlTp)) then 
c              IURlt = IReslt(IUVal)
c              UValBF(j) = FURlt
c            EndIf
c          EndDo
c        EndIf 
c      EndDo 
      
        Do I = 1, NAtoms
         If(IAn(I).eq.Num) then
           Write(IOut,*) 'I = ', I
            Do J = InToWp(V(IMapBs+I-1)),
     $        InToWp(V(IMapBs+I))-1
             Write(IOut,*) 'J = ', J 
             If(InToWP(V(IMapL+J-1)).eq.IReslt(IShlTp)) then
               Write(IOut,*) 'J after if = ', J
                    IURlt = IReslt(IUVal)
                    UValBF(j) = FURlt
              EndIf
            EndDo
          EndIf
        EndDo
c
      ElseIf(IReslt(IAtomI).gt.0) then
c
c        Do 200 I = 1,NBasis
c          If(IIGet(V(IMapAt),i).eq.IReslt(IAtomI)) then 
c           If(IIGet(V(IMapL),i).eq.IReslt(IShlTp)) then
c            IURlt = IReslt(IUVal)
c             UValBF(i) = FURlt
c              endIf
c           endIf
c  200     Continue
c
CAAT+
c      IBsStr = IIGet(V(IMapBs),IReslt(IAtomI))
c      IBsEnd = IIGet(V(IMapBs),IReslt(IAtomI)+1)-1
c      Do I = IBsStr,IBsEnd
c        If(IIGet(V(IMapL),i).eq.IReslt(IShlTp)) then
c          IURlt = IReslt(IUVal)
c          UValBF(i) = FURlt
c        EndIf
c      EndDo
c
CAAT-
      Do  69 I = InToWp(V(IMapBs+IReslt(IAtomI)-1)),
     $      InToWp(V(IMapBs+IReslt(IAtomI)))-1
          If(InToWP(V(IMapL+I-1)).eq.IReslt(IShlTp)) then
            Do 70 J = InToWp(V(IMapBs+IReslt(IAtomI)-1)), I
              If(InToWP(V(IMapL+J-1)).eq.IReslt(IShlTp)) then
                IURlt = IReslt(IUVal) 
                UMat(I*(I-1)/2+J) = FURlt
              EndIf
  70        continue
          EndIf
  69    continue     
c
      Else
        Call GauErr('DFTURd cannot determine atom index')
        endIf
C
C     Read next line and loop back unless end of file           
C
      Call RdWEOF(1,In,Line,LenLin,IEof,Len)
      If(Len.gt.1.and.IEof.eq.0) Goto 10 

      If(DEBUG) then
        write(IOut,*) 'AAT---- NBasis..... : ', Nbasis
CAAT+
C        Call GauErr('Debugging in progress...AAT')
CAAT-
        Call OutMtS(IOut,'Final UValBF:',0,0,UValBF,NBasis,1,NBasis,1)
        Write(IOut,*)
        Write(IOut,*) ' UMat  .......:', UMat
        Write(IOut,*)' ****************************** ALI DID IT !! =) '
        Write(IOut,*)
    
       Write(IOut,*)
C        Call GauErr('STOP')
        endIf

      Return
      End
