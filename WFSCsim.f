      PROGRAM WFSCsim
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Simulate a timeline for WFS&C evolution to determine the optimal control  
C cadence.
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
      IMPLICIT REAL*8 (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCC
C 
C Define arrays used by the program     
C
      DIMENSION t(40000), rms(40000), rmslin(40000)     
      DIMENSION is(40000), ic(40000)     
C      
CCCCCCCCCCCCCCCCCCCC
C
C Common block related to the slope with which RMS wavefront error    
C deteriorates                                          
C
      COMMON /slopevals/ slopemin, slopemax, slopecur
C
C Common block related to (regular small) wavefront error jumps
C
      COMMON /jumpvals/ wjumpmax, probjump, wjumpcur 
C
C Common block related to (unusual big) wavefront error jumps
C      
      COMMON /jumpBvals/ wBjumpmin, wBjumpmax, probBjump, wBjumpcur
C      
CCCCCCCCCCCCCCCCCCCC
C
C Open output files
C
C A single simulated timeline with 4-day cadence:
      OPEN (UNIT=11,FILE='WFSCsim.out1',STATUS='UNKNOWN')
C Overall results as function of cadence:
      OPEN (UNIT=12,FILE='WFSCsim.out2',STATUS='UNKNOWN')
C Statistics with Monte-Carlo RMS of key quantities      
      OPEN (UNIT=13,FILE='WFSCsim.out3',STATUS='UNKNOWN')
C Analytical predictions      
      OPEN (UNIT=14,FILE='WFSCsim.out4',STATUS='UNKNOWN')
C Linear segments for the timeline in .out1 file       
      OPEN (UNIT=15,FILE='WFSCsim.out5',STATUS='UNKNOWN')
C Many simulated timelines for all cadences:      
      OPEN (UNIT=16,FILE='WFSCsim.out6',STATUS='UNKNOWN')
C
CCCCCCCCCCCCCCCCCCCC
C 
C Set constants
C      
C Number of days per year
C
      year = 365.24D0
C
CCCCCCCCCCCCCCCCCCCC
C      
C Set Simulation parameters
C
C Number of Monte-Carlo timelines to assess
C      
      NMC = 8000
C
C Initialize the random seed 
C
      WRITE (*,*) 'Initialize Random Seed IDUM (e.g., 68)'
      READ (*,*) idum
      WRITE (*,*)
      idum = (-1*ABS(idum))-1
C
C Set time of the simulations in days. The start is always at t=0.
C In this case we simulate for 10 years.
C      
      tend = 3652.0D0
C
C Set the step size for timeline
C      
      dt   = 0.1D0
C
C Calculate the number of timesteps to evaluate
C
      Nt   = INT(tend/dt)
C
C Set the Gaussian scatter in wavefront error in nm
C that is used for reporting results, but this is
C not used for any control decisions, which are always based on the underlying
C trends.
C      
      drms = 0.63D0
C      
CCCCCCCCCCCCCCCCCCCC
C
C Set properties of WFS&C operations
C
C Set the duration in hours taken up by a WFS visit, based on 53 minutes
C presented by Lallo at 09/20/2024 meeting, and 13 minute savings from cutting
C dithering + LOS      
C      
      durwfs = 40.0D0 / 60.0D0    
C      
C Set the duration in hours taken up by a WFS&C visit, based on approximate
C 2hr statement by Lallo at 09/20/2024 meeting, agreed by Comeau, and
C subtracting 2*13 minutes. 
C      
      durwfsc = 94.0D0 / 60.0D0 
C      
C Set the delay in days which with a wavefront control can be executed 
C after the last wavefront sensing. Per discussion with Comeau. 
C      
      dtc  = 2.0D0
C
C Set the probability with which an attempted WFS visit fails
C
      probfail = 0.01D0
C      
CCCCCCCCCCCCCCCCCCCC      
C
C Set the parameters that describe the guessed wavefront evolution 
C of JWST
C
C Set the minimum and maximum slope in nm/day with which RMS wavefront error
C deteriorates over time. This is set to av +/- RMS of the values reported
C by ../slope_analysis/slopefits_stats.com.      
C      
      slopemin = 0.0862222D0 - 0.0513105D0
      slopemax = 0.0862222D0 + 0.0513105D0
C
C Set the maximum size of sudden wavefront error jumps. Based on my
C analysis in jump_analysis/. Excludes big jumps added below that happen       
C over longer timescales. The probability distribution of jumps is set to
C be flat between 0 and 60.0D0.
C      
      wjumpmax = 60.0D0
C
C Set the average time between wavefront error jumps in days.
C      
      tjumpav  = year / 4.50D0
C
C Calculate the probability of a jump at any given timestep
C
      probjump = 1.0D0 / (tjumpav/dt)
C      
C Set the minimum and maximum size of big wavefront error jumps. Based
C on my analysis in jump_analysis/, with an assumed flat probability
C distribution in between. 
C      
      wBjumpmin = 100.0D0
      wBjumpmax = 500.0D0
C
C Set the average time between big wavefront error jumps in days.
C Set here at a rough guess of once every 1.63 years. This is based on the 
C calendar days between the 2 big jumps seen during science operations. 
C
      tBjumpav  = 595.0D0
C
C Calculate the probability of a bigjump at any given timestep
C
      probBjump = 1.0D0 / (tBjumpav/dt)
C
C Set the baseline wavefront error in nm that is achieved after a control.
C This is for the OTE only. Based on the most recent WFC activity (02/01/2025).
C      
      rmsbase = 62.83D0
C
CCCCCCCCCCCCCCCCCCCC
C
C Set parameters that define what wavefront error is acceptable
C      
C Set the baseline wavefront error that triggers a control at the next  
C opportunity, which is based per Comeau on a threshold of 8 nm above
C the baseline.
C    
      rmslim  = rmsbase + 8.0D0
C
C Set the maximum wavefront error that is acceptable for science observations,
C and otherwise some of them need to be repeated. This is based on
C analysis in repeat_analysis/ that shows that no observations were 
C ever repeated when the delta RMS was below 30 nm
C      
      rmssci  = (0.5D0*(rmsbase+rmslim)) + 30.0D0
C
C Set the fraction of the total time at wavefront error above rmssci that
C needs to be repeated. Based on analysis in jump_analysis/. Given the 
C shot noise, with only 4 visits ever repeated above 30nm, there isn 't a
C single rigorous answer here.
C
      fracsci = 0.0416D0     
C
C Set the maximum wavefront error that above which requirements are violated.
C This needs to exceed rmssci.
C      
      rmsreq  = 131.0D0
C
C Set the fraction of the total time at wavefront error above rmsreq that  
C needs to be repeated. 
C
      fracreq = 0.83D0
C
CCCCCCCCCCCCCCCCCCCCC 
C
C For reference, set the analytical approximation described in the Tech Memo.
C The total number of jumps per year equals
C     N_j = (year/tjumpav) + (year/tBjumpav) 
C The fraction of science observations that needs to be repeated is
C     f_r = N_r/N_j 
C where
C     N_r = 0.0     * (30/50) * (year/tjumpav) +      
C           fracsci * (20/50) * (year/tjumpav) +    
C           fracreq           * ((year/tBjumpav)          
C
      xNj = (year/tjumpav) + (year/tBjumpav)
      xNr = ((1.0D0-((rmssci-(0.5D0*(rmsbase+rmslim)))/wjumpmax))*
     &             fracsci*(year/tjumpav)) +
     &            (fracreq*(year/tBjumpav))    
      fractot = xNr/xNj
C
C Now evaluate the parameters defined in my Tech Memo
C      
      apar = (year*(durwfs/24.0D0))
      cpar = ((year*(slopemin+slopemax)/(2*(rmslim-rmsbase)))+xNj)*
     &         ((durwfsc-durwfs)/24.0D0)
      bpar = xNj*fractot/2.0D0
      dpar = (2.0D0*dtc*bpar)
C
      WRITE (*,*) 'Analytical Maximum', SQRT(apar/bpar)
      WRITE (*,*) 'fractot', fractot
      WRITE (*,*) ' '
C      
CCCCCCCCCCCCCCCCCCCC
C     
C Set the parameters that we wish to vary to find the most efficient way 
C to operate the observatory.
C      
C Set the wavefront sensing cadence in days, always starting after the
C last control. We start with the minimum value under consideration.
C
      its  = 2
C
C Set a marker that we can jump to to run the program with different values
C      
 201  dts = DBLE(its)
C      
CCCCCCCCCCCCCCCCCCCC
C
C Initialize Monte-Carlo sums
C
      Njump      = 0      
      NBjump     = 0      
      sum1fok1   = 0.0D0
      sum2fok1   = 0.0D0
      sum1fok2   = 0.0D0
      sum2fok2   = 0.0D0
      sum1fok3   = 0.0D0
      sum2fok3   = 0.0D0
      sum1fok4   = 0.0D0
      sum2fok4   = 0.0D0
      sum1wav    = 0.0D0  
      sum2wav    = 0.0D0
      sum1wrms   = 0.0D0  
      sum2wrms   = 0.0D0
      sum1durwf  = 0.0D0
      sum2durwf  = 0.0D0
      sum1dursci = 0.0D0
      sum2dursci = 0.0D0
      sum1durnet = 0.0D0
      sum2durnet = 0.0D0
C
C Loop over Monte-Carlo timelines
C
      DO iMC=1,NMC
C
C Initialize the RMS wavefront error at tstart of the simulation
C
        rmscur    = rmsbase
        rmslincur = rmsbase
C
C Write starting value to output file with linear slopes
C        
        IF ((iMC.EQ.1).AND.(its.EQ.4)) WRITE (15,*) 0.0D0, rmslincur
C        
C Initialize the time of the last wavefront sensing 
C
        tslast = 0.0D0
C        
C Initialize the time of the last wavefront control 
C
        tclast = 0.0D0 
C        
C Initialize the slope at the start of the simulation
C        
        CALL DRAWSLOPE(idum)
C
C Initialize the flag that determines whether a control needs to happen at
C the next opportunity       
C
        icontrol = 0
C
C Initialize the flag that determines whether a sensing is forced at 
C the next opportunity       
C
        isforce  = 0
C
C Initialize Monte-Carlo various sums:
C
C Number of timesteps at which the wavefront error is below the control        
C threshold (1); above, but below the the threshold acceptable for science (2);
C above, but below the requirements threshold (3); above the requirements
C threshold (4).
C        
        Nok1 = 0
        Nok2 = 0
        Nok3 = 0
        Nok4 = 0
C
C Number of times that a wavefront sensing is executed without a control 
C or with a control
C
        Nwfs  = 0
        Nwfsc = 0
C
C First and second moment sums of wavefront error
C        
        sumw1 = 0.0D0
        sumw2 = 0.0D0
C
CCCCCCCCCCCCCCCCCCCC
C        
C Loop over the timeline
C         
        DO j=1,Nt
C
C Set the new time
C           
          tcur = DBLE(j) * dt
C
C Allow for the possibility of a (regular small) wavefront error jump
C          
          CALL DRAWJUMP(idum)
C          
C Allow for the possibility of an (unusual big) wavefront error jump
C          
          CALL DRAWBIGJUMP(idum)
C          
C Set the new wavefront error and Gaussian scatter
C
          rmslincur = rmslincur + (slopecur*dt) 
          rmscur    = rmscur + (slopecur*dt) + wjumpcur + wBjumpcur
          drmscur   = drms * GASDEV(IDUM)
C
          IF (wjumpcur.GT.1.0D-10) THEN
             Njump    = Njump + 1
             ijumpnew = 1
C             WRITE (*,*) 'Small', wjumpcur
          END IF   
C          
          IF (wBjumpcur.GT.1.0D-10) THEN
             NBjump = NBjump + 1
             ijumpnew = 1
C             WRITE (*,*) 'Big', wBjumpcur
          END IF   
C          
C Defaults are that no wavefront sensing or control will happen, for which 
C flags are stored in arrays
C          
          is(j)  = 0
          ic(j)  = 0
C          
C If it is necessary and possible to do a wavefront control, then execute it,
C so that the wavefront error and control flag are reset. Remember in an 
C array that the control occurred. Also force a sensing to occur. If a control
C occurred, then also redraw the slope with which the wavefront error will
C increase going forward.
C
C Note: we assume that a planned WFC is always successful, even if the
C accompanying WFS need not be.          
C          
          IF ((icontrol.EQ.1).AND.((tcur-tslast).GE.dtc)) THEN
C Write the end of linear segments
            IF ((iMC.EQ.1).AND.(its.EQ.4)) WRITE (15,*) tcur, rmslincur
            rmscur    = rmsbase
            rmslincur = rmsbase
            drmscur  = 0.0D0 
            icontrol = 0
            ic(j)    = 1 
            isforce  = 1
            CALL DRAWSLOPE(idum)
C Write the beginning of linear segments
            IF ((iMC.EQ.1).AND.(its.EQ.4)) WRITE (15,*) tcur, rmslincur
          END IF   
C          
C Evaluate whether a sensing is going to be executed at this timestep.
C If it is forced, then reset the forcing flag. Also update the counters
C of how many sensings are performed with and without control
C          
          IF (((tcur-tslast).GE.dts).OR.(isforce.EQ.1)) THEN
C
            IF (isforce.EQ.1) THEN
              isforce = 0
              Nwfsc = Nwfsc + 1
            ELSE
              Nwfs  = Nwfs + 1
            END IF   
C
C If is(j)=1 it indicates that a sensing happened. If it is 2, then 
C the sensing reveals a jump that just happened.
C            
            tslast = tcur
            IF (ijumpnew.EQ.0) THEN
               is(j)  = 1
            ELSE
               is(j)  = 2
               ijumpnew = 0
            END IF
C            
C Based on a random drawing, decide whether the sensing visit failed or not.
C             
            ifail = 0
            IF (RAN1(idum).LE.probfail) ifail = 1 
C
          END IF   
C          
C If there was a WFS at this time, and it was successful,
C evaluate whether the wavefront error trend exceeds the control limit
C
          IF ((is(j).GE.1).AND.(ifail.EQ.0).AND.
     &        (rmscur.GE.rmslim)) THEN          
            icontrol = 1
          END IF   
C
C Store state after this timestep in arrays          
C
          t(j)      = tcur
          rms(j)    = rmscur + drmscur
          rmslin(j) = rmslincur
C
C Update sums
C
          sumw1 = sumw1 + rms(j)
          sumw2 = sumw2 + (rms(j)**2)
          IF (rms(j).LT.rmslim) THEN
            Nok1 = Nok1 + 1
          ELSE IF (rms(j).LT.rmssci) THEN     
            Nok2 = Nok2 + 1
          ELSE IF (rms(j).LT.rmsreq) THEN   
            Nok3 = Nok3 + 1
          ELSE
            Nok4 = Nok4 + 1            
          END IF
C            
C End the loop over the timeline           
C
        END DO
C
C Do a sanity check
C        
        IF (Nok1+Nok2+Nok3+Nok4.NE.Nt) STOP 'summing error'
C      
CCCCCCCCCCCCCCCCCCCC
C
C Write results of the simulated timeline only for the first
C simulation, and only for 4 day cadence. When a control is done, 
C always write the value both before and after control, at the same timestep.
C
        IF ((iMC.LE.1).AND.(its.EQ.4)) THEN
          DO j=1,Nt
            IF (ic(j).EQ.1) THEN     
              WRITE (11,'(2I5,2F12.5,2I5,F12.5)')
     &            its, iMC, t(j),
     &            MIN(rms(j-1),140.0D0), 1, 0, rmslin(j)
              WRITE (11,'(2I5,2F12.5,2I5,F12.5)')
     &            its, iMC, t(j),
     &            MIN(rms(j),140.0D0), is(j), ic(j), rmslin(j)
            ELSE
              WRITE (11,'(2I5,2F12.5,2I5,F12.5)')
     &            its, iMC, t(j),
     &            MIN(rms(j),140.0D0), is(j), ic(j), rmslin(j)
            END IF 
          END DO
        END IF 
C
C Write results of the simulated timelines for all cadences, possibly only for
C a subset of all Monte-carlos. Used to create WFE histograms outside of
C this program with histograms1.ipynb. Write results only on a 1-day cadence.
C A finer cadence yields larger output files, without meaningfully reducing
C the noise in the histograms that I produce from these files.         
C
        Ncad = INT(1.00001D0/dt)
C        
        IF (iMC.LE.NMC) THEN
          DO j=1,Nt
C              
            iwrite = 0              
            IF (Ncad*(j/Ncad).EQ.j) iwrite = 1
C            
            IF ((rms(j).GE.rmssci).AND.(rms(j).LT.rmsreq).AND.
     &          (RAN1(IDUM).LE.fracsci)) THEN
              iwrite = 0                 
            ELSE IF ((rms(j).GE.rmsreq).AND.
     &          (RAN1(IDUM).LE.fracreq)) THEN
              iwrite = 0                 
            END IF              
C
            IF (iwrite.EQ.1) THEN  
              WRITE (16,'(2I5,2F12.5)')
     &            its, iMC, t(j), MIN(rms(j),140.0D0)
            END IF 
C                
          END DO
        END IF 
C
C Evaluate statistics for this timeline
C Fraction of the time spent in various states
C        
        fok1 = DBLE(Nok1)/DBLE(Nt)
        fok2 = DBLE(Nok2)/DBLE(Nt)
        fok3 = DBLE(Nok3)/DBLE(Nt)
        fok4 = DBLE(Nok4)/DBLE(Nt)
C
C Overall average and RMS wavefront error
C        
        wav  = sumw1 / DBLE(Nt)
        wrms = SQRT((sumw2/DBLE(Nt))-(wav**2))
C
C Total duration in days of all WFS and WFSC activities, Express in
C units of days/year.                                       
C
        durwf = ( (DBLE(Nwfs)*(durwfs/24.0D0)) +
     &            (DBLE(Nwfsc)*(durwfsc/24.0D0)) ) * (year/tend)
C
C Net cost of this timeline, which also adds in the time spent on
C science observations that will need to be repeated. Express these in
C units of days/year.
C        
        dursci = ((fracsci*fok3) + (fracreq*fok4)) * year
        durnet = durwf + dursci
C        
C Write to the output file
C        
        WRITE (12,'(2I5,4F9.4,5F9.2)')
     &       its, iMC, fok1, fok2, fok3, fok4,
     &       wav, wrms, durwf, dursci, durnet
C
C Update Monte-Carlo sums
C
        sum1fok1   = sum1fok1   + fok1
        sum2fok1   = sum2fok1   + (fok1**2)
        sum1fok2   = sum1fok2   + fok2
        sum2fok2   = sum2fok2   + (fok2**2)
        sum1fok3   = sum1fok3   + fok3
        sum2fok3   = sum2fok3   + (fok3**2)
        sum1fok4   = sum1fok4   + fok4
        sum2fok4   = sum2fok4   + (fok4**2)
        sum1wav    = sum1wav    + wav  
        sum2wav    = sum2wav    + (wav**2)
        sum1wrms   = sum1wrms   + wrms  
        sum2wrms   = sum2wrms   + (wrms**2)
        sum1durwf  = sum1durwf  + durwf
        sum2durwf  = sum2durwf  + (durwf**2)
        sum1dursci = sum1dursci + dursci
        sum2dursci = sum2dursci + (dursci**2)
        sum1durnet = sum1durnet + durnet
        sum2durnet = sum2durnet + (durnet**2)
C      
CCCCCCCCCCCCCCCCCCCC
C
C End the Monte-carlo loop
C
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Evaluate Monte-carlo statistics
C
      avfok1    = sum1fok1/DBLE(NMC)  
      rmsfok1   = SQRT((sum2fok1/DBLE(NMC))-(avfok1**2))  
      avfok2    = sum1fok2/DBLE(NMC)  
      rmsfok2   = SQRT((sum2fok2/DBLE(NMC))-(avfok2**2)) 
      avfok3    = sum1fok3/DBLE(NMC)  
      rmsfok3   = SQRT((sum2fok3/DBLE(NMC))-(avfok3**2))  
      avfok4    = sum1fok4/DBLE(NMC)  
      rmsfok4   = SQRT((sum2fok4/DBLE(NMC))-(avfok4**2))  
      avwav     = sum1wav/DBLE(NMC)  
      rmswav    = SQRT((sum2wav/DBLE(NMC))-(avwav**2))  
      avwrms    = sum1wrms/DBLE(NMC) 
      rmswrms   = SQRT((sum2wrms/DBLE(NMC))-(avwrms**2))
      avdurwf   = sum1durwf/DBLE(NMC)
      rmsdurwf  = SQRT((sum2durwf/DBLE(NMC))-(avdurwf**2))
      avdursci  = sum1dursci/DBLE(NMC)
      rmsdursci = SQRT((sum2dursci/DBLE(NMC))-(avdursci**2))
      avdurnet  = sum1durnet/DBLE(NMC)
      rmsdurnet = SQRT((sum2durnet/DBLE(NMC))-(avdurnet**2))
C      
C Write to screen.
C The first line also writes the average number of jumps per year.
C      
      WRITE (*,'(I5,2F9.2)')
     &     its, (year/tend) * DBLE(Njump)/DBLE(NMC),
     &          (year/tend) * DBLE(NBjump)/DBLE(NMC)                    
      WRITE (*,'(A5,4F9.4,5F9.2)')
     &     'av', avfok1, avfok2, avfok3, avfok4, avwav, avwrms,
     &     avdurwf, avdursci, avdurnet
      WRITE (*,'(A5,4F9.4,5F9.2)')
     &     'rms', rmsfok1, rmsfok2, rmsfok3, rmsfok4, rmswrms, rmswrms,
     &     rmsdurwf, rmsdursci, rmsdurnet0
C
C Write to the output file
C       
      WRITE (13,'(I5,2(4F9.4,5F9.2))')
     &     its, avfok1, avfok2, avfok3, avfok4, avwav, avwrms,
     &          avdurwf, avdursci, avdurnet,
     &          rmsfok1, rmsfok2, rmsfok3, rmsfok4, rmswrms, rmswrms,
     &          rmsdurwf, rmsdursci, rmsdurnet
C
C Write analytical approximation
C      
      WRITE (14,'(I5,2(3F9.4,5F9.2))')
     &     its, (apar/dts) + cpar, (bpar*dts) + dpar,
     &          (apar/dts) + cpar + (bpar*dts) + dpar
C
CCCCCCCCCCCCCCCCCCCC
C
C Repeat with an increased WFS cadence      
C
      IF (its.LT.14) THEN
        its = its + 1
        GOTO 201 
      END IF   
C      
C Close output files
C
      CLOSE (UNIT=11)
      CLOSE (UNIT=12)
      CLOSE (UNIT=13)
      CLOSE (UNIT=14)
      CLOSE (UNIT=15)
      CLOSE (UNIT=16)
C      
C End program
C
      END


      SUBROUTINE DRAWSLOPE(idum)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Draw a random slope with which RMS wavefront error              
C deteriorates given parameters in the common block /slopevals/.
C Return the slope in that same common block as slopecur.      
C IDUM is the random seed. 
C
C Present implementation uses a Gaussian deviate.      
C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      
      IMPLICIT REAL*8 (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCC
C 
      COMMON /slopevals/ slopemin, slopemax, slopecur
C
CCCCCCCCCCCCCCCCCCCC
C
C      slopecur = slopemin + (RAN1(idum)*(slopemax-slopemin))
C
C Draw random Gaussian deviate, but reject outside of [-1,1]
C      
 701  scur = GASDEV(idum)
      IF (ABS(scur).GT.1.0D0) GOTO 701
C
      slopecur = (0.5D0*(slopemax+slopemin)) +
     &           (scur*0.5D0*(slopemax-slopemin))
C     WRITE (*,'(F12.3)') slopecur 
C       
CCCCCCCCCCCCCCCCCCCC
C
      END 


      SUBROUTINE DRAWJUMP(idum)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Draw a random RMS wavefront error jump (usually zero)             
C given parameters in the common block /jumpvals/.
C IDUM is the random seed. 
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      
      IMPLICIT REAL*8 (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCC
C 
      COMMON /jumpvals/ wjumpmax, probjump, wjumpcur 
C
CCCCCCCCCCCCCCCCCCCC
C
      IF (RAN1(idum).LE.probjump) THEN
        wjumpcur = wjumpmax * RAN1(idum)
      ELSE         
        wjumpcur = 0.0D0
      END IF
C
CCCCCCCCCCCCCCCCCCCC
C
      END       


      SUBROUTINE DRAWBIGJUMP(idum)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Draw a random big RMS wavefront error jump (usually zero)             
C given parameters in the common block /jumpBvals/.
C Return the slope in that same common block as slopecur.      
C IDUM is the random seed. 
C      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C      
      IMPLICIT REAL*8 (a-h,o-z)
C
CCCCCCCCCCCCCCCCCCCC
C 
      COMMON /jumpBvals/ wBjumpmin, wBjumpmax, probBjump, wBjumpcur 
C
CCCCCCCCCCCCCCCCCCCC
C
      IF (RAN1(idum).LE.probBjump) THEN
        wBjumpcur = wBjumpmin + ((wBjumpmax-wBjumpmin)*RAN1(idum))
      ELSE         
        wBjumpcur = 0.0D0
      END IF
C
CCCCCCCCCCCCCCCCCCCC
C
      END       


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Numerical Recipes Routines
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,
     *AM=1.0D0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2D-7,RNMX=1.0D0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      REAL*8 FUNCTION gasdev(idum)
      INTEGER idum
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.0D0*ran1(idum)-1.0D0
        v2=2.0D0*ran1(idum)-1.0D0
        rsq=v1**2+v2**2
        if(rsq.ge.1.0D0.or.rsq.eq.0.0D0)goto 1
        fac=sqrt(-2.0D0*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
