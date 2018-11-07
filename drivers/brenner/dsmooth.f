                       function Dsmooth(x)
C   --------------------------------------------------------------------------
C      derivative of the cut-off function

       implicit real*8(a-h,o-z)
      include 'common_files.inc'

       if(x.lt.2.60d0) then
             if(x.lt.2.45d0) then
                  Dsmooth = -((ro/x)**2)*2.0*xnc*((x/rctb)**xnc)/x
     .            *exp(2.0*(aa2 - (x/rctb)**xnc))
             else
                  xx = x - 2.45d0
                  Dsmooth = ttb1 + xx*(2.0d0*ttb2 + xx*3.0d0*ttb3)
             endif
       else
             Dsmooth = 0.0d0
       endif
c       Dsmooth = -Dsmooth
       return
       end
