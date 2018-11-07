                          function smooth(x)
C   --------------------------------------------------------------------------
C      cut-off function to smoothly remove the long range interactions

       implicit real*8(a-h,o-z)

      include 'common_files.inc'
        if(x.lt.2.60d0) then
             if(x.lt.2.45d0) then
                  smooth = ((ro/x)**2)*exp(2.0*(aa2 - (x/rctb)**xnc))
             else
                  xx = x - 2.45d0
                  smooth = ttb0 + xx*(ttb1 + xx*(ttb2 + xx*ttb3))
             endif
        else
             smooth = 0.0d0
        endif

       return
       end


