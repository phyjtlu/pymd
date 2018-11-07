      COMMON/eatint/iv(nmabig),jv(nmabig), kliend
      COMMON/eatd /eps(ntypes,ntypes),sig(ntypes,ntypes)
     &,RMAXLJ(ntypes,ntypes), RSLJ(ntypes,ntypes)
     &,RSPL(ntypes,ntypes),RSPLS(ntypes,ntypes)
     &,XMM(ntypes,ntypes),XMMS(ntypes,ntypes) 
     &,XM(ntypes,ntypes),C2(ntypes,ntypes)
     &,C3(ntypes,ntypes),pvdw,rrs(3) 
      COMMON/csurf/epss(ntypes,ntypes),sigs(ntypes,ntypes)
     &,surf ,tau(ntypes) 
      COMMON/ncsurf/ndir 
      COMMON/xlook/dellj,vlook(10000,ntypes,ntypes)
     &,dlook(10000,ntypes,ntypes) 

      save /eatint/,/eatd /,/csurf/,/ncsurf/,/xlook/