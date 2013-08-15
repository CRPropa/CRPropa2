PRO plot_skymap,FILENAME1, FILENAME2
  
  nlines = FILE_LINES(FILENAME1)
  lun1 = 122
  str = replicate({CR, crl:0.0, crb:0.0, crls:0.0, crbs:0.0}, nlines)

  nlinessource = FILE_LINES(FILENAME2)
  lun2 = 123
  source = replicate({SC, l:0.0, b:0.0}, nlinessource)

  openr, lun1, FILENAME1
  readf,lun1,str
  close,lun1
  free_lun,lun1

  openr, lun2, FILENAME2
  readf,lun2,source
  close,lun2
  free_lun,lun2

  vl = str.crl
  vb = str.crb
  vls = str.crl
  vbs = str.crb

  xcr = DBLARR(nlines)
  ycr = DBLARR(nlines)
  xs = DBLARR(nlines)
  ys = DBLARR(nlines)


  xscat = DBLARR(nlinessource)
  yscat = DBLARR(nlinessource)

  set_plot, 'x'

;  outfilename = name + '.eps'
;  set_plot,'ps'
;  DEVICE, /ENCAPSUL, FILE=outfilename, $ ;XSIZE=20, $
;          YSIZE=18, BITS_PER_PIXEL=4, /COLOR
;          /COLOR

;  !x.title = '!3Log!D10!N(E/eV)'
;  !y.title = '!3E!U2!NdN/dE (eV cm!U-2!N s!U-1!N sr!U-1!N)'

  red = FSC_COLOR('Red')
  magenta = FSC_COLOR('Magenta')
  green = FSC_COLOR('Forest Green')
  black = FSC_COLOR('Black')
  blue = FSC_COLOR('Blue')
  
;  aitoff, str.crvl, str.crvb, xcr,ycr     ;Convert to X,Y coordinates  ; Test all particles
  aitoff, source.l-180, source.b-90, xscat,yscat ;Convert to X,Y coordinates  
  aitoff, vl, vb, xcr,ycr     ;Convert to X,Y coordinates     ; Test just above GZK particles
  aitoff, vls, vbs, xs,ys ;Convert to X,Y coordinates  
;  aitoff, scat.l, scat.b, xscat,yscat     ;Convert to X,Y coordinates  

  aitoff_grid,/label,/new              ;Create labeled grid
  plots,xcr,ycr,color=red,psym=2       ;Overlay "star" positions  
  plots,xs,ys,color=green,psym=4
  plots,xscat,yscat,color=blue,psym=5


  legend,['!3UHECRs','!3Sources', '!3Other sources'],psym=[2,4,5],colors=[red,green,blue],/right,charsize=1.1,thick=1,charthick=1
;  DEVICE, /CLOSE
end
