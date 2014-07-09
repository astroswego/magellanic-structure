;##############################################################################
;
; Copyright (C) 2012-2013, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/software
;
; If you have found this software useful for your research, I would
; appreciate an acknowledgment to the use of the "LTS_PLANEFIT program
; described in Cappellari et al. (2013, MNRAS, 432, 1709), which
; combines the Least Trimmed Squares robust technique of Rousseeuw &
; Driessen (2006) into a least-squares fitting algorithm which allows
; for errors in all variables variables and intrinsic scatter." 
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;##############################################################################
;+
; NAME:
;       LTS_PLANEFIT
; PURPOSE:
;       Best plane *robust* fit to data with errors in all three 
;       coordinates and fitting for the intrinsic scatter.
; EXPLANATION:
;       Linear Least-squares approximation in two-dimension (z = a + b*x + c*y),
;       when x, y and z data have errors, using the Least Trimmed Squares
;       regression technique for robustness.
;       See Sec.3.2 here http://adsabs.harvard.edu/abs/2013MNRAS.432.1709C
;
;       Uses the robust LTS-FAST technique of Rousseeuw & Driessen (2006)
;       http://www.springerlink.com/content/06k45m57x01028x6/
;       See also http://books.google.co.uk/books?id=woaH_73s-MwC&pg=PA15
;
; CALLING EXAMPLE:
;       lts_planefit, x, y, z, sigx, sigy, sigz, par, sig_par, chi_sq, $
;           FRAC=frac, PLOT=plot, GOOD=good, BAD=bad, RMS=rms, $
;           NRAND=nrand, /BAYES, OVERPLOT=overplot, PIVOTX=pivotx, PIVOTY=pivoty, $
;           TEXT=text, EPSZ=epsz
;
; MODIFICATION HISTORY:
;       V1.0: Michele Cappellari, Oxford, 21 March 2011
;       V2.0: Converted from lts_linefit. MC, Oxford, 06 April 2011
;       V2.01: Added PIVOT keyword, MC, Oxford, 1 August 2011
;       V2.02: Fixed program stop affecting earlier IDL versions.
;           Thanks to Xue-Guang Zhang for reporting the problem
;           and the solution. MC, Turku, 10 July 2013
;       V2.03: Scale line spacing with character size in text output.
;           MC, Oxford, 19 September 2013
;       V2.04: Check that all input vectors have the same size.
;           MC, Baltimore, 8 June 2014
;       V2.05: Text plotting changes. MC, Oxford, 26 June 2014
;-
;------------------------------------------------------------------------------
function ltsp_residuals, abc, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz
compile_opt idl2, hidden

res = (abc[0] + abc[1]*x + abc[2]*y - z) / sqrt((abc[1]*sigx)^2 + (abc[2]*sigy)^2 + sigz^2)

return, res
end
;----------------------------------------------------------------------------
pro ltsp_fitting, x, y, z, sigx, sigy, sigz, abc, sig_abc, chi2
compile_opt idl2, hidden

fcnargs = {X:x, Y:y, Z:z, SIGX:sigx, SIGY:sigy, SIGZ:sigz}
abc = MPFIT('ltsp_residuals', abc, FUNCTARGS=fcnargs, $
        /QUIET, BESTNORM=chi2, PERROR=sig_abc)

end
;----------------------------------------------------------------------------
pro ltsp_algorithm, x, y, z, sigx, sigy, sigz, h, abc, good
compile_opt idl2, hidden

; Robust least trimmed squares regression.
; Pg. 38 of Rousseeuw & Driessen (2006)
; http://www.springerlink.com/content/06k45m57x01028x6/
;
n = n_elements(x)
m = 500 ; Number of random starting points
abcV = dblarr(m,3)
chi2V = dblarr(m)
for j=0,m-1 do begin ; Draw m random starting points
    w = (sort(randomu(s,n)))[0:2] ; Reshuffle the list without repetition
    abc = planefit(x[w], y[w], z[w]) ; Find a plane going trough three random points
    for k=0,2 do begin ; Run C-steps up to H_3
        res = ltsp_residuals(abc, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz)
        good = (sort(abs(res)))[0:h-1] ; Fit the h points with smallest errors
        ltsp_fitting, x[good], y[good], z[good], sigx[good], sigy[good], sigz[good], abc, sig_abc, chi_sq
    endfor    
    abcV[j,*] = abc
    chi2V[j] = chi_sq
endfor

; Perform full C-steps only for the 10 best results
;
w = sort(chi2v)
nbest = 10
chi_sq = 1d30
for j=0,nbest-1 do begin
    abc1 = abcV[w[j],*]
    repeat begin ; Run C-steps to convergence
        abcOld = abc1
        res = ltsp_residuals(abc1, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz)
        good1 = (sort(abs(res)))[0:h-1] ; Fit the h points with smallest errors
        ltsp_fitting, x[good1], y[good1], z[good1], sigx[good1], sigy[good1], sigz[good1], abc1, sig_abc1, chi1_sq
    endrep until array_equal(abcOld, abc1)
    if chi_sq gt chi1_sq then begin
        abc = abc1  ; Save best solution
        good = good1
        chi_sq = chi1_sq
    endif
endfor

end
;------------------------------------------------------------------------------
function ltsp_outliers, sig_int, $
    X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz1, ABC=abc, SIG_ABC=sig_abc, $
    H=h, GOOD=good, BAD=bad, CHI_SQ=chi_sq, OFFS=offs, DEBUG=debug, CLIP=clip
compile_opt idl2, hidden

sigz = sqrt(sigz1^2 + sig_int^2) ; Gaussian intrinsic scatter

if h eq n_elements(x) then begin ; No outliers detection

    abc = planefit(x, y, z) ; quick initial guess
    ltsp_fitting, x, y, z, sigx, sigy, sigz, abc, sig_abc, chi_sq
    bad = -1 ; No outliers

endif else begin ; Robust fit and outliers detection

    ; Initial estimate using the maximum breakdown of
    ; the method of 50% but minimum efficiency
    ;
    ltsp_algorithm, x, y, z, sigx, sigy, sigz, h, abc, good

    ; 99% confidence sigma clipping
    ;
    while 1 do begin
       
        res = ltsp_residuals(abc, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz)
        
        ;#######################plotting
        
        if KEYWORD_SET(debug) then begin
        
            zfit = abc[0] + abc[1]*x + abc[2]*y
            sigzfit = sqrt((sigx*abc[1])^2 + (sigy*abc[2])^2)
            ploterror, z, zfit, sigz, sigzfit, PSYM=4, TITLE='sig_int: '+string(sig_int,FORMAT='(g0.3)'), /NOHAT
            if h lt n_elements(res) then begin
                bad = (sort(abs(res)))[h:*]
                oploterror, z[bad], zfit[bad], sigz[bad], sigzfit[bad], PSYM=4, COLOR='red', ERRCOLOR='red', /NOHAT
            endif
            cgplot, z, z, COLOR='red', /OVERPLOT
            wait, 0.1
        
        endif

        ;#######################plotting
        
        sig = stddev(res[good])
        goodOld = good
        good = where(abs(res) lt clip*sig, h, COMPLEMENT=bad)
        if array_equal(good, goodOld) then break
        ltsp_fitting, x[good], y[good], z[good], sigx[good], sigy[good], sigz[good], abc, sig_abc, chi_sq
              
    endwhile 

endelse

; To determine 1sigma error on the intrinsic scatter the chi2
; is decreased by 1sigma=sqrt(2(h-3)) while optimizing (a,b,c) [NB: DOF=3]
;
if offs eq 1 then dchi = sqrt(2d*(h-3d)) else dchi = 0d 

return, (chi_sq+dchi)/(h-3d) - 1d
end
;------------------------------------------------------------------------------
pro ltsp_single_fit, x, y, z, sigx, sigy, sigz, par, sig_par, chi_sq, h, good, bad, clip, DEBUG=debug
compile_opt idl2, hidden

sig1 = 0d ; minimum sig_int is zero
f1 = ltsp_outliers(sig1, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz, ABC=abc, H=h, GOOD=good, OFFS=0, DEBUG=debug, CLIP=clip)
res = abc[0] + abc[1]*x + abc[2]*y - z ; Total residuals ignoring measurement errors
sig2 = stddev(res[good]) ; Maximum sig_int when all scatter is intrinsic
f2 = ltsp_outliers(sig2, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz, H=h, OFFS=0, DEBUG=debug, CLIP=clip)

fargs = {X:x, Y:y, Z:z, SIGX:sigx, SIGY:sigy, SIGZ:sigz, H:h, OFFS:0, DEBUG:debug, CLIP:clip} ; requires Chi2/DOF=1
sig_int = cap_zbrent(sig1, sig2, f1, f2, FUNC_NAME='ltsp_outliers', FUNCTARGS=fargs, SUCCESS=ok, /QUIET)
if ok eq 0 then begin
    sig_int = 0d
    sig_int_err = 0d
endif else begin
    sig1 = sig_int
    f1 = sqrt(2d/(h-3d)) ; By construction chi2=h-3 when sig=sig_int
    f2 += sqrt(2d/(h-3d)) 
    fargs = {X:x, Y:y, Z:z, SIGX:sigx, SIGY:sigy, SIGZ:sigz, H:h, OFFS:1, DEBUG:debug, CLIP:clip} ; chi2 can always decrease
    sigMax_int = cap_zbrent(sig1, sig2, f1, f2, FUNC_NAME='ltsp_outliers', FUNCTARGS=fargs, SUCCESS=ok, /QUIET)
    sig_int_err = sigMax_int - sig_int
endelse    

; Repeat fit at best-fitting location
;
tmp = ltsp_outliers(sig_int, X=x, Y=y, Z=z, SIGX=sigx, SIGY=sigy, SIGZ=sigz, $
    ABC=abc, SIG_ABC=sig_abc, CHI_SQ=chi_sq, H=h, GOOD=good, BAD=bad, OFFS=0, DEBUG=debug, CLIP=clip)
par = [abc[*],sig_int]
sig_par = [sig_abc[*],sig_int_err] ; No formal errors in sig_int

end
;------------------------------------------------------------------------------

pro lts_planefit, x0, y0, z, sigx, sigy, sigz, par, sig_par, chi_sq, $
    FRAC=frac, PLOT=plot, GOOD=good, BAD=bad, RMS=rms, $
    BOOTSTRAP=bootstrap, CLIP=clip1, NRAND=nrand, BAYES=bayes, $
    _EXTRA=ex, OVERPLOT=overplot, PIVOTX=pivotx, PIVOTY=pivoty, $
    TEXT=text, DEBUG=debug1, EPSZ=epsz
compile_opt idl2
on_error, 2

if total(n_elements(x0) eq [n_elements(y0), n_elements(z), $
          n_elements(sigx), n_elements(sigy), n_elements(sigz)]) ne 5 then $
    message, '[X, Y, Z, SIGX, SIGY, SIGZ] must have rthe same size'

debug = keyword_set(debug1)

if n_elements(clip1) eq 0 then clip = 2.6 else clip = clip1 ; Default clip=2.6*sig-->99% for Gaussian

if n_elements(pivotx) gt 0 then begin
    x = x0-pivotx 
endif else begin
    x = x0
    pivotx = 0d
endelse    

if n_elements(pivoty) gt 0 then begin
    y = y0-pivoty 
endif else begin
    y = y0
    pivoty = 0d
endelse    

p = 2 ; two dimensions
n = N_ELEMENTS(x)
if n_elements(frac) eq 0 then h = (n+p+1)/2 else h = round(frac*n) > (n+p+1)/2

ltsp_single_fit, x, y, z, sigx, sigy, sigz, par, sig_par, chi_sq, h, good, bad, clip, DEBUG=debug
rms = stddev(par[0] + par[1]*x[good] + par[2]*y[good] - z[good])

print, par, sig_par

end
;------------------------------------------------------------------------------
