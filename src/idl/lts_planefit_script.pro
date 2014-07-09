pro lts_planefit_script

args = command_line_args()
input = args[0]

readcol, input, x,y,z,sigx,sigy,sigz, f='d,d,d,d,d,d'

abc = [1, 2, 3]
sig_abc = [0.1, 0.2, 0.3]

lts_planefit, x, y, z, sigx, sigy, sigz, abc, sig_abc, chi2, $
    /PLOT, /TEXT, /EPSZ, XTITLE='z', YTITLE=textoidl('a + b (x-x_0) + c (y-y_0)'), $
    PIVOTX=median(x), PIVOTY=median(y)

end
