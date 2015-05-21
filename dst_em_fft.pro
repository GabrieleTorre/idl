pro dst_em_fft , str , cpsf , it=it , level=level 

; NAME:
;	dst_em_fft
; PURPOSE:
;	Deconvolution/desaturation routine based on EM method implemented by means the FFT.
; EXPLANATION:
;	This routine has a multiple usage in the desaturation method:
;		1) deconvolution method routine for the not saturated image;
;		2) compute the correlation product for determining primary saturated pixels; 
;		3) desaturation routine in the case of large saturated areas ( > 1000 );
; CALLING SEQUENCE:
;	dst_em_fft , str , cpsf , it=it , level=level , pow=pow
; INPUTS:
;	str		= dst_str global structure for the desaturation routine.	
;	cpsf	= diffraction component of the PSF
; OPTIONAL:
;	it		= maximum iteration number for EM method (default 300)
;	level	= EM stopping rule parameter (default 0.01)
; OUTPUT:	
;	str.x 		= retrieved intensities for saturated pixels
;	str.y 		= array corresponding to the diffraction fringes intensities
;	str.c_exp	= expectation computed by the method (cpsf*x + bg at the last iteration)
;	str.c_stat 	= c_statistic computed between y and c_exp					
; CALLS:
;	CALLED BY:
;		DST_X_Y
;		DECONV
;		DESATURATION
;	CALLS:
;		CONVOLVE, C_SATATISTIC, AVG
;
	
default , it , str.it
default , level , str.lev
default , pow , 1

dim = size(/dim, *str.im)

s   = *str.s
g	= *str.g

mask_g  = fltarr(dim) & mask_g[g] = 1
x		= fltarr(dim) & x[s] = TAG_EXIST(str,'bg_bkup') ? (*str.bg_bkup)[s] : 1

mask = *str.im * 0 & mask[dim[0]/2,dim[1]/2] = 1
cpsf = convolve(mask , cpsf)>0
cpsf2 = cpsf * cpsf

emp_kkt = fltarr(it+1)
exp_kkt = fltarr(it+1)

y 	= ( *str.im * mask_g )>0.
bkg = ( *str.bg * mask_g )>0.				 ;;; GT modification

V = convolve(mask_g , cpsf , /corr )>0.

for I = 0, it-1 do begin

	y_i = (convolve( x , cpsf ) + bkg)>0
	
	Z = (f_div(y , y_i)) * mask_g

	U = convolve( z , cpsf , /corr )>0.
	x[s] *= ( f_div( U , V ) )[s]

	emp_kkt[i] = total( ( x * (V - U))^2.)
	exp_kkt[i] = total( x^2. * (convolve( f_div(1. , y_i * mask_g)>0 , cpsf2 , /corr )>0.) )

	if emp_kkt[i] le level * exp_kkt[i] then break
	
	print, i ,' --- ', emp_kkt[i]/exp_kkt[i] , level
	
endfor

C_stat = c_statistic( y , y_i ) 
 
print , i , '  R_flx:' , total( y_i[g]) , '  Exp_flx:' , total( y[g] )  , '  C_stat:' , c_stat , ' ' 

str.x = ptr_new(x[s])
str.y = ptr_new(y)
str.c_exp  = ptr_new(y_i[g])
str.c_stat = c_stat					

end
