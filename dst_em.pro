function dst_sqm , str , cpsf 

	sqm = fltarr(size(/n_e, *str.g) , str.ns)
	ij_sat = array_indices(*str.im , *str.s)

	for i=0, str.ns-1 do sqm[*,i] = (shift(cpsf , ij_sat[0,i] - str.npix/2 , ij_sat[1,i] - str.npix/2 ))[*str.g] 

	str.sqm = ptr_new(sqm)
	str.y	= ptr_new((*str.im)[*str.g] > 0.)
	str.bg	= ptr_new((*str.bg)[*str.g] > 0.)

	return, sqm

end

pro dst_em, str, cpsf, it=it, level=level, fft_comp=fft_comp

default , it		, str.it
default , level		, str.lev
default , fft_comp	, 0

emp_kkt = fltarr(it)
exp_kkt = fltarr(it)

if keyword_set(fft_comp) then begin 

	mask_g = *str.im*0 & mask_g[*str.g] = 1
	cpsf2  = cpsf * cpsf
	
endif else begin 
	
	sqm  = dst_sqm (str , cpsf)
	sqm2 = sqm * sqm

endelse

;;; Initialization
y 	= keyword_set(fft_comp) ? ( *str.im * mask_g )>0. : *str.y 
bg	= keyword_set(fft_comp) ? ( *str.bg * mask_g )>0. : *str.bg 
x	= keyword_set(fft_comp) ? fltarr(size(/dim, *str.im))+1 : fltarr(str.ns) + 1
V	= keyword_set(fft_comp) ? convolve(mask_g, cpsf, /corr)>0. : make_array(size(/dim, *str.g), value=1) # sqm 

for i = 0 , it-1 do begin 

	if keyword_set(fft_comp) then	y_i = mask_g * (convolve( x , cpsf ) + bg)>0 else $
									y_i = isarray(x) ? sqm # x + bg : sqm * x + bg
	
	z = f_div( y , y_i )

	U = keyword_set(fft_comp) ? (convolve( z , cpsf , /corr )>0.) : reform( z # sqm )

	if keyword_set(fft_comp) then	x[*str.s] *= ( f_div( U , V ) )[*str.s] else $ 
									x *= f_div( U , V )

	;;; STOPPING RULE
	emp_kkt[i] = total( ( x * (V - U) )^2. )
	exp_kkt[i] = keyword_set(fft_comp) ? total( x^2. * (convolve( f_div(1. , y_i * mask_g)>0 , cpsf2 , /corr )>0.) ) : $
										 total( x^2. * (f_div(1.,y_i) # sqm2) ) 
	
	if emp_kkt[i] le level * exp_kkt[i] then break
	
	print, i ,' --- ', emp_kkt[i]/exp_kkt[i] , level
	
endfor

c_stat = c_statistic( y , y_i ) 

print , i , '  R_flx:' , total( y_i ) , '  Exp_flx:' , total( y )  , '  C_stat:' , c_stat , ' ' 
 
str.x		= keyword_set(fft_comp) ? ptr_new(x[*str.s]) : ptr_new(x)
str.c_exp	= keyword_set(fft_comp) ? ptr_new(y_i[*str.g]) : ptr_new(y_i)
str.c_stat	= C_stat

end