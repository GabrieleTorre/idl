pro dst_disp_factor_2, str , mult_fact , psf

; NAME:
;	dst_disp_factor_2
; PURPOSE:
;	compute the diffraction component to subtract from the image, the convolution between 
;	the retrieved intensties, inside the saturation region, and the central core of the 
;	PSF and the intensities for the bloomed pixels.
; EXPLANATION:
;	1) Computes the convolution product between the intensities obtained from the EM method 	
;		and the central core of the psf. 
;	2) The diffraction fringes generated by the retrieved intensities are computed and 
;		subtracted from the original image
;	3) The bloomed intensities are substituted by the estimated background values.
; CALLING SEQUENCE:
;	dst_disp_factor_2, str , mult_fact , psf
; INPUTS:
;	str		  =	dst_str global structure for the desaturation routine.	
;	mult_fact = multiplicative factor obtained in desaturation routine necessary to 
;				renormalize the estimated background, minimizing the interpolation error 
;				made during the interpolation of the zero-th frequency component in 
;				dst_bkg_fourier_interp.pro
;	psf 	  = structure containing all the psf components, diffraction, dispersion and 
;				the complete one. This structure is returned as output from dst_psf_gen.pro
;				routine.	
;	OUTPUT:
;	It fills the x tag of str
; CALLS:
;	CALLED BY:
;		DESATURATION
;	CALLS:
;		CONVOLVE
; PROCEDURE:
;	

im = *str.im
s  = *str.s 
b  = *str.b

indx = findgen(size(/n_e , im))

result = im 

check = where(b gt 0 , ct)

if ct gt 1 then begin
	
	f =  where( histogram(indx	,min=0,max=size(/n_e , im)) 		- $
				histogram(b		,min=0,max=size(/n_e , im))			- $
				histogram(s   	,min=0,max=size(/n_e , im)) eq 1) 

endif else begin 

	f =  where( histogram(indx	,min=0,max=size(/n_e , im)) 		- $
				histogram(s   	,min=0,max=size(/n_e , im)) eq 1) 

endelse

if ct gt 0 then result[b] = mult_fact*(convolve(*str.bg_bkup,psf.opsf)>0)[b]
	
;;; Necessary step to avoid border artifacts during fringes subtraction step

mask_temp = im*0+1 & mask_temp[f] = 0

ct_off = 0.01

f = where(convolve(mask_temp , psf.opsf)>0 lt ct_off)

tmp = im*0 & tmp[s] = *str.x

result[f] -= (convolve(tmp,psf.cpsf)>0)[f]
result[s] = (convolve(tmp,psf.opsf)>0)[s]

*str.x = result>0.

end
