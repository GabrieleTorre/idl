;pro demo_for_desaturation
;;; Demo for the desaturation procedure for AIA images

;;; Cutout data request line.
;ssw_cutout_service,'2011/09/06 22:00:00','2011/09/06 22:30:00',ref_helio='N15W18',fovx=350,fovy=350,email='torre@dima.unige.it',waves=[94,131,171,193,211,304,335],max_frames=1000,instrument='aia',aec=1

; NOTE: in the code images will be resized to 499 * 499 pixels FOV by ctrl_data function.
; So, a bigger FOV for the input data is needed.

tstart = '15-feb-2011 01:45:00'						; start time selection
tend   = '15-feb-2011 01:50:00'						; end time selection

path = '15feb'										; data storage folder path

;wav	= ['94','131','171','193','211','304','335']	; wavelength to process
wav	= ['131'] ;,'171','193','211','304','335']	; wavelength to process

obj = obj_new('desat')

result = obj -> desaturation( wav , tstart , tend  , path , save_fts = 1 , aec = 1, loud = 1, lev=1.0 )

end