-- Simple area search based on the range of RA and Dec getting all object
-- i band information (psf and kron flux) in the eFEDS area
-- Query website: https://hsc-release.mtk.nao.ac.jp/datasearch/?
-- Make sure to check the 'fits' radio button!


SELECT
    object_id
    , pdr3_wide.forced.ra
    , pdr3_wide.forced.dec
    , pdr3_wide.forced.a_i
	, pdr3_wide.forced2.i_psfflux_flux
	, pdr3_wide.forced2.i_psfflux_fluxerr
    , pdr3_wide.forced2.i_kronflux_flux
    , pdr3_wide.forced2.i_kronflux_fluxerr
	, pdr3_wide.forced.i_cmodel_flux
	, pdr3_wide.forced.i_cmodel_fluxerr
	, pdr3_wide.meas.i_filterfraction_weighted
FROM
    pdr3_wide.forced
    LEFT JOIN pdr3_wide.forced2 USING (object_id)
    LEFT JOIN pdr3_wide.meas USING (object_id)
WHERE
    isprimary
    AND boxSearch(coord, 126.5, 145.5, -3.0, 6.0)
          /* is equivalent to
                 ra  BETWEEN 34.0 AND 36.0
             AND dec BETWEEN -5.0 AND -4.5
             but boxSearch() is much faster
          */
;