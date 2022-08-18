-- VHS query for the eFEDS field
-- -> YHJKs information
-- Query website: http://horus.roe.ac.uk:8080/vdfs/VSQL_form.jsp
-- (Note: Select VHS: VISTA Hemisphere Survey in the survey box!)
SELECT
    ra, dec, pStar, pGalaxy,
    yAperMag4, yAperMag4Err, yAperMag6, yAperMag6Err,
    jAperMag4, jAperMag4Err, jAperMag6, jAperMag6Err,
    hAperMag4, hAperMag4Err, hAperMag6, hAperMag6Err,
    ksAperMag4, ksAperMag4Err, ksAperMag6, ksAperMag6Err,
    eBV, aY, aJ, aH, aKs
FROM
    vhsSource
WHERE
    ra BETWEEN 126.5 AND 146.2
    AND dec BETWEEN -3.2 AND 6.2
    AND priOrSec=0