var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MPISphericalHarmonics","category":"page"},{"location":"#MPISphericalHarmonics","page":"Home","title":"MPISphericalHarmonics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MPISphericalHarmonics.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MPISphericalHarmonics]","category":"page"},{"location":"#Base.:+-Tuple{MagneticFieldCoefficients, MagneticFieldCoefficients}","page":"Home","title":"Base.:+","text":"+(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)\n\nforce = true adds the coefficients even if the radius or center are not equal (set to values of the first coefficients).\n\n\n\n\n\n","category":"method"},{"location":"#Base.:--Tuple{MagneticFieldCoefficients, MagneticFieldCoefficients}","page":"Home","title":"Base.:-","text":"-(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)\n\nforce = true subtracts the coefficients even if the radius or center are not equal (set to values of the first coefficients).\n\n\n\n\n\n","category":"method"},{"location":"#MPISphericalHarmonics.findFFP!-Tuple{MagneticFieldCoefficients}","page":"Home","title":"MPISphericalHarmonics.findFFP!","text":"ffp = findFFP!(coeffsMF::MagneticFieldCoefficients)\n\nFind FFP and set it as coeffsMF.ffp.\n\n\n\n\n\n","category":"method"},{"location":"#MPISphericalHarmonics.findFFP-Tuple{MagneticFieldCoefficients}","page":"Home","title":"MPISphericalHarmonics.findFFP","text":"findFFP(coeffsMF::MagneticFieldCoefficients; \n        returnasmatrix::Bool=true)\n\nDescription:  Newton method to find the FFPs of the magnetic fields\n \nInput:\n\ncoeffsMF   - MagneticFieldCoefficients\n\nkwargs:\n\nreturnasmatrix - Boolean\n      true  -> return FFPs as Matrix with size (3,#Patches) (default)\n      false -> return FFPs as Array of NLsolve.SolverResults with size #Patches\n\nOutput:\n\nffp - FFPs of the magnetic field\n\n\n\n\n\n","category":"method"},{"location":"#MPISphericalHarmonics.magneticField-Union{Tuple{T}, Tuple{MPIFiles.SphericalTDesign, Union{AbstractArray{T, 3}, AbstractMatrix{T}}}, Tuple{MPIFiles.SphericalTDesign, Union{AbstractArray{T, 3}, AbstractMatrix{T}}, Int64}, Tuple{MPIFiles.SphericalTDesign, Union{AbstractArray{T, 3}, AbstractMatrix{T}}, Int64, Bool}} where T<:Real","page":"Home","title":"MPISphericalHarmonics.magneticField","text":"magneticField(tDesign::SphericalTDesign, field::Union{AbstractArray{T,2},AbstractArray{T,3}};\n\t       L::Int=Int(tDesign.T/2),\n\t       calcSolid::Bool=true) where T <: Real\n\nDescription:  Calculation of the spherical harmonic coefficients based on the measured t-design\n \nInput:\n\ntDesign\t- Measured t-design (type: SphericalTDesign)\nfield       - Measured field (size = (J,N,C)) with J <= 3\n\nkwargs:\n\nL           - Order up to which the coeffs be calculated (default: t/2)\ncalcSolid   - Boolean (default: true)\n  false -> spherical coefficients\n  true -> solid coefficients\n\nOutput:\n\ncoeffs    - spherical/solid coefficients, type: Array{SphericalHarmonicCoefficients}(3,C)\n\n\n\n\n\n","category":"method"},{"location":"#MPISphericalHarmonics.shiftFFP!-Tuple{MagneticFieldCoefficients}","page":"Home","title":"MPISphericalHarmonics.shiftFFP!","text":"shiftFFP!(coeffsMF::MagneticFieldCoefficients)\n\nShift magnetic-field coefficients into FFP (and calculate it if not available).\n\n\n\n\n\n","category":"method"}]
}
