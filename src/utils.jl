"""
    loadMagneticFieldMeasurementData(filename::String)

Load magnetic field measurement data by choosing the correct function for the version of the file and return all required data.
"""
function loadMagneticFieldMeasurementData(filename::String)

  measurementData = h5open(filename, "r") do file

    # test if the file contains the correct data
    if haskey(file, "positions") && file["positions"] isa HDF5.Group # newest file structure
      measurementData = loadMagneticFieldMeasurementDataV2(filename)
    elseif haskey(file, "positions") && !(file["positions"] isa HDF5.Group)# oldest file structure
      measurementData = loadMagneticFieldMeasurementDataV1(filename)
    else # file structure not known
      throw(ErrorException("Unknown file structure."))
    end
  end

  # get required data
  field = measurementData["fields"]
    # positions
    tDes = measurementData["positions"]["tDesign"]
    radius = tDes["radius"]
    N = tDes["N"]
    t = tDes["t"]
    center = tDes["center"]

  if haskey(measurementData, "sensor")
    correction =  measurementData["sensor"]["correctionTranslation"]
  else
    correction = Matrix{Float64}(I, 3, 3)
    @warn "No correction of the sensor translations found. Using identity matrix."
  end

  return field, radius, N, t, center, correction
end

"""
    loadMagneticFieldMeasurementDataV2(filename::String)

Load magnetic field measurement data from hdf5-file (version v2).
"""
function loadMagneticFieldMeasurementDataV2(filename::String)

  measurementData = h5open(filename, "r") do file
    measurementData = Dict{String, Any}()

    # magnetic field
    measurementData["fields"] = read(file, "/fields")   # measured field (size: 3 x #points x #patches)

    # positions
    measurementData["positions"] = Dict{String, Any}()
    measurementData["positions"]["tDesign"] = Dict{String, Any}()
    tDes = measurementData["positions"]["tDesign"]
      tDes["radius"] = read(file, "/positions/tDesign/radius") # radius of the measured ball
      tDes["N"] = read(file, "/positions/tDesign/N")           # number of points of the t-design
      tDes["t"] = read(file, "/positions/tDesign/t")           # t of the t-design
      tDes["center"] = read(file, "/positions/tDesign/center") # center of the measured ball
    
    # sensor
    measurementData["sensor"] = Dict{String, Any}()
    measurementData["sensor"]["correctionTranslation"] = read(file, "/sensor/correctionTranslation")

    # optional data
    if haskey(file, "currents")
      measurementData["currents"] = read(file, "currents") # currents used for the measurement
    end

    return measurementData
  end

  return measurementData
end

"""
    loadMagneticFieldMeasurementDataV1(filename::String)

Load magnetic field measurement data from hdf5-file (version v1).
"""
function loadMagneticFieldMeasurementDataV1(filename::String)

  measurementData = h5open(filename, "r") do file
    measurementData = Dict{String, Any}()

    measurementData["fields"] = read(file, "fields")   # measured field (size: 3 x #points x #patches)
    measurementData["fieldsError"] = read(file,"fieldsError") 	# field error stemming from the measurement device

    # positions
    measurementData["positions"] = Dict{String, Any}()
    measurementData["positions"]["tDesign"] = Dict{String, Any}()
    tDes = measurementData["positions"]["tDesign"]
      tDes["positions"] = read(file,"positions") 		# measured positions (shifted and scaled t-design)
      tDes["radius"] = read(file, "positionsTDesignRadius") # radius of the measured ball
      tDes["N"] = read(file, "positionsTDesignN")           # number of points of the t-design
      tDes["t"] = read(file, "positionsTDesignT")           # t of the t-design
      tDes["center"] = read(file, "positionsCenter") # center of the measured ball

    # optional data
    if haskey(file, "currents")
      measurementData["currents"] = read(file, "currents") # currents used for the measurement
    end

    return measurementData
  end
  
  return measurementData
end