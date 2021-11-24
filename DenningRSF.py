#----------------------------------------------------------------------------------------
#                                       fRI RESEARCH
#                                  GRIZZLY BEAR PROGRAM
#----------------------------------------------------------------------------------------
# SCRIPT NAME:  DenningRSF.py
#               v.2019.0702
#
# PURPOSE:      The tool regenerates the grizzly bear Denning RSF
#
# ARGUMENTS:    args    name        Input Description
#               [1]   deliv_yr  -   Deliverable Year
#               [2]   dest      -   Output destination folder <Required>
#               [3]   userBND   -   Area of Interest: polygon feature for which models
#                                   will be generated <Required>
#               [4]   newRoads  -   New Roads: linear features representing new
#                                   roads {Optional}
#               [5]   rplcRoads -   Replace roads: If user wants to use their complete 
#                                   road layer as opposed to supplied
#               [6]   decRoads  -   Decommissioned Roads: linear features to be removed
#                                   from landscape {Optional}
#               [7]   habDel    -   Habitat Deletions: polygon features to be removed
#                                   from landscape {Optional}
#               [8]   detailsON -   Tool option to message details on or off
#               [9]   cleanupON -   Tool option to keep or delete intermediate data
#
# OUTPUTS:      Updated 7-class Denning RSF (raster) clipped to the AOIs and table with 
#               mean RSF value
# 
# EDITORS:      Joshua Crough, GIS Analyst (based off Karine Pigeon's methods)
#                   fRI Research
#               Julie Duval, GIS Program Lead
#                   fRI Research
#
# CREATED ON:   November 2016
#
# LAST UPDATES: Jan 2018 (J.Duval) - Added deliverable year so user can select
#                                    a year to run
#               Sept 2018 (J.Duval) - Updated script so that it can be run by more than 1
#                                       user at a time (fixed file lock issue)
#               Oct 2018 (J.Duval) - Added function to check path and folder name
#               Oct 2018 (J.Duval) - Updated temp.mdb to temp.gdb. Moved shared functions
#                                   to GBToolsSharedFunctions module
#               May 2019 (J.Duval) - Created separate module to handle tool versions
#                                    Moved more shared functions to
#                                    GBToolsSharedFunctions module
#               June 2019 (J.Duval) - Updates to improve code efficiency
#                                   - Changed FeatureToRaster to PolygonToRaster in code.
#                                   - Added log file to display option
#               July 2019 (J.Duval) - Paths concatenated with os.path.join() as fp
# 
# NOTES:        written for ArcGIS 10.x (Will not work in earlier versions)
#               tested with ArcGIS 10.3.1
#               tested with ArcGIS 10.7
#               requires ArcInfo with Spatial Analyst extension.
# ---------------------------------------------------------------------------------------
# SCRIPT FUNCTIONS:
#       main(args)          - Main function that calls all other functions
#       prepData(userBND, popunit) - Preprocessing of inputs
#       clipData()          - Clip updateable inputs
#       addRoads(newRoads)  - Optional Parameter: incorporate new roads 
#       dcmRoads(decRoads)  - Optional Parameter: remove decommissioned roads
#       rmvHab(habDel)      - Optional Paremeter: remove anthropogenic features (polygon)
#       calRSF()            - Calculate Denning RSF.
#       Pdisplay(msg, log)  - Displays primary messages
#       display(msg, detailsON, log)  - Displays detailed messages (optional) 
#       error(msg)          - Displays error messages to user
# ---------------------------------------------------------------------------------------
#     * GBTpath: GBTools path, derived by sys python module in script
#     * O.: output tables and raster paths --> from Class object inside script
#     * S.: shapefile paths --> from BaseInputs module
#     * FC.: Feature Class paths --> from BaseInputs module
#     * tFC.: Temporary feature Class paths --> from Class objectg inside script
#     * R.: Raster paths --> from BaseInputs module
#     * TR.: Temporary/Intermediate Raster paths --> from BaseInputs module
#========================================================================================
import sys, string, os, datetime, traceback
from CheckLicenses import CheckArcInfo, CheckSpatialExt
from GBToolsSharedFunctions import isNameOk, isPathOk, chktmpOP, del_temp_outputs
from GBToolsSharedFunctions import Pdisplay, display, error, f
from ToolVersion import toolVersion
import arcpy
from arcpy import env
from arcpy.sa import * # For Map Algebra (Raster function)
from os.path import join as fp

class exitError(Exception):
    pass

#========================================================================================
# Import Base Inputs from User-Selected deliverable year.  Available years are
# determined on the calling tool itself.
#----------------------------------------------------------------------------------------
deliv_yr = sys.argv[1]
exec(toolVersion.delivDICT[long(deliv_yr)])

#========================================================================================
#Processing Environment
#----------------------------------------------------------------------------------------
arcpy.env.overwriteOutput = True

#========================================================================================
# Set Global Variables and tool paths
#----------------------------------------------------------------------------------------
global GBTpath, detailsON, cleanupON, tmpIP, tmpgdb, log
GBTpath = os.path.split(sys.path[0])[0]
tmpIP = fp(GBTpath, "results", "DEN_RSF", sys.argv[2], "tmpIP")
tmpgdb = fp(tmpIP, "temp.gdb")

#========================================================================================
# Test naming conventions
#----------------------------------------------------------------------------------------
try:
    isNameOk(sys.argv[2])
    isPathOk(GBTpath[0:-2])
except:
    sys.exit(0)  #terminate tool
   
#========================================================================================
# Global Tool Parameters
#----------------------------------------------------------------------------------------
# A checkbox on form for user to decide whether or not to view extra messages
# Default is yes.
detailsON = sys.argv[8]

# A checkbox on form for user to decide whether or not to delete temporary files.
# Default is yes.
cleanupON = sys.argv[9]

class O(object):
    #====================================================================================
    # OUTPUTS
    # These set the path of all the results, which are saved in user specified folders.
    #====================================================================================
    dest = sys.argv[2]
    OP = fp(GBTpath, "results", "DEN_RSF", dest)

    if os.path.exists(OP) == False and isNameOk(dest) == True:
        os.mkdir(OP)
    else:
        if os.listdir(OP):
            error("\n!!! Output folder '" + dest + "' already exists\n" +
                  "!!! Tool Exited\n")
            exit()

    # Output Rasters
    denRSF = fp(OP, "den_rsf")

    # Output Tables
    denRSF_tbl = fp(OP, "den_rsf_tbl.dbf")

class tFC(object):
    #====================================================================================
    # Temporary Feature Classes
    bnd1k = fp(tmpgdb, "bnd1k")
    roadsclip = fp(tmpgdb, "roadsclip")
    newroadsbuf = fp(tmpgdb, "newroadsbuf")
    rdrclm_fp = fp(tmpgdb, "rdrclm_fp")
    roads_upd = fp(tmpgdb, "roads_upd")
    removals = fp(tmpgdb, "removals")

#========================================================================================
# Set up file to save display log for user
#----------------------------------------------------------------------------------------
try:
    starttime = datetime.datetime.now()
    logFile = fp(os.path.split(tmpIP)[0], "tool_log.txt")
    log = open(logFile, "w+")
    log.write("Tool Started on: " + str(starttime.strftime("%d-%b-%Y %I:%M %p")) + "\n")
except:
    error("Problem opening or writing to the tool_log.txt file")
    
#----------------------------------------------------------------------------------------
# Main functions
#========================================================================================
# The main() function reads all input arguments and calls the appropriate functions.
# All paths are also set in the BaseInputs module and are prefixed with function calls
# R. (rasters), TR. (temporary rasters) S. (shapefiles), T. (tables) and
# FC. (feature classes).
def main(args):
    try:
        # Only proceed if the required ArcGIS licenses are available
        if CheckArcInfo() == "yes" and CheckSpatialExt() == "yes":
            Pdisplay("*** GBtools path detected as:\n*** " + GBTpath, log)
            #---------------------------------------------------------------------------
            # Clean up old files (if exist) and create folder called 'hstate' in temp
            # folder
            #===========================================================================
            chktmpOP(tmpgdb, tmpIP)

            # SET ENVIRONMENT VARIABLES
            #--------------------------
            arcpy.env.workspace = tmpIP
            arcpy.env.cellSize = 30
            arcpy.env.snapRaster = R.lcov_gbz
            arcpy.env.rasterStatistics = "STATISTICS 1 1"

            # User defined arguments
            #--------------------------
            userBND = sys.argv[3]
            newRoads = sys.argv[4]
            rplcRoads = sys.argv[5]
            decRoads = sys.argv[6] 
            habDel = sys.argv[7]

            # If user wants to use their own roads then replace the input roads with
            # theirs
            userRD = "N"
            ipRoads = FC.roads
            if rplcRoads == True or rplcRoads == -1 or rplcRoads == "true":
                if newRoads != "#":
                    userRD = "Y"
                    ipRoads = newRoads

            # Output dataset
            finalOut = sys.argv[8]
            
            #----------------------------------------------------------------------------
            # Section 1: Preprocessing of inputs, set analysis extent.
            #============================================================================
            Pdisplay(f.sep60 +
                     "SECTION 1: PREPROCESSING OF INPUTS" +
                     "\n     Area of Interest:\n     " + userBND, log)
            prepData(userBND)

            #----------------------------------------------------------------------------
            # Section 2: Clip updateable inputs.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 2: CLIP OUT UPDATEABLE INPUTS" +
                    "\n     ...extracting rasters within analysis extent", log)
            clipData(newRoads, decRoads, ipRoads)

            #----------------------------------------------------------------------------
            # Section 3: OPTIONAL: Incorporate new roads/add users road
            #============================================================================
            # If data user wants to add new roads then move ahead with it
            if newRoads != "#":
                if userRD == "N":
                    # If appending to current roads
                    Pdisplay(f.sep60 +
                            "SECTION 3:  OPTIONAL PARAMETER: NEW ROADS" 
                            "\n     ... incorporating new roads:"
                            "\n     " + newRoads, log)
                else:
                    # If replacing all roads
                    Pdisplay(f.sep60 +
                        "SECTION 3:  OPTIONAL PARAMETER: NEW ROADS" 
                        "\n     ... replacing current roads with:"
                        "\n     " + newRoads, log)
                addRoads(newRoads, userRD)
            else:
                Pdisplay(f.sep60 +
                        "SECTION 3:  OPTIONAL PARAMETER: NEW ROADS" 
                        "\n     ... no new roads to incorporate <section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 4: OPTIONAL: Remove decommissioned roads from clipped roads layer. 
            #============================================================================
            # If data user wants to remove reclaimed roads then move ahead with it
            if decRoads != "#":            
                Pdisplay(f.sep60 +
                        "SECTION 4:  OPTIONAL PARAMETER: DECOMMISSIONED ROADS" 
                        "\n     ... You have chosen to remove: " +
                        "\n     " + decRoads, log)
                dcmRoads(decRoads)
            else:
                Pdisplay(f.sep60 +
                        "SECTION 4:  OPTIONAL PARAMETER: DECOMMISSIONED ROADS" 
                        "\n     ... no decommissioned roads to remove <section skipped>",
                         log)

            #----------------------------------------------------------------------------
            # Section 5: OPTIONAL: Add other anthropogenic features (polygon)         
            #============================================================================
            # If data user wants to add a forecast age then move ahead with it
            if habDel != "#":            
                Pdisplay(f.sep60 +
                        "SECTION 5:  OPTIONAL PARAMETER: ADD OTHER ANTHRO FEATURES" 
                        "\n     ... you have chosen to remove: " + habDel, log)
                rmvHab(habDel)
            else:
                Pdisplay(f.sep60 +
                         "SECTION 5:  OPTIONAL PARAMETER: ADD OTHER ANTHRO FEATURES" 
                         "\n     ... no anthropogenic features to remove "
                         "<section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 6: Calculate Denning RSF.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 6:  CALCULATE DEN RSF", log)
            calRSF()

            #----------------------------------------------------------------------------
            # Section 7:  Clean temp files.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 7: CLEANING TEMP OUTPUTS", log)
            if cleanupON == True or cleanupON == -1 or cleanupON == "true":
                display("7.0 Deleting all temporary files.", detailsON, log)
                del_temp_outputs(tmpgdb, tmpIP)
            else:
                display("7.0 Keeping temporary and files in " + tmpIP + ".", detailsON,
                        log)

            #----------------------------------------------------------------------------
            # Final Message to user upon success:
            #============================================================================
            Pdisplay(f.sep60 + "FINAL OUTPUTS", log)
            Pdisplay("\nOutput folder: " + O.OP + "\n", log)
            Pdisplay("*** RASTERS ***", log)
            Pdisplay("Den RSF:      den_rsf", log)
            Pdisplay("\n*** TABLES ***", log)
            Pdisplay("Den RSF:      den_rsf_tbl.dbf", log)
            Pdisplay(f.sep60 + "DEN RSF TOOL COMPLETED\n\n", log)

            # Refresh Catalog
            arcpy.RefreshCatalog(GBTpath)

        else:
            Pdisplay("\n### The tool has not completed - please check " 
                    "that you have the required ArcGIS license levels ###\n", log)   

    except exitError:
        error("!!! ...Exiting Tool", log)
        sys.exit(0)  #terminate tool

    except arcpy.ExecuteError:
        # Get the geoprocessing error messages
        error("\nGP ERRORS:\n" + 
              arcpy.GetMessage(0) + "\n" +
              arcpy.GetMessages(2), log)
    except:
        # sys.exc_info() provides a tuple (type, value, traceback)
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in main():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
    finally:
        # Release Extension
        Pdisplay("*** Released Spatial Analyst Extension\n", log)
        arcpy.CheckInExtension("spatial")

#----------------------------------------------------------------------------------------
# 1 Preprocess inputs and set analysis extent. Accept shapefile or feature class as
#   parameter for analysis extent. Create bounding grid to extract updateable input grids
#========================================================================================
def prepData(userBND):
    try:
        #--------------------------------------------------------------------------------
        # 1.0: Clip area by model extents to ensure that areas not covered by models are
        #       excluded.
        #================================================================================
        display("1.0 Clipping of the Area of Interest (AOI) to Model extents.",detailsON)
        
        # Clip user defined area (AOI) by Caribou-Grizzly extent to ensure that areas
        # not covered by models are excluded.
        arcpy.Clip_analysis(userBND, FC.denBnd, S.AOI)

        #--------------------------------------------------------------------------------
        # 1.1: Create GRIDCODE field to populate and convert shapefile to raster
        #================================================================================
        # See if there is a GRIDCODE field in AOI, if so delete first
        fields = arcpy.ListFields(S.AOI)
        for field in fields:
            if field.name == "GRIDCODE":
                arcpy.DeleteField_management(S.AOI, "GRIDCODE")

        # Add GRIDCODE field and set value to 1
        arcpy.AddField_management(S.AOI, "GRIDCODE", "short")
        arcpy.CalculateField_management(S.AOI, "GRIDCODE", 1)

        # Convert AOI to raster with a cell size of 30
        display("1.1 Converting AOI to Raster.", detailsON, log)
        arcpy.PolygonToRaster_conversion(S.AOI, "GRIDCODE", TR.bndgrid, "CELL_CENTER")

        #--------------------------------------------------------------------------------
        # 1.2: Buffer analysis boundary by 1km for edge processing
        #================================================================================
        display("1.2 Buffering AOI by 1km.", detailsON, log)
        arcpy.Buffer_analysis(S.AOI, tFC.bnd1k, "1000", "", "", "ALL", "")

        # Set raster environment variable
        arcpy.env.mask = tFC.bnd1k
        arcpy.env.extent = tFC.bnd1k

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in prepData():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 2 Clip out updateable variables from Base Inputs, put into Temporary Inputs (tmpIP).
#   These may get overwritten.                       
#========================================================================================
def clipData(newRoads, decRoads, ipRoads):
    try:
        #--------------------------------------------------------------------------------
        # 2.0-2.1: Extract rasters to AOI.
        #================================================================================
        # Create arrays to clip data to AOI
        oMessage = ["road density","mask"]
        cRaster = [R.rdDen600, R.mask_den]
        sRaster = [TR.ned_rd, TR.mask]

        # Loop through arrays to process outputs
        for i in range(len(oMessage)):    
            display("2." + str(i) + " Extracting " + oMessage[i] + " raster.", detailsON,
                    log)
            sRast = Raster(TR.bndgrid) * Raster(cRaster[i])
            sRast.save(sRaster[i])

        # Only clip roads if needed
        if newRoads != "#" or decRoads != "#":
            display("2.2 Clipping roads to analysis extent.", detailsON, log)
            # Clip roads layer to AOI
            arcpy.Clip_analysis(ipRoads, tFC.bnd1k, tFC.roadsclip)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in clipData():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 3 OPTIONAL: Incorporate new roads
#========================================================================================
def addRoads(ipRoads, userRD):
    try:
        #--------------------------------------------------------------------------------
        # 3.0: Buffer roads by 30m.
        #================================================================================
        display("3.0 Buffering roads to create 60m right-of-way.", detailsON, log)
        arcpy.Buffer_analysis(ipRoads, tFC.newroadsbuf, "30", "", "", "ALL", "")

        #--------------------------------------------------------------------------------
        # 3.1: Create new field to populate.
        #================================================================================
        display("3.1 Adding GRIDCODE field to layer.", detailsON, log)
        # See if there is a GRIDCODE field in roads, if so delete first
        fields = arcpy.ListFields(tFC.newroadsbuf)
        for field in fields:
            if field.name == "GRIDCODE":
                arcpy.DeleteField_management(tFC.newroadsbuf, "GRIDCODE")

        # Add GRIDCODE field and set value to 1 (for mask)
        arcpy.AddField_management(tFC.newroadsbuf, "GRIDCODE", "short")
        arcpy.CalculateField_management(tFC.newroadsbuf, "GRIDCODE", 1)

        # Append new roads to existing roads (if not adding all)
        if userRD <> "Y":  
            display("   -> Appending new roads to existing road feature class.",
                    detailsON)
            arcpy.Append_management(ipRoads, tFC.roadsclip, "NO_TEST")

        #--------------------------------------------------------------------------------
        # 3.2: Union updated roads with boundary.  The field [GRIDCODE] has a value of 0
        #       to indicate it is a Right-of-Way (RoW), a value of 1 everywhere else.
        #================================================================================
        display("3.2 Unioning new roads with boundary.", detailsON, log)
        arcpy.Union_analysis([S.AOI, tFC.newroadsbuf] , S.update, "ALL", "", "GAPS")

        #--------------------------------------------------------------------------------
        # 3.3: Convert updated roads to raster.  A value of 1 to indicateit is a
        #       Right-of-Way (RoW), and a value of 0 everywhere else.
        #================================================================================
        display("3.3 Converting unioned roads to raster.", detailsON, log)
        arcpy.PolygonToRaster_conversion(S.update, "GRIDCODE_1", TR.addMask,
                                         "MAXIMUM_COMBINED_AREA")

        display("3.4 Creating road mask raster.", detailsON, log)
        # Swap values
        rAMask = Reclassify(TR.addMask, "Value", RemapValue([[0,1],[1,0]]),"DATA")
        rAMask.save(TR.change)

        display("3.5 Updating mask raster.", detailsON, log)
        rngRSF = (Raster(TR.change) * Raster(TR.mask))
        rngRSF.save(TR.maskupd)

        # Overwrite existing mask data
        arcpy.CopyRaster_management(TR.maskupd, TR.mask)

        #--------------------------------------------------------------------------------
        # 3.7: Calculate road density (600m radius).
        #================================================================================
        display("3.7 Calculating road density (within 600m).", detailsON, log)
        outLDens = LineDensity(tFC.roadsclip, "NONE", 30, 600, "SQUARE_KILOMETERS")
        outLDens.save(TR.ned_rd)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in addRoads():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 4 OPTIONAL: Remove reclaimed roads from clipped roads layer
#========================================================================================
def dcmRoads(decRoads):
    try:
        #--------------------------------------------------------------------------------
        # 4.0: Buffer roads by 30m.
        #================================================================================
        display("4.0 Buffering roads to create 60m right-of-way.", detailsON, log)
        arcpy.Buffer_analysis(decRoads, tFC.rdrclm_fp, "30", "", "FLAT", "ALL", "")
        
        #--------------------------------------------------------------------------------
        # 4.1: Create new field to populate.
        #================================================================================
        display("4.1 Adding GRIDCODE field to layer.", detailsON, log)
        # See if there is a GRIDCODE field in roads, if so delete first
        fields = arcpy.ListFields(tFC.rdrclm_fp)
        for field in fields:
            if field.name == "GRIDCODE":
                arcpy.DeleteField_management(tFC.rdrclm_fp, "GRIDCODE")

        display("   ...Removing road.", detailsON, log)
        arcpy.Erase_analysis(tFC.roadsclip, tFC.rdrclm_fp, tFC.roads_upd)

        display("   ...Updating road layer.", detailsON, log)
        arcpy.FeatureClassToFeatureClass_conversion(tFC.roads_upd, tmpgdb, "roadsclip")

        #--------------------------------------------------------------------------------
        # 4.1: Calculate road density (600m radius).
        #================================================================================
        display("4.1 Calculating road density (within 600m).", detailsON, log)
        outLDens = LineDensity(tFC.roadsclip, "NONE", 30, 600, "SQUARE_KILOMETERS") 
        outLDens.save(TR.ned_rd)
        
    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in rclmRoads():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise


def rmvHab(habDel):
    try:
        #--------------------------------------------------------------------------------
        # 5.0: Erase other anthropogenic features.
        #================================================================================
        display("5.0 Clipping " + habDel + " to analysis extent.", detailsON, log)
        # Clip anthto layer to AOI
        arcpy.Clip_analysis(habDel, tFC.bnd1k, tFC.removals)

        #--------------------------------------------------------------------------------
        # 5.1: Create new field to populate.
        #================================================================================
        display("5.1 Adding GRIDCODE field to layer.", detailsON, log)
        # See if there is a GRIDCODE field in anthro, if so delete first
        fields = arcpy.ListFields(tFC.removals)
        for field in fields:
            if field.name == "GRIDCODE":
                arcpy.DeleteField_management(tFC.removals, "GRIDCODE")

        # Add GRIDCODE field and set value to 1 (for mask)
        arcpy.AddField_management(tFC.removals, "GRIDCODE", "short")
        arcpy.CalculateField_management(tFC.removals, "GRIDCODE", 1)

        #--------------------------------------------------------------------------------
        # 5.2: Union updated anthro with boundary.  The field [GRIDCODE] has a value of 0
        #       to indicate it is a Right-of-Way (RoW), a value of 1 everywhere else.
        #================================================================================
        display("5.2 Unioning new anthro with boundary.", detailsON, log)
        arcpy.Union_analysis([S.AOI, tFC.removals] , S.update, "ALL", "", "GAPS")

        #--------------------------------------------------------------------------------
        # 5.3: Convert updated anthro to raster.  A value of 1 to indicateit is a
        #       anthro, and a value of 0 everywhere else.
        #================================================================================
        display("5.3 Converting antrho features to raster.", detailsON, log)
        arcpy.PolygonToRaster_conversion(S.update, "GRIDCODE_1", TR.addMask,
                                         "MAXIMUM_COMBINED_AREA")

        display("5.4 Creating anthro mask raster.", detailsON, log)
        # Swap values
        rAMask = Reclassify(TR.addMask, "Value", RemapValue([[0,1],[1,0]]),"DATA")
        rAMask.save(TR.change)

        display("5.5 Updating mask raster.", detailsON, log)
        rngRSF = (Raster(TR.change) * Raster(TR.mask))
        rngRSF.save(TR.maskupd)

        # Overwrite existing mask data
        arcpy.CopyRaster_management(TR.maskupd, TR.mask)
         
    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in rmvHab():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise
    
#----------------------------------------------------------------------------------------
# 6 Calculate Denning RSF                    
#========================================================================================
def calRSF():
    try:
        #--------------------------------------------------------------------------------
        # 6.0: Calculate raw denning RSF.
        #================================================================================
        display("6.0 Calculating raw denning RSF.", detailsON, log)
        # Coefficients
        c1 = 0.05;  c2 = -1.69;  c3 = -0.02
        c4 = -0.84;   c5 = 0.02;   c6 = -0.02

        # Reclassification break points
        bp1 = 0.0001; bp2 = 0.007843; bp3 = 0.015686
        bp4 = 0.023529; bp5 = 0.031373; bp6 = 0.043137

        # RSF equation
        tmpRSF = Exp((Raster(R.slope) * (c1)) + (Raster(TR.ned_rd) * (c2)) +
                     (Raster(R.fallFood) * (c3)) + (Raster(R.wet9600m) * (c4)) +
                     (Raster(R.springFood) * (c5)) + (Raster(R.decid4800m) * (c6)))
        tmpRSF.save(TR.rsf_raw)

        #-----------------------------------------------------------------------------
        # 6.1 Rescale to between 0 and 1,the rescaled RSF is (rawRSF - min) / (max - min)
        #=============================================================================
        display("6.1 Rescaling RSF between 0 and 1.", detailsON, log)
        # Get raster min/max properties
        #ipRawRSF = arcpy.Raster(TR.rsf_raw)
        # NOTE: The rescale equation used for this doesn't work well when area clipped,
        #        therefore, the min and max values are taken directly from the
        #        full raw rsf raster, and hardcoded in the code.
        #        WHEN UPDATING BACKGROUND LAYERS IT IS ADVISED TO TAKE THE MIN AND
        #        MAX VALUES FROM UPDATED RAW RSF AND ADD NEW VALUES BELOW
        #rMIN = float(arcpy.GetRasterProperties_management(ipRawRSF, "MINIMUM").getOutput(0))
        #rMAX = float(arcpy.GetRasterProperties_management(ipRawRSF, "MAXIMUM").getOutput(0))
        rMIN = 5.2881943268801E-10
        rMAX = 67.2001571655273

        rsTmpRSF = (Raster(TR.rsf_raw) - rMIN)/(rMAX - rMIN)
        rsTmpRSF.save(TR.rsf_rscld)

        #-----------------------------------------------------------------------------
        # 6.2 Classify seasonal RSF output into 7 classes
        #=============================================================================
        display("6.2 Reclassifying rescaled RSF to 7 classes.", detailsON, log)
        sRSF = Reclassify(TR.rsf_rscld, "VALUE",
                          RemapRange([[0,bp1,0], [bp1,bp2,1], [bp2,bp3,2], [bp3,bp4,3],
                                      [bp4,bp5,4], [bp5,bp6,5], [bp6,1,6]]), "DATA")
        sRSF.save(TR.rsf_rng)

        #-----------------------------------------------------------------------------
        # 6.3 Apply masks: Slopes < 4.5 degree and anthropogenic features
        # (roads, wells, etc)
        #=============================================================================
        display("6.3 Removing non-candidate areas for denning.", detailsON, log)
        rngRSF = (Raster(TR.rsf_rng) * Raster(TR.mask))
        rngRSF.save(O.denRSF)

        #--------------------------------------------------------------------------------
        # 6.4: Calculate mean rsf_max.
        #================================================================================
        display("6.4 Calculating mean rsf_max.", detailsON, log)
        mRSFm = ZonalStatisticsAsTable(TR.bndgrid, "Value", O.denRSF, O.denRSF_tbl,
                                       "DATA", "MEAN")

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in calRSF():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise


#========================================================================================
#Beginning of output and code
#========================================================================================
try:
    #------------------------------------------------------------------------------------
    # Script Identifier
    #------------------------------------------------------------------------------------
    Pdisplay(f.sep60 +
             f.sep25 + "       fRI Research" +
             f.sep25 + "   Grizzly Bear Program" +
             f.sep25 + "      TOOL: DEN RSF" +
             f.sep25 +
             f.sep60, log)

    main(sys.argv)

except (KeyboardInterrupt, SystemExit):
    error("!!! Tool Exited\n", log)

finally:
    del arcpy
    endtime = datetime.datetime.now()
    log.write("Tool Completed on: " + str(endtime.strftime("%d-%b-%Y %I:%M %p")) +
              "\nElapsed Time: " +
              str(round(((endtime - starttime).total_seconds()/60.0), 2)) + " minutes\n")
    log.close()
#========================================================================================
