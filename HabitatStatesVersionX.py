#----------------------------------------------------------------------------------------
#                                      fRI RESEARCH
#                                  GRIZZLY BEAR PROGRAM
#                                       GBTOOLS
#----------------------------------------------------------------------------------------
# SCRIPT NAME:  HabitatStates.py
#               v.X
# 
# PURPOSE:      The tool regenerates the grizzly bear Habitat States model (combination
#               of the RSF and mortality risk models) from base landscape layers.
#
# ARGUMENTS:    args    name        Input Description
#               [1]   deliv_yr  -   Deliverable Year
#               [2]   dest      -   Output destination folder <Required>
#               [3]   userBND   -   Area of Interest: polygon feature for which models
#                                   will be generated <Required>
#               [4]   popunit   -   Population Unit: selected from a list.  Determines
#                                   which coefficients to use for RSF <Required>
#               [5]   newClearing - New vegetation clearing such as cutblocks, well
#                                   pads and other clearing types. {Optional}
#               [6]   clearingyr -  Year of clearing defined in newClearing attribute
#                                   table {Optional}
#               [6]   newHerb   -   New linear features categorized as "herb" such as
#                                   pipelines or trails {Optional} 
#               [7]   newRoads  -   New Roads: linear features representing new
#                                   roads {Optional}
#               [8]   recRoads  -   Reclaimed Roads: linear features to be removed
#                                   from landscape.  User has option to change footprint
#                                   to barren, upland herb, upland treed or same as
#                                   surrounding habitat {Optional}
#               [9]   recRdType
#               [10]  habDel    -   Habitat Deletions: polygon features to be removed
#                                   from landscape.  RSF score = 0 {Optional}
#               [11]  fcstyr   -    Forecast year: used to model succession of new stands.
#                                   Blocks are no longer considered regeneration after
#                                   50 years. {Optional}
#               [12]  detailsON -   Tool option to message details on or off
#               [13]  cleanupON -   Tool option to keep or delete intermediate data
#               [14]  process   -   Choice of process: Habitat States, RSF or Risk
#
# OUTPUTS:      Updated habitat models (rasters) clipped to the AOI - mortality risk, 
#               RSF for spring/summer/fall, maximum RSF and habitat states.
#               Statistics (tables) for the AOI, including mean value of habitat states.
# 
# EDITORS:      Jerome Cranston, GIS Specialist
#                   Arctos Ecological Services
#               Julie Duval, GIS Program Lead
#                   fRI Research
#               Josh Crough, GIS Analyst
#                   fRI Research
#               Dan Wismer, GIS Analyst
#                   fRI Research
#
# CREATED ON:   December 2011
#
# LAST UPDATES: Mar 2012 (J.Cranston)  -   amalgamation of earlier scripts
#                                          RSFcalc_p7_93.py and RiskCalc_p6_93.py 
#               May 2012 (J.Cranston)  -   changed femalerng, protctd6mi, and whitez6mi 
#                                          from float to integer, adjusted coefficients.
#               Sep 2013 (J.Duval)  - updated script formatting and split all
#                                     sections into functions.
#               Apr 2014 (J.Duval)  - testing and troubleshooting completed.
#               Apr 2015 (J.Crough) - update script formatting to arcpy
#               May 2017 (J.Duval)  - added layer files for outputs
#                                   - users have a choice to run only RSF or Risk
#               Jan 2018 (J.Duval)  - Added deliverable year so user can select
#                                     a year to run
#               Sep 2018 (J.Duval)  - Updated script so that it can be run by more than 1
#                                     user at a time (fixed file lock issue)
#               Oct 2018 (J.Duval)  - Added function to check path and folder name
#               Oct 2018 (J.Duval)  - Updated temp.mdb to temp.gdb. Moved shared functions
#                                     to GBToolsSharedFunctions module
#               Feb 2019 (J.Duval)  - Update how linear features are treated in the model.
#                                     The user has a choice on how to define reclaimed
#                                     roads.
#               May 2019 (J.Duval)  - Created separate module to handle tool versions
#                                     Moved more shared functions to
#                                     GBToolsSharedFunctions module
#               Jun 2019 (J.Duval)  - Updates to improve code efficiency
#                                   - Changed FeatureToRaster to PolygonToRaster in code.
#                                   - Removed cellsize of 10m that was used in previous
#                                     versions of the tool.
#                                   - Added log file to display option
#               Jul 2019 (J.Duval)  - Paths concatenated with os.path.join() as fp
#               Nov 2019 (D.Wismer) - Fixed appending new lines to roads and trails
#                                   - Fixed raster calculation for regen age
#               Mar 2020 (D.Wismer) - Added & commented out Chin coefficents
#               Aug 2020 (D.Wismer) - Updated newClearing function to accept clearing 
#                                     years and forecast year
#                                   - Refactored forecast age function to accept a
#                                     forecast year to predict output results
#               Dec 2020 (D.Wismer) - Converted arcpy.env.mask to TR.bndgrid
#                                   - Converted CON function to Raster Calc syntax
#                                   - Refactored section 6, reclaimed roads match option
#                                   - Deleted clearinylyr from in_memory 
#
# NOTES:        written for ArcGIS 10.x (Will not work in earlier versions)
#               tested with ArcGIS 10.3.1
#               tested with ArcGIS 10.7
#               tested with ArcGIS 10.8
#               requires ArcInfo with Spatial Analyst extension.
# ---------------------------------------------------------------------------------------
# SCRIPT FUNCTIONS:
#       main(args)          - Main function that calls all other functions
#       prepData(userBND, popunit) - Preprocessing of inputs
#       clipData()          - Clip updateable inputs
#       addClearing(newClearing, clearingyr, forecastyr)
#       Optional Parameter: incorporate new clearings, clearing year and forecast year
#       addLnHerb(newHerb)  - Optional Parameter: incorporate new linear features (herb)
#       addRoads(newRoads)  - Optional Parameter: incorporate new roads 
#       rclmRoads(recRoads, recRdType) - Optional Parameter: remove reclaimed roads
#       fAge(fcstyr)        - Optional Parameter: forecast year
#       secInputs()         - Derive Secondary Inputs
#           -> distFEdge(ipLC,opCont,aBorder,opEuc,opMA) - within secInputs function
#       calRSF(popUnit, habDel) - Calculate RSF. Optional Parameter: habitat deletions
#       calRisk()           - Calculate RISK
#       calHState()         - Calculate HABITAT STATE
#       createLayerFiles()  - Created output layer files for the output rasters
#       Pdisplay(msg, log)  - Displays primary messages
#       display(msg, detailsON, log)  - Displays detailed messages (optional) 
#       error(msg, log)     - Displays error messages to user
#       prepLines()         - Preps new input lines
# ---------------------------------------------------------------------------------------
#     * GBTpath: GBTools path, derived by sys python module in script
#     * O.: output tables and raster paths --> from Class object inside script
#     * S.: shapefile paths --> from BaseInputs module
#     * FC.: Feature Class paths --> from BaseInputs module
#     * tFC.: Temporary feature Class paths --> from Class objectg inside script
#     * R.: Raster paths --> from BaseInputs module
#     * TR.: Temporary/Intermediate Raster paths --> from BaseInputs module
#     * f. : Formatting for displaying output separator lines
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
# Processing Environment
#----------------------------------------------------------------------------------------
arcpy.env.overwriteOutput = True
arcpy.Delete_management("in_memory") # Clear in_memory workspace

#========================================================================================
# Set Global Variables tool paths
#----------------------------------------------------------------------------------------
global GBTpath, detailsON, cleanupON, tmpIP, tmpgdb, log
GBTpath = os.path.split(sys.path[0])[0]
tmpIP = fp(GBTpath, "results", "HAB_STATE", sys.argv[2], "tmpIP")
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
detailsON = sys.argv[13]

# A checkbox on form for user to decide whether or not to delete temporary files.
# Default is yes.
cleanupON = sys.argv[14]

# A choice on form for user to decide which process to run.
# Default is 1. habitat states.
process = int(sys.argv[15][0])

class O(object):
    #====================================================================================
    # OUTPUTS
    # These set the path of all the results, which are saved in user specified folders.
    dest = sys.argv[2]
    OP = fp(GBTpath, "results", "HAB_STATE", dest)

    if os.path.exists(OP) == False:
        os.mkdir(OP)
    else:
        if os.listdir(OP):
            error("\n!!! Output folder '" + dest + "' already exists\n" +
                  "!!! Tool Exited\n")
            exit()
    
    # Output Rasters
    risk = fp(OP, "risk")
    risk_class = fp(OP, "risk_class")
    rsf_class = fp(OP, "rsf_class")
    rsf_max = fp(OP, "rsf_max")
    rsf_s1 = fp(OP, "rsf_s1")
    rsf_s2 = fp(OP, "rsf_s2")
    rsf_s3 = fp(OP, "rsf_s3")
    habState = fp(OP, "hab_state")

    # Output Tables
    drd_fut_tbl = fp(OP, "drd_fut.dbf")
    risk_tbl =  fp(OP, "risk.dbf")
    rsf_max_tbl = fp(OP, "rsf_max.dbf")
    habState_tbl = fp(OP, "habState.dbf")

class tFC(object):
    #====================================================================================
    # Temporary Feature Classes
    bnd1k = fp(tmpgdb, "bnd1k")
    bnd1k_arc = fp(tmpgdb, "bnd1k_arc")
    border1 = fp(tmpgdb, "border1")
    border2 = fp(tmpgdb, "border2")
    border3 = fp(tmpgdb, "border3")
    border4 = fp(tmpgdb, "border4")
    border5 = fp(tmpgdb, "border5")
    #harvestclip = fp(tmpgdb, "harvestclip")
    clearingclip = fp(tmpgdb, "clearingclip")
    edgnvg = fp(tmpgdb, "edgnvg")
    edgrgn = fp(tmpgdb, "edgrgn")
    edguph = fp(tmpgdb, "edguph")
    edgutree = fp(tmpgdb, "edgutree")
    edgwtree = fp(tmpgdb, "edgwtree")
    foredge = fp(tmpgdb, "foredge")
    newroadsbuf = fp(tmpgdb, "newroadsbuf")
    #newpipesbuf = fp(tmpgdb, "newpipesbuf")
    newherbbuf = fp(tmpgdb, "newherbbuf")
    rdrclm_fp = fp(tmpgdb, "rdrclm_fp")
    removals = fp(tmpgdb, "removals")
    roadsclip = fp(tmpgdb, "roadsclip")
    roads_upd = fp(tmpgdb, "roads_upd")
    trailsclip = fp(tmpgdb, "trailsclip")
    rclmRdBuf = fp(tmpgdb, "rclmRdBuf")

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
            #============================================================================
            # Set which habitat states process to follow
            Pdisplay("*** GBtools path detected as:\n*** " + GBTpath, log)
            if process == 2:
                Pdisplay("*** You have selected to process RSF ouptuts only\n", log)
            if process == 3:
                Pdisplay("*** You have selected to process Risk ouptuts only\n", log)

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
            popunit = sys.argv[4] 
            newClearing = sys.argv[5]
            clearingyr = sys.argv[6]
            newHerb = sys.argv[7] 
            newRoads = sys.argv[8] 
            recRoads = sys.argv[9]
            recRdType = sys.argv[10]
            habDel = sys.argv[11] 
            fcstyr = sys.argv[12]  

            #----------------------------------------------------------------------------
            # Section 1: Preprocessing of inputs, set analysis extent.
            #============================================================================
            Pdisplay(f.sep60 +
                     "SECTION 1: PREPROCESSING OF INPUTS" +
                     "\n     Area of Interest:\n     " + userBND +
                     "\n     Population Unit: " + popunit, log)
            prepData(userBND, popunit)

            # Set raster environment variable
            arcpy.env.mask = TR.bndgrid

            #----------------------------------------------------------------------------
            # Section 2: Clip updateable inputs.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 2: CLIP OUT UPDATEABLE INPUTS" +
                    "\n     ...extracting rasters within analysis extent", log)
            clipData()

            #----------------------------------------------------------------------------
            # Section 3: OPTIONAL: Incorporate new vegetation clearings
            # (harvest blocks, well pads)
            #============================================================================
            # If user wants to include vegetation clearings, then move ahead with it
            if newClearing != "#":            
                Pdisplay(f.sep60 +
                        "SECTION 3:  OPTIONAL PARAMETER: VEGETATION CLEARINGS" 
                        "\n     ... Adding new vegetation clearings to regen raster:"
                        "\n     " + newClearing, log)
                addClearing(newClearing, clearingyr, fcstyr)
            else:
                Pdisplay(f.sep60 +
                         "SECTION 3:  OPTIONAL PARAMETER: NEW VEGETATION CLEARINGS" 
                         "\n     ... no new vegetation clearings to incorporate "
                         "<section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 4: OPTIONAL: Incorporate new line features as herb
            # (pipelines, trails)
            #============================================================================
            # If user wants to add new linear features as herb then move ahead with it
            if newHerb != "#":            
                Pdisplay(f.sep60 +
                        "SECTION 4:  OPTIONAL PARAMETER: NEW LINE FEATURES AS HERB" 
                        "\n     ... incorporating new line features as herb:"
                        "\n     " + newHerb, log)
                addLnHerb(newHerb)
            else:
                Pdisplay(f.sep60 +
                        "SECTION 4:  OPTIONAL PARAMETER: NEW LINE FEATURES AS HERB" 
                        "\n     ... no new line features as herb to incorporate "
                        "<section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 5: OPTIONAL: Incorporate new roads
            #============================================================================
            # If user wants to add new roads then move ahead with it
            if newRoads != "#":            
                Pdisplay(f.sep60 +
                        "SECTION 5:  OPTIONAL PARAMETER: NEW ROADS" 
                        "\n     ... incorporating new roads:"
                        "\n     " + newRoads, log)
                addRoads(newRoads)
            else:
                Pdisplay(f.sep60 +
                        "SECTION 5:  OPTIONAL PARAMETER: NEW ROADS" 
                        "\n     ... no new roads to incorporate <section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 6: OPTIONAL: Remove reclaimed roads from clipped roads layer.  
            #============================================================================
            # If user wants to remove reclaimed roads then move ahead with it
            if recRoads != "#" and (process == 1 or process == 3):            
                Pdisplay(f.sep60 +
                        "SECTION 6:  OPTIONAL PARAMETER: RECLAIMED ROADS" 
                        "\n     ... You have chosen to remove: " +
                        "\n     " + recRoads + " as " + recRdType, log)
                rclmRoads(recRoads, recRdType)
            else:
                Pdisplay(f.sep60 +
                        "SECTION 6:  OPTIONAL PARAMETER: RECLAIMED ROADS" 
                        "\n     ... no reclaimed roads to remove <section skipped>", log)
                if process == 2:
                    Pdisplay("     ... reclaimed roads are not considered for RSF model",
                             log)

            #----------------------------------------------------------------------------
            # Section 7: OPTIONAL: Forecast age.          
            #============================================================================
            # If user wants to add a forecast age then move ahead with it
            if int(fcstyr) != int(deliv_yr): # Relative to BaseInputs
                age = int(fcstyr) - int(deliv_yr)
                Pdisplay(f.sep60 +
                        "SECTION 7:  OPTIONAL PARAMETER: FORECAST YEAR", log)
                if age > 1:
                    Pdisplay("     ... assigning forecast age of " + str(age), log)
                    fAge(age)
                else:
                   # if the forecasted age is 1 year, the CC equation will give a
                   # negative value.
                   # therefore, it will be assumed 0 (neglibeable effect)
                   Pdisplay("\n     ... no changes when forecast age = 1 ", log)

            else:
                Pdisplay(f.sep60 +
                        "SECTION 7:  OPTIONAL PARAMETER: FORECAST YEAR" 
                        "\n     ... no forecast age assigned <section skipped>", log)

            #----------------------------------------------------------------------------
            # Section 8: Derive secondary inputs.
            #============================================================================
            Pdisplay(f.sep60 +
                    "SECTION 8:  DERIVE SECONDARY INPUTS", log)
            secInputs()

            # Set raster environment variable
            arcpy.env.mask = TR.bndgrid

            #----------------------------------------------------------------------------
            # Section 9: Calculate RSF.
            #            OPTIONAL: Habitat deletions.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 9:  CALCULATE RSF", log)
            if process == 1 or process == 2:
                calRSF(popunit, habDel)
            else:
                Pdisplay("     ... user selected to run for Risk only <section skipped>",
                         log)

            #----------------------------------------------------------------------------
            # Section 10: Calculate Risk.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 10:  CALCULATE RISK", log)
            if process == 1 or process == 3:
                calRisk()

            #----------------------------------------------------------------------------
            # Section 11: Calculate Habitat States.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 11:  CALCULATE HABITAT STATES", log)
            if process == 1:
                calHState()
            else:
                Pdisplay("     ... user selected to run for Risk or RSF only.  Habitat\n"
                         "          States will not be calculated. <section skipped>",
                         log)

            #----------------------------------------------------------------------------
            # Section 12:  Clean temp files.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 12:  CLEAN TEMPORARY OUTPUTS", log)
            if cleanupON == True or cleanupON == -1 or cleanupON == "true":
                display("12.0 Deleting all temporary files.", detailsON, log)
                del_temp_outputs(tmpgdb, tmpIP)
            else:
                display("12.0 Keeping intermediate outputs in " + tmpIP + ".", detailsON,
                        log)

            #----------------------------------------------------------------------------
            # Section 13: Create output layer files.
            #============================================================================
            Pdisplay(f.sep60 + "SECTION 13:  CREATE OUTPUT LAYER FILES", log)
            createLayerFiles()
                
            #----------------------------------------------------------------------------
            # Final Message to user upon success:
            #============================================================================
            Pdisplay(f.sep60 + "FINAL OUTPUTS", log)
            Pdisplay("\nOutput folder: " + O.OP + "\n", log)
            Pdisplay("*** OUTPUT RASTERS ***", log)
            if process == 1:
                Pdisplay("Habitat States: hab_state", log)
            if process == 1 or process == 2:
                Pdisplay("RSF Class:      rsf_class", log)
                Pdisplay("RSF Season 1:   rsf_s1", log)
                Pdisplay("RSF Season 2:   rsf_s2", log)
                Pdisplay("RSF Season 3:   rsf_s3", log)
                Pdisplay("RSF Max:        rsf_max", log)
            if process == 1 or process == 3:
                Pdisplay("Risk:           risk", log)
                Pdisplay("Risk Class:     risk_class", log)
            Pdisplay("\n*** OUTPUT TABLES ***")
            if process == 1:
                Pdisplay("Habitat States: habState.dbf", log)
            if process == 1 or process == 2:
                Pdisplay("RSF Max:        rsf_max.dbf", log)
            if process == 1 or process == 3:
                Pdisplay("Risk:           risk.dbf", log)
            if process == 1:
                Pdisplay("Future Mean Distance to Roads:         drd_fut.dbf", log)
            Pdisplay(f.sep60 + "HABITAT STATES COMPLETED\n\n", log)

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
def prepData(userBND, popunit):
    try:
        #--------------------------------------------------------------------------------
        # 1.0: Clip area by model extents to ensure that areas not covered by models are
        #       excluded.
        #================================================================================
        display("1.0 Clipping Area of Interest (AOI) to Model extents.", detailsON, log)
        
        # Clip user defined area (AOI) by Caribou-Grizzly extent to ensure that areas
        # not covered by models are excluded.
        arcpy.Clip_analysis(userBND, FC.GBZ_bnd, S.AOI)

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

        #--------------------------------------------------------------------------------
        # 1.3: Convert the buffered AOI to a line
        #================================================================================
        display("1.3 Converting the buffered AOI to a line.", detailsON, log)
        arcpy.PolygonToLine_management(tFC.bnd1k, tFC.bnd1k_arc)

        #--------------------------------------------------------------------------------
        # Section 1.4 The buffered boundary is converted to 5 line feature classes
        #             (borders).  These will later be unioned with the edges between
        #             land cover classes to ensure the distance-to-edge grids extend to
        #             the full analysis extent.
        #================================================================================
        display("1.4 Creating borders for the distance to edge analysis", detailsON, log)
        for x in ['1','2','3','4','5']:
            arcpy.FeatureClassToFeatureClass_conversion(tFC.bnd1k_arc,tmpgdb,"border" +x)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in prepData():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        error("\nAn error may occur if the AOI falls outside the GB Tool Boundary\n",log)
        raise
    
#----------------------------------------------------------------------------------------
# 2 Clip out updateable variables from Base Inputs, put into Temporary Inputs (tmpIP).
#   These may get overwritten.                       
#========================================================================================
def clipData():
    try:
        #--------------------------------------------------------------------------------
        # 2.0-2.6: Extract rasters to AOI.
        #================================================================================
        # Create arrays to clip data to AOI
        oMessage = ["landcover", "regen", "regenage", "percent conifer",
                    "default negative decay (100) distance to trails", 
                    "default negative decay (500) distance to trails",
                    "default negative decay distance to roads"]
        cRaster = [R.lcov_gbz, R.regen_gbz, R.regenage, R.pctcon_gbz,
                   R.ned_tr100, R.ned_tr500, R.ned_rd]
        sRaster = [TR.lcov, TR.rgn, TR.rgnage, TR.pctcon,
                   TR.ned_tr100, TR.ned_tr500, TR.ned_rd]

        # Loop through arrays to process outputs
        if process == 2: j = [0, 1, 2, 3] #RSF Only
        else: j = range(len(oMessage))
        
        for i in j:    
            display("2." + str(i) + " Extracting " + oMessage[i] + " raster.", detailsON,
                    log)
            sRast = Raster(TR.bndgrid) * Raster(cRaster[i])
            sRast.save(sRaster[i])

        #--------------------------------------------------------------------------------
        # 2.7: Extract crown closure raster to regen and treed pixels
        #================================================================================
        display("2.7 Extracting crown closure.", detailsON, log)
        
        display("   ... to regen", detailsON, log)
        cc_rgn = Con(Raster(TR.rgn) == 1, R.cc_gbz, 0)
        cc_rgn.save(TR.cc_rgn)

        display("   ... to treed", detailsON, log)
        cc_trd = Con(Raster(TR.rgn) == 0, R.cc_gbz, 0)
        cc_trd.save(TR.cc_trd)

        # Create arrays to clip road/trail data to AOI
        oMessage = ["roads", "trails"]
        ipFC = [FC.roads, FC.trails]
        opFC = [tFC.roadsclip, tFC.trailsclip]

        #--------------------------------------------------------------------------------
        # 2.8: Clip roads and trails to AOI
        #================================================================================

        # Loop through arrays to process outputs
        n = 8
        for i in range(len(oMessage)):    
            display("2." + str(n) + " Clipping " + oMessage[i] + " to analysis extent.",
                    detailsON, log)

            # Clip layer to AOI
            arcpy.Clip_analysis(ipFC[i], tFC.bnd1k, opFC[i])

            # Append to border 2 data
            display("   2." + str(n) + ".1 Appending border to " + oMessage[i] + ".",
                    detailsON, log)
            arcpy.Append_management(tFC.border2, opFC[i], "NO_TEST")          
            n = n+1

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in clipData():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        error("\nAn error may occur if an input raster is missing data within the AOI\n",
              log)
        raise
    
#----------------------------------------------------------------------------------------
# 3 OPTIONAL: Timber harvest scenario
#========================================================================================
def addClearing(newClearing, clearingyr, fcstyr):
    try:
        #--------------------------------------------------------------------------------
        # 3.0: Clip harvest blocks by model extents
        #================================================================================
        display("3.0 Clipping input clearings to analysis extent.", detailsON, log)
        arcpy.Clip_analysis(newClearing, S.AOI, tFC.clearingclip)

        #--------------------------------------------------------------------------------
        # 3.1: Create AGE, REGEN, CC & GRIDCODE fields in clipped feature to populate
        #================================================================================
        display("3.1 Calculating clearing AGE, REGEN, CC and GRIDCODE.", detailsON, log)

        # Add new fields
        arcpy.AddField_management(tFC.clearingclip, "AGE", "short")
        arcpy.AddField_management(tFC.clearingclip, "REGEN", "short")
        arcpy.AddField_management(tFC.clearingclip, "CC", "short")
        arcpy.AddField_management(tFC.clearingclip, "GRIDCODE", "short")
        arcpy.AddField_management(tFC.clearingclip, "_EXISTS", "short")

        if clearingyr != "#":
            fields = [clearingyr, "AGE", "REGEN", "CC", "GRIDCODE", "_EXISTS"]
            with arcpy.da.UpdateCursor(tFC.clearingclip,fields) as cursor:
                for row in cursor:
                    # AGE: relative to 2018 BaseInputs
                    age = 2018 - row[0]
                    row[1] = age
                    
                    # REGEN: relative to forecast year
                    if (int(fcstyr) - row[0]) >= 0 and (int(fcstyr) - row[0]) <= 50: 
                        row[2] = 1
                    else:
                        row[2] = 0 # clearing does not yet exist or is now treed
                        
                    # CC: relative to age
                    CC = -3.0977 + (2.4784 * abs(age)) - (0.0218 * abs(age) * abs(age))
                    if CC < 0:
                        row[3] = 0
                    else:
                        row[3] = CC
                    
                    # GRIDCODE (landcover): relative to age
                    if age <= 3:
                        row[4] = 7 # barren
                    elif age > 3 and age <= 17:
                        row[4] = 3 # upland herb
                    elif age > 17 and age <= 30:
                        row[4] = 3 # shrub
                    elif age > 30:
                        row[4] = 1 # upland tree

                    # EXISTS: Record if user-clearings exists
                    # i.e. if fcstyr = 2025 and clearingyr = 2030: user-clearing != exist
                    if (int(fcstyr) - row[0]) >= 0:
                        row[5] = 1
                    else:
                        row[5] = 0
                        
                    # Update row    
                    cursor.updateRow(row)
                # Clean cursor    
                del cursor, row
                    
        else:
            # If no year is provided, clearing is assumed to be a 2018 addition
            arcpy.CalculateField_management(tFC.clearingclip, "AGE", 0)
            arcpy.CalculateField_management(tFC.clearingclip, "REGEN", 1)
            arcpy.CalculateField_management(tFC.clearingclip, "CC", 0)
            arcpy.CalculateField_management(tFC.clearingclip, "GRIDCODE", 7)
            arcpy.CalculateField_management(tFC.clearingclip, "_EXISTS", 1)
                        
        #--------------------------------------------------------------------------------
        # 3.2: Convert clearing attributes to raster
        #================================================================================
        display("3.3 Converting clearing attributes to raster ", detailsON, log)
        arcpy.MakeFeatureLayer_management(tFC.clearingclip,"clearinglyr","_EXISTS = 1")

        display ("   ... regen")
        clearing_rgn = arcpy.PolygonToRaster_conversion("clearinglyr", "REGEN",
                                  "in_memory//clr_rgn","MAXIMUM_COMBINED_AREA")

        display ("   ... regen age")
        clearing_age = arcpy.PolygonToRaster_conversion("clearinglyr", "AGE",
                                  "in_memory//clr_rgnage","MAXIMUM_COMBINED_AREA")
        
        display ("   ... crown closure")
        clearing_cc = arcpy.PolygonToRaster_conversion("clearinglyr", "CC",
                                 "in_memory//clr_cc","MAXIMUM_COMBINED_AREA")
        display ("   ... land cover")
        clearing_lcov = arcpy.PolygonToRaster_conversion("clearinglyr", "GRIDCODE",
                                       "in_memory//clr_lcov","MAXIMUM_COMBINED_AREA")
        # Delete feature layer
        arcpy.Delete_management("clearinglyr")
        
        #--------------------------------------------------------------------------------
        # 3.3: Update BaseInputs with clearing attributes
        #================================================================================
        display("3.3 Updating inputs with clearing attribute rasters ", detailsON, log)

        display ("   ... regen ")
        arcpy.Mosaic_management([TR.rgn, clearing_rgn], TR.rgn, "LAST")
        arcpy.BuildRasterAttributeTable_management(TR.rgn, "Overwrite")

        display ("   ... regen age")
        arcpy.Mosaic_management([TR.rgnage, clearing_age], TR.rgnage, "LAST")

        display ("   ... crown closure in regen")
        arcpy.Mosaic_management([TR.cc_rgn, clearing_cc], TR.cc_rgn, "LAST")

        display ("   ... crown closure in treed")
        cc_trd = Con(Raster(TR.rgn) == 0, R.cc_gbz, 0)
        cc_trd.save(TR.cc_trd)

        display ("   ... land cover")
        arcpy.Mosaic_management([TR.lcov, clearing_lcov], TR.lcov, "LAST")
        
    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in addClearing():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 4 OPTIONAL: Incorporate new line features as herb (pipelines, trails)
#========================================================================================
def addLnHerb(newHerb):
    try:
        prepLines("4", newHerb, tFC.newherbbuf, tFC.trailsclip, 3)

        #--------------------------------------------------------------------------------
        # 4.5: Reclassify the updated raster to create a Pick order.
        #      Where change = 3 (new  RoW), priority = 1. 
        #      Where change = 0 (everywhere other than RoW), priority = 2 so picks 2nd
        #       raster in list = lcov
        #================================================================================
        display("4.5 Creating priority grid.", detailsON, log)
        rOrder = Reclassify(TR.change, "Value", RemapValue([[3,1],[0,2]]),"DATA")
        rOrder.save(TR.order)

        #--------------------------------------------------------------------------------
        # 4.6-.10: Update background rasters to account for new features.
        #================================================================================
        # Create arrays to update rasters
        # NOTE: Fixed error in previous years tools where the pctcon output was actually
        #       the cc_trd output -> gp.copyraster_management(cc_trdupd, newpctcon)
        #       (J. Crough)
        oMessage = ["landcover", "crown closure in regen", "crown closure in treed",
                    "regen", "percent conifer"]
        rPick = [[TR.change, TR.lcov], [0,TR.cc_rgn], [0,TR.cc_trd],
                 [0,TR.rgn], [0,TR.pctcon]]
        sRaster = [TR.lcovupd, TR.cc_rgnupd, TR.cc_trdupd, TR.rgnupd, TR.pctconupd]
        owRaster = [TR.lcov, TR.cc_rgn, TR.cc_trd, TR.rgn, TR.pctcon]
        n = 6

        # Loop through arrays to process outputs
        for i in range(len(oMessage)):    
            display("4." + str(n) + " Updating " + oMessage[i] + " raster.", detailsON,
                    log)
            sRast = Pick(TR.order, rPick[i])
            sRast.save(sRaster[i])

            # Overwrite existing data
            arcpy.CopyRaster_management(sRaster[i], owRaster[i])
            n = n+1
            
        #--------------------------------------------------------------------------------
        # 4.11: Calculate cost distance to linear herb features. Only do this if new
        #       features are added, otherwise use clipped ned_tr (100, 500).
        #================================================================================
        if process == 1 or process == 3:
            display("4.11 Calculating proximity to linear features classed "
                    "as upland herb.", detailsON, log)
            cdTrails = CostDistance(tFC.trailsclip, R.tri_cost)
            cdTrails.save(TR.dtrails)

            display("   ...decay factor of 100.", detailsON, log)
            # Decay factor of 100
            n100 = Exp(-0.01 * Raster(TR.dtrails))
            n100.save(TR.ned_tr100)

            display("   ...decay factor of 500.", detailsON, log)
            # Decay factor of 500
            n500 = Exp(-0.002 * Raster(TR.dtrails))
            n500.save(TR.ned_tr500)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in addLnHerb():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 5 OPTIONAL: Incorporate new roads
#========================================================================================
def addRoads(newRoads):
    try:
        prepLines("5", newRoads, tFC.newroadsbuf, tFC.roadsclip, 7) 

        #--------------------------------------------------------------------------------
        # 5.5: Reclassify the updated roads raster to create a Pick order.
        #      Where change = 7 (new  RoW), priority = 1. 
        #      Where change = 0 (everywhere other than RoW), priority = 2 so picks 2nd
        #       raster in list = lcov
        #================================================================================
        display("5.5 Creating priority grid.", detailsON, log)
        rOrder = Reclassify(TR.change, "Value", RemapValue([[7,1],[0,2]]),"DATA")
        rOrder.save(TR.order)

        #--------------------------------------------------------------------------------
        # 5.6-.9: Update background rasters to account for new roads.
        #================================================================================
        # Create arrays to update rasters
        oMessage = ["landcover", "crown closure in regen",
                    "crown closure in treed", "regen"]
        rPick = [[TR.change,TR.lcov], [0,TR.cc_rgn], [0,TR.cc_trd], [0,TR.rgn]]
        sRaster = [TR.lcovupd, TR.cc_rgnupd, TR.cc_trdupd, TR.rgnupd]
        owRaster = [TR.lcov, TR.cc_rgn, TR.cc_trd, TR.rgn]
        n = 6

        # Loop through arrays to process outputs
        for i in range(len(oMessage)):    
            display("5." + str(n) + " Updating " + oMessage[i] + " raster.", detailsON,
                    log)
            sRast = Pick(TR.order,rPick[i])
            sRast.save(sRaster[i])

            # Overwrite existing data
            arcpy.CopyRaster_management(sRaster[i], owRaster[i])
            n = n+1

        #--------------------------------------------------------------------------------
        # 5.10: Calculate cost distance to roads.
        #================================================================================
        if process == 1 or process == 3:
            display("5.10 Calculating proximity to roads.", detailsON, log)
            cdRoads = CostDistance(tFC.roadsclip, R.tri_cost)
            cdRoads.save(TR.droads)

            display("   ...decay factor of 100.", detailsON, log)
            # Decay factor of 100
            n100 = Exp(-0.01 * Raster(TR.droads))
            n100.save(TR.ned_rd)        

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in addRoads():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 6 OPTIONAL: Remove reclaimed roads from clipped roads layer
#========================================================================================
def rclmRoads(recRoads, recRdType):
    try:
        if recRdType == "#":
            recRdType = "Barren"

        #--------------------------------------------------------------------------------
        # 6.0: Assign land cover type to type of reclaimed road selected by user.
        #================================================================================
        display("6.0 Assigning land cover type " + recRdType, detailsON, log)
        if recRdType == "Barren": lc = 7; ulc = 0
        elif recRdType == "Herb": lc = 3; ulc = 0 # upland herb
        else: lc = 99; ulc = 0 # match neighbours

        #--------------------------------------------------------------------------------
        # 6.1: Buffer reclaimed roads
        #================================================================================
        display("6.1 Buffering reclaimed roads by 30m.", detailsON, log)
        arcpy.Buffer_analysis(recRoads, tFC.rclmRdBuf, "30", "", "FLAT", "ALL", "")
            
        # Add GRIDCODE field and set value to lc
        arcpy.AddField_management(tFC.rclmRdBuf, "GRIDCODE", "short")
        arcpy.CalculateField_management(tFC.rclmRdBuf, "GRIDCODE", lc)

        #--------------------------------------------------------------------------------
        # 6.2: Union updated features with boundary.  The field [GRIDCODE] has a value of
        #      "lc" and a value of 0 everywhere else.
        #================================================================================
        display("6.2 Unioning new features with boundary.", detailsON, log)
        arcpy.Union_analysis([S.AOI, tFC.rclmRdBuf] , S.update, "ALL", "", "GAPS")

        #--------------------------------------------------------------------------------
        # 6.3: Convert updated features to raster.  A value of "lc" to indicate it is a
        #       Right-of-Way (RoW), and a value of 0 everywhere else.
        #================================================================================
        display("6.3 Converting features to raster.", detailsON, log)
        arcpy.PolygonToRaster_conversion(S.update, "GRIDCODE_1", TR.change,
                                         "MAXIMUM_COMBINED_AREA")

        #--------------------------------------------------------------------------------
        # 6.4: Update background rasters to account for new features.
        #================================================================================
        if lc == 99:
            display("6.4 Matching neighboring cells in rasters", detailsON, log)
            oMessage = ["landcover", "crown closure in regen", "crown closure in treed",
                        "regen", "percent conifer"]
            sRaster = [TR.lcovupd, TR.cc_rgnupd, TR.cc_trdupd, TR.rgnupd, TR.pctconupd]
            owRaster = [TR.lcov, TR.cc_rgn, TR.cc_trd, TR.rgn, TR.pctcon]
            n = 1

            TR.mask = Reclassify(TR.change,"Value",RemapValue([[0,1],[99,"NoData"]]))

            # Loop through arrays to process outputs
            # Matching surrounding cells
            for i in range(len(oMessage)):
                arcpy.ClearEnvironment("mask") # avoids .vat schema locks 
                display("6.4." + str(n) + " Updating " + oMessage[i] + " raster.",
                        detailsON, log)
                if i in  [0, 3]:
                    stat_type = "MAJORITY"
                else:
                    stat_type = "MEAN"
                    
                sRast = Con(Raster(TR.change) == 99,FocalStatistics(Times(owRaster[i],
                        TR.mask),NbrAnnulus(1,3,"CELL"), stat_type, "DATA"),owRaster[i])
                
                # Back fill no-data pixels in focal majority
                if i in [0]:
                   fillRast = Con(IsNull(sRast),3,sRast) # lc
                   fillRast.save(sRaster[i])
                elif i in [3]:
                    fillRast = Con(IsNull(sRast),0,sRast) # rgn
                    fillRast.save(sRaster[i])
                else:
                    sRast.save(sRaster[i])

                # Overwrite existing data
                arcpy.CopyRaster_management(sRaster[i], owRaster[i])
                n = n + 1
            # Reset mask variable
            arcpy.env.mask = TR.bndgrid         
        else:
            #----------------------------------------------------------------------------
            # 6.4: Reclassify the updated raster to create a Pick order.
            #      Where change = "lc", priority = 1. 
            #      Where change = 0, priority = 2 so picks 2nd raster in list = lcov
            #============================================================================
            display("6.4 Creating priority grid.", detailsON, log)
            rOrder = Reclassify(TR.change, "Value", RemapValue([[lc,1],[0,2]]),"DATA")
            rOrder.save(TR.order)

            #----------------------------------------------------------------------------
            #  Update background rasters to account for new features.
            #============================================================================
            # Create arrays to update rasters
            # NOTE: Fixed error in previous years tools where the pctcon output was
            #       actually the cc_trd output -> gp.copyraster_management(cc_trdupd,
            #       newpctcon)(J. Crough)
            oMessage = ["landcover", "crown closure in regen", "crown closure in treed",
                        "regen", "percent conifer"]
            rPick = [[TR.change, TR.lcov], [ulc, TR.cc_rgn], [ulc, TR.cc_trd],
                     [ulc, TR.rgn], [ulc, TR.pctcon]]
            sRaster = [TR.lcovupd, TR.cc_rgnupd, TR.cc_trdupd, TR.rgnupd, TR.pctconupd]
            owRaster = [TR.lcov, TR.cc_rgn, TR.cc_trd, TR.rgn, TR.pctcon]
            n = 1

            # Loop through arrays to process outputs
            for i in range(len(oMessage)):    
                display("6.4." + str(n) + " Updating " + oMessage[i] + " raster.",
                        detailsON, log)
                sRast = Pick(TR.order, rPick[i])
                sRast.save(sRaster[i])

                # Overwrite existing data
                arcpy.CopyRaster_management(sRaster[i], owRaster[i])
                n = n + 1        

        #--------------------------------------------------------------------------------
        # 6.5: Erase reclaimed roads from roads layer.
        #================================================================================
        display("6.5 Erasing reclaimed roads from road layer.\n"
                "   ...Removing road from road layer.", detailsON, log)
        arcpy.Erase_analysis(tFC.roadsclip, tFC.rclmRdBuf, tFC.roads_upd)
        arcpy.FeatureClassToFeatureClass_conversion(tFC.roads_upd, tmpgdb, "roadsclip")

        #--------------------------------------------------------------------------------
        # 6.6: Generate cost distance to updated roads.
        #================================================================================
        display("6.6 Generating cost distance to updated roads.", detailsON, log)
        if int(arcpy.GetCount_management(tFC.roadsclip)[0]) > 1:
            cdURoads = CostDistance(tFC.roadsclip, R.tri_cost)
            cdURoads.save(TR.drd_fut)

            display("   ...decay factor of 100.", detailsON, log)
            # Decay factor of 100
            n100 = Exp(-0.01 * Raster(TR.drd_fut))
            n100.save(TR.ned_rd)  

        else: # All roads in AOI have been removed
            display("   ...all roads have been reclaimed.")
            cdURoads = Reclassify(TR.bndgrid, "Value", RemapValue([[1,0]]),"DATA")
            cdURoads.save(TR.drd_fut)
            arcpy.CopyRaster_management(TR.drd_fut, TR.ned_rd)
            
        #--------------------------------------------------------------------------------
        # 6.7: Summarize the values of a the cost distance raster within the analysis
        #      boundary. Output to table.  All of the statistics will be calculated.  
        #================================================================================
        display("6.7 Summarize cost distance values within the analysis extent\n"
                "    -> output table: " + O.drd_fut_tbl, detailsON, log)
        zsURoads = ZonalStatisticsAsTable(TR.bndgrid, "Value", TR.drd_fut,
                                          O.drd_fut_tbl, "DATA")

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in rclmRoads():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 7 OPTIONAL: Increment based on forecast age all regen areas will be tested; if their
#               forecast age > 50, they will revert to upland tree Crown closure in regen
#               will be incremented
#========================================================================================
def fAge(age):
    try:
        #--------------------------------------------------------------------------------
        # 7.0: Add forecast age to current block age. Multiply by regen to ensure that
        #      forecast age is not added to non-regen areas, i.e. where regen is 0.
        #================================================================================
        display("7.0 Adding forecast age to current block age.", detailsON, log)
        aFA = age + Raster(TR.rgnage) 
        aFAcon = Con(aFA > 9999, 0, aFA)
        aFAcon.save(TR.rgnage_fcst)

        #--------------------------------------------------------------------------------
        # 7.1: Update regen: remove blocks > 50 years from regen as they are no longer
        #       regen.
        #================================================================================
        display("7.1 Removing harvest blocks greater than 50 years from regen raster.",
                detailsON, log)
        rBlock = Con((Raster(TR.rgnage_fcst) > 50)|(Raster(TR.rgnage_fcst)) == 0, 0, 1)
        rBlock.save(TR.rgn_fcst)
        arcpy.CopyRaster_management(TR.rgn_fcst, TR.rgn)

        #--------------------------------------------------------------------------------
        # 7.2: Update crown closure from regen and from treed: increment cc
        #       Note: 0.0142857 is inverse of 70
        #================================================================================
        display("7.2 Increment crown closure in future regen.", detailsON, log)

        CC = -3.0977 + (2.4784 * age) - (0.0218 * age * age)  # from cc vs age curve
        display("      ...canopy of new blocks at forecast age is " + str(CC) + "\n" +
                "      ...incrementing crown closure in regen raster.", detailsON, log)
        cellStats = CellStatistics([(70 - Raster(TR.cc_rgn)), 0], "MAXIMUM", "DATA")
        cClsure = Raster(TR.cc_rgn) + (CC * cellStats * 0.0142857)
        cClsure.save(TR.cc_rgnupd)
        #================================================================================
        # JEROMES NOTES FROM ~2009
        #--------------------------------------------------------------------------------
        # 1) for cc = 0 (all new blocks), cc is incremented by strCC (second term is 1)
        # 2) for cc > 0 and cc < 70, cc is incremented by a fraction of strCC
        #       (as cc_rgn >> 70, increment >> 0)
        # 3) the MAX function ensures that pixel values over 70, which make the second
        #       term negative, do not cause CC to decrease.
        # 4) cc_rgn is 0 outside of regen too, so all areas get the increment...
        #================================================================================

        # order: all regen areas >50 will have value=1, all other areas = 2
        # assign incremented cc to regen areas > 50
        display("   7.2.0 Creating position raster.", detailsON, log)
        pRast = Con(Raster(TR.rgnage_fcst) > 50, 1, 2)
        pRast.save(TR.order)

        # for all regen areas > 50, new cc_trd value will be taken from cc_rgnupd;
        # elsewhere, from cc_trd
        display("   7.2.1 Updating crown closure in treed raster.", detailsON, log)
        upTcc = Pick(TR.order, [TR.cc_rgnupd, TR.cc_trd])
        upTcc.save(TR.cc_trdupd)
        arcpy.CopyRaster_management(TR.cc_trdupd, TR.cc_trd)

        # Reduce regen extent
        display("   7.2.2 Updating crown closure regen extent.", detailsON, log)
        upRE = Raster(TR.cc_rgnupd) * Raster(TR.rgn)
        upRE.save(TR.cc_rgn)

        #--------------------------------------------------------------------------------
        # 7.3: Update percent conifer raster.  For all regen areas > 50, new  pctcon
        #       value will be 90; elsewhere, from p7pctcon.
        #      NOTE: Updated by J. Crough (April 2015): Original was overwriting clipped
        #            percent conifer raster with original (R.pctcon_gbz/pctcon_gbz).
        #            This would overwrite any modifications to it made in Section 4.10
        #            (adding pipelines) if both optional parameters were used.
        #               -> Use TR.pctcon in Pick rather than R.pctcon_gbz
        #================================================================================
        display("7.3 Updating percent conifer raster.", detailsON, log)
        upPC = Con(Raster(TR.order) == 1, 90, TR.pctcon) 
        upPC.save(TR.pctconupd)
        arcpy.CopyRaster_management(TR.pctconupd, TR.pctcon)

        #--------------------------------------------------------------------------------
        # 7.4: Update landcover: Change blocks > 30 years to upland tree. Blocks > 30
        #       and < 50 are still regen but are upland tree.
        #      NOTE (Jerome): Could also increment cc first and use a cc cutoff value.
        #                     Should check whether there is a better correlation of cc
        #                     and upland tree blocks
        #================================================================================
        display("7.4 Updating landcover: change blocks > 30 years to upland tree.",
                detailsON, log)
        pRast = Con(Raster(TR.rgnage_fcst) > 30, 1, 2)
        pRast.save(TR.order)

        #  order: all regen areas > 30 will have value = 1, all other areas = 2
        #  intermediate classes (barren, upland herb, shrub dont matter - they all go
        #  to non-forest for risk)
        upLC = Con(Raster(TR.order) == 1, 1, TR.lcov) 
        upLC.save(TR.lcovupd)
        arcpy.CopyRaster_management(TR.lcovupd, TR.lcov)

        #================================================================================
        # JEROMES NOTES
        #--------------------------------------------------------------------------------
        # sc_x_utree has always been zero everywhere in regen, because utree is 0. But
        #   now that old regen is reverting back to utree, sc needs a new value. So
        #   p7pctcon will have to be clipped by the AOI, and wherever regen has
        #   reverted to not regen (>50), sc will be replaced with a new value (90?)
        #   This new value should depend on elevation, natural subregion, etc.
        #================================================================================

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in fAge():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 8 Derive secondary inputs                    
#========================================================================================
def secInputs():
    try:
        #--------------------------------------------------------------------------------
        # 8.0: Distance to forest edge: dfor_rscld (risk).  Reclassify landcover to
        #       forest ({class 1 and 2} --> 1) vs non-forest ({class 3 - 10} --> 0).
        #================================================================================
        display("8.0 Reclassifying landcover.", detailsON, log)
        # Create array's for Reclassifying the landcover
        lblReclass = ["forest/nonforest", "regen", "nonveg", "upland herb", "wet treed", 
                      "upland treed", "wet herb", "shrub"]
        rvReclass = [[[1,2,1],[3,10,0]],[[1,0],[0,1]],[[1,5,0],[6,10,1]],
                     [[1,2,0],[3,3,1],[4,10,0]],[[1,1,0],[2,2,1],[3,10,0]],
                     [[1,1,1],[2,10,0]],[[1,3,0],[4,4,1],[5,10,0]],
                     [[1,4,0],[5,5,1],[6,10,0]]]
        opReclass = [TR.forest, TR.erase, TR.noveg, TR.uherb, TR.wtree, TR.utree, 
                     TR.wherb, TR.shrub]

        if process != 3: processRange = range(len(lblReclass))
        else: processRange = [0]
        for i in processRange:
            n = i + 1
            if lblReclass[i] == "regen":
                display("   8.0." + str(n) + " Removing landcover classes within " +
                        lblReclass[i] + ".", detailsON, log)

                # Reclassify landcover
                display("      ...reclassifying.", detailsON, log)
                reclass = Reclassify(TR.rgn, "VALUE", RemapValue(rvReclass[i]), "DATA")
                reclass.save(opReclass[i])

                display("      ...extracting.", detailsON, log)
                lRSF = Raster(opReclass[i]) * Raster(TR.lcov)
                lRSF.save(TR.rsflcov)
            else:
                display("   8.0." + str(n) + " Extracting " + lblReclass[i] +
                        " from landcover.", detailsON, log)
                if i < 2:
                    rcLayer = TR.lcov
                else:
                    rcLayer = TR.rsflcov 

                # Reclassify landcover
                reclass = Reclassify(rcLayer, "VALUE", RemapRange(rvReclass[i]), "DATA")
                reclass.save(opReclass[i])

        #--------------------------------------------------------------------------------
        # 8.1: Process landcover
        #================================================================================
        display("8.1 Processing landcover.", detailsON, log)
        # Build arrays for further processing
        lblProcess = ["regen", "upland herb", "wet treed", "forest/nonforest",
                      "upland treed", "non-veg"]
        ipLC = [TR.rgn, TR.uherb, TR.wtree, TR.forest, TR.utree, TR.noveg]
        opCont = [tFC.edgrgn, tFC.edguph, tFC.edgwtree, tFC.foredge,
                  tFC.edgutree, tFC.edgnvg]
        aBorder = [tFC.border2, tFC.border3, tFC.border4, tFC.border1, tFC.border5,
                   tFC.border1]
        opEuc = [TR.dedge2, TR.dedge3, TR.dedge4, TR.dfor, TR.dedge5, TR.dedge1]
        opMA = [TR.dfor_ergn, TR.dfor_uhrb, TR.dfor_wtree, TR.dfor_rscld, TR.dfor_utree,
                TR.dfor_envg]

        # Go through and process values in arrays
        for i in range(len(ipLC)):
            n = i + 1
            display("   8.1." + str(n) + " Processing " + lblProcess[i] + ".", detailsON,
                    log)
            # for first 3 records need to ensure landtype exists in analysis extent
            if i < 3 and process != 3:
                numrecs = arcpy.GetCount_management(ipLC[i])
                if int(numrecs.getOutput(0)) > 1:
                    display("      ..." + lblProcess[i] + " has " + str(numrecs) +
                            " records.", detailsON, log)
                    # If there is more than one record process distance to forest
                    # edge (function)
                    distFEdge(ipLC[i], opCont[i], aBorder[i], opEuc[i], opMA[i])
                else:
                    display("      ..." + lblProcess[i] + " has " + str(numrecs) +
                            " record.", detailsON, log)

                    # If there is only one record reclassify raster
                    if i == 0:
                        # if regen is empty, use as dfor_ergn
                        display("      ...using regen grid for output.", detailsON, log)
                        arcpy.CopyRaster_management(ipLC[i], opMA[i])
                    else:
                        # reclassify to 0
                        display("      ...reclassifying " + lblProcess[i] + " to zero.",
                                detailsON, log)
                        reclass = Reclassify(ipLC[i],"VALUE",RemapValue([[1,0]]), "DATA")
                        reclass.save(opMA[i])
            else:
                # after the first 3 records calculate distance to forest edge (function)
                if process != 3:
                    distFEdge(ipLC[i], opCont[i], aBorder[i], opEuc[i], opMA[i])
                else:
                    if i == 3:
                        distFEdge(ipLC[i], opCont[i], aBorder[i], opEuc[i], opMA[i])

        #--------------------------------------------------------------------------------
        # 8.2: Derive species composition in upland treed classes.
        #================================================================================
        display("8.2 Deriving species composition in upland treed classes.", detailsON,
                log)
        # Calculate by multiplying upper treed raster with the percent conifer raster
        if process != 3:
            calSUT = Raster(TR.utree) * Raster(TR.pctcon)
            calSUT.save(TR.sc_x_utree)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in secInputs():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 9 Calculate RSF                    
#========================================================================================
def calRSF(popunit, habDel):
    try:
        #--------------------------------------------------------------------------------
        # 9.0: Required parameter: choose population unit.
        #================================================================================
        display("9.0 Applying RSF coefficients for population unit: " + popunit + ".\n"
                "            ...LOOPING...", detailsON, log)
        season = 1

        while season < 4:
            if popunit == "Yellowhead" and season == 1:
                c1 = -0.091;  c2 = -0.265;  c3 = -0.275;  c4 = 0.513; c5 = 0.503
                c6 = 0.342;   c7 = -0.008;  c8 = -0.001;  c9 = 0.001
                c10 = -0.00006 #divided by 100
                c11 = -1.950; c12 = 0.509;  c13 = -0.246
                c14 = -0.768; c15 = -4.130; c16 = -0.00572 #divided by 100
                # class breakpoints
                bp1 = 0.051496;  bp2 = 0.072563;  bp3 = 0.088948
                bp4 = 0.107674;  bp5 = 0.135763;  bp6 = 0.170874
                bp7 = 0.205985;  bp8 = 0.241096;  bp9 = 0.304296
                op = O.rsf_s1

            elif popunit == "Yellowhead" and season == 2:
                c1 = -0.111;  c2 = -0.214;  c3 = -0.530; c4 = -0.333;  c5 = 0.497
                c6 = -0.015;  c7 = -0.012;  c8 = -0.01;  c9 = -0.001
                c10 = -0.00029 #divided by 100
                c11 = -1.347; c12 = -0.008; c13 = -0.581
                c14 = -1.308; c15 = -4.104; c16 = -0.00703 #divided by 100
                # class breakpoints
                bp1 = 0.034278;  bp2 = 0.048969;  bp3 = 0.062027
                bp4 = 0.076718;  bp5 = 0.091408;  bp6 = 0.109363
                bp7 = 0.132215;  bp8 = 0.161597;  bp9 = 0.204036
                op = O.rsf_s2

            elif popunit == "Yellowhead" and season == 3:
                c1 = -0.025; c2 = -0.182; c3 = -0.421; c4 = -0.487; c5 = 0.398
                c6 = 0.054;  c7 = -0.011; c8 = -0.01;  c9 = -0.001
                c10 = -0.00055 #divided by 100
                c11 = -1.407; c12 = -0.008; c13 = -0.547
                c14 = -1.315; c15 = -4.171; c16 = -0.00717  #divided by 100
                # class breakpoints
                bp1 = 0.028076;  bp2 = 0.040554;  bp3 = 0.053032
                bp4 = 0.065510;  bp5 = 0.079548;  bp6 = 0.096706
                bp7 = 0.116983;  bp8 = 0.143499;  bp9 = 0.184053
                op = O.rsf_s3

            elif popunit == "Grande Cache" and season == 1:
                c1 = -1.468; c2 = -1.528; c3 = -1.622; c4 = -5.159; c5 = 1.224
                c6 = -3.955; c7 = -0.015; c8 = -0.002; c9 = -0.005
                c10 = 0.00005 #divided by 100
                c11 = -0.772; c12 = -0.619; c13 = -0.973
                c14 = 0.083;  c15 = -0.746; c16 = -0.00445 #divided by 100
                # class breakpoints
                bp1 = 0.011980;  bp2 = 0.035939;  bp3 = 0.044924
                bp4 = 0.059898;  bp5 = 0.077868;  bp6 = 0.092842
                bp7 = 0.107817;  bp8 = 0.131776;  bp9 = 0.161725
                op = O.rsf_s1

            elif popunit == "Grande Cache" and season == 2:
                c1 = -1.135; c2 = -2.188; c3 = -2.745; c4 = -5.762; c5 = 0.713
                c6 = 0.823;  c7 = -0.02;  c8 = -0.003; c9 = -0.014
                c10 = 0.00088 #divided by 100
                c11 = -1.629; c12 = -1.687; c13 = -0.286
                c14 = 1.035;  c15 = -2.809; c16 = -0.01261 #divided by 100
                # class breakpoints
                bp1 = 0.010309;  bp2 = 0.017181;  bp3 = 0.024053
                bp4 = 0.030926;  bp5 = 0.041234;  bp6 = 0.051543
                bp7 = 0.068724;  bp8 = 0.092777;  bp9 = 0.151192
                op = O.rsf_s2

            elif popunit == "Grande Cache" and season == 3:
                c1 = -0.148; c2 = -2.374; c3 = -3.342; c4 = -4.305; c5 = -1.864
                c6 = -3.468; c7 = -0.029; c8 = -0.018; c9 = -0.01
                c10 = 0.00208 #divided by 100
                c11 = -1.5;   c12 = -1.473;  c13 = -0.957
                c14 = -1.439; c15 = -1.69;   c16 = -0.00997 #divided by 100
                # class breakpoints
                bp1 = 0.017764;  bp2 = 0.031975;  bp3 = 0.042634
                bp4 = 0.056845;  bp5 = 0.071056;  bp6 = 0.08882
                bp7 = 0.110137;  bp8 = 0.142112;  bp9 = 0.206062
                op = O.rsf_s3

            elif popunit == "Clearwater" and season == 1:
                c1 = -10.574; c2 = -4.233; c3 = 4.11;   c4 = -10.125; c5 = 4.539
                c6 = 0.221;   c7 = -0.038; c8 = -0.137; c9 = 0.057
                c10 = 0.0006 #divided by 100
                c11 = -2.371; c12 = 0; c13 = -1.277
                c14 = -0.273; c15 = -2.784; c16 = -0.00042 #divided by 100
                bp1 = 0.505993;  bp2 = 0.606419;  bp3 = 0.641182
                bp4 = 0.795683;  bp5 = 0.880659;  bp6 = 0.915422
                bp7 = 0.938597;  bp8 = 0.961773;  bp9 = 0.977223
                op = O.rsf_s1

            elif popunit == "Clearwater" and season == 2:
                c1 = -1.616; c2 = -2.266; c3 = -3.148; c4 = -8.017; c5 = -2.523
                c6 = -2.283; c7 = -0.008; c8 = 0.01;   c9 = -0.024
                c10 = 0.0009 #divided by 100
                c11 = -1.942;  c12 = 0; c13 = -0.55
                c14 = -4.0853; c15 = -3.847; c16 = -0.00476 #divided by 100
                bp1 = 0.01589;  bp2 = 0.02270;  bp3 = 0.03178
                bp4 = 0.04086;  bp5 = 0.04994;  bp6 = 0.05902
                bp7 = 0.07037;  bp8 = 0.08172;  bp9 = 0.09761
                op = O.rsf_s2

            elif popunit == "Clearwater" and season == 3:
                c1 = -4.347; c2 = -2.061; c3 = -2.782; c4 = -2.139; c5 = -3.619
                c6 = -3.297; c7 = -0.014; c8 = 0.012;  c9 = -0.027
                c10 = 0.00043 #divided by 100
                c11 = -3.359; c12 = 0; c13 = -0.837
                c14 = -3.838; c15 = -2.707; c16 = -0.00009 #divided by 100
                bp1 = 0.00174;  bp2 = 0.00522;   bp3 = 0.008701
                bp4 = 0.013921; bp5 = 0.019141;  bp6 = 0.026102
                bp7 = 0.034803; bp8 = 0.045244;  bp9 = 0.071346
                op = O.rsf_s3

            elif popunit == "Livingstone-Castle" and season == 1:
                c1 = -10.0; c2 = -1.945; c3 = -0.651; c4 = -10; c5 = -0.133; c6 = -1.648
                c7 = -0.009; c8 = 0.035; c9 = -0.007; c10 = 0.00053 #divided by 100
                c11 = -0.087; c12 = 0; c13 = -1.566; c14 = -1.463; c15 = -1.304
                c16 = 0.00619 #divided by 100
                bp1 = 0.109723; bp2 = 0.152947; bp3 = 0.199496
                bp4 = 0.242720; bp5 = 0.275969; bp6 = 0.302569
                bp7 = 0.325843; bp8 = 0.355767; bp9 = 0.392342
                op = O.rsf_s1

            elif popunit == "Livingstone-Castle" and season == 2:
                c1 = -10.0; c2 = -5.792; c3 = -1.289; c4 = -10; c5 = -1.24; c6 = -1.382
                c7 = -0.016; c8 = 0.064;  c9 = -0.013; c10 = 0.0009 #divided by 100
                c11 = 0.345; c12 = 0; c13 = -2.301
                c14 = -1.042; c15 = -3.178; c16 = -0.00357  #divided by 100
                bp1 = 0.060329; bp2 = 0.097675; bp3 = 0.120658
                bp4 = 0.137894; bp5 = 0.155131; bp6 = 0.175241
                bp7 = 0.203969; bp8 = 0.238442; bp9 = 0.290153
                op = O.rsf_s2

            elif popunit == "Livingstone-Castle" and season == 3:
                c1 = -10.0;  c2 = -5.544; c3 = -1.064; c4 = -10; c5 = -0.916; c6 = 1.127
                c7 = -0.009; c8 = 0.004;  c9 = -0.015; c10 = -0.00008 #divided by 100
                c11 = -0.875; c12 = 0; c13 = -3.518; c14 = 2.385; c15 = -4.005
                c16 = -0.01061  #divided by 100
                bp1 = 0.023971; bp2 = 0.038952; bp3 = 0.047941
                bp4 = 0.05693;  bp5 = 0.06916;  bp6 = 0.083897
                bp7 = 0.107868; bp8 = 0.137831; bp9 = 0.182776
                op = O.rsf_s3

            elif popunit == "Swan Hills" and season == 1:
                c1 = -1.279; c2 = 3.455; c3 = -1.098; c4 = -1.807; c5 = 4.1; c6 = 2.682
                c7 = -0.008; c8 = -0.05; c9 = -0.006; c10 = 0.00083 #divided by 100
                c11 = -2.069; c12 = 0.473; c13 = -4.601
                c14 = -5.81;  c15 = 2.406; c16 = -0.00381  #divided by 100
                bp1 = 0.020635;  bp2 = 0.038691;  bp3 = 0.056747
                bp4 = 0.072224;  bp5 = 0.087700;  bp6 = 0.105756
                bp7 = 0.134130;  bp8 = 0.196036;  bp9 = 0.286315
                op = O.rsf_s1

            elif popunit == "Swan Hills" and season == 2:
                c1 = -1.123; c2 = -1.439; c3 = -2.948; c4 = -3.415; c5 = -1.708
                c6 = -3.141; c7 = -0.027; c8 = -0.017; c9 = -0.014
                c10 = 0.00104 #divided by 100
                c11 = -1.058; c12 = 0.18;   c13 = -2.655
                c14 = -0.272; c15 = -2.748; c16 = -0.02172  #divided by 100
                bp1 = 0.003485; bp2 = 0.006971; bp3 = 0.010456
                bp4 = 0.013941; bp5 = 0.017436; bp6 = 0.022654
                bp7 = 0.029625; bp8 = 0.040081; bp9 = 0.057507
                op = O.rsf_s2

            elif popunit == "Swan Hills" and season == 3:
                c1 = -0.177; c2 = -0.393; c3 = -1.307; c4 = -1.594; c5 = -1.643
                c6 = -1.669; c7 = -0.034; c8 = -0.023; c9 = 0.0001
                c10 = -0.00014 #divided by 100
                c11 = -1.59;  c12 = -1.635; c13 = 0.653
                c14 = -1.021; c15 = 1.503;  c16 = -0.00864 #divided by 100
                bp1 = 0.002367; bp2 = 0.004734; bp3 = 0.007101
                bp4 = 0.010651; bp5 = 0.014202; bp6 = 0.020119
                bp7 = 0.028404; bp8 = 0.040238; bp9 = 0.061541
                op = O.rsf_s3

             # Chinchaga models are were not devloped using GPS data and are unreliable 
##            elif popunit == "Chinchaga" and season == 1:
##                c1 = -1.468; c2 = -1.528; c3 = -1.622; c4 = -5.159; c5 = 1.224; c6 = -3.955
##                c7 = -0.015; c8 = -0.002; c9 = -0.005; c10 = 0.005 #divided by 100
##                c11 = -0.772; c12 = -0.619; c13 = -0.973
##                c14 = 0.083;  c15 = -0.746; c16 = -0.445  #divided by 100
##                bp1 = 0.015422;  bp2 = 0.073998;  bp3 = 0.101745
##                bp4 = 0.120242;  bp5 = 0.135657;  bp6 = 0.154155
##                bp7 = 0.175736;  bp8 = 0.206565;  bp9 = 0.243561
##                op = O.rsf_s1
##
##            elif popunit == "Chinchaga" and season == 2:
##                c1 = -1.135; c2 = -2.188; c3 = -2.745; c4 = -5.762; c5 = 0.713
##                c6 = 0.823; c7 = -0.02; c8 = -0.003; c9 = -0.014
##                c10 = 0.088 #divided by 100
##                c11 = -1.629; c12 = -1.687;   c13 = -0.286
##                c14 = 1.035; c15 = -2.809; c16 = -1.261  #divided by 100
##                bp1 = 0.003524; bp2 = 0.024641; bp3 = 0.049278
##                bp4 = 0.070395; bp5 = 0.091513; bp6 = 0.11615
##                bp7 = 0.144306; bp8 = 0.183021; bp9 = 0.242854
##                op = O.rsf_s2
##
##            elif popunit == "Chinchaga" and season == 3:
##                c1 = -0.148; c2 = -2.374; c3 = -3.342; c4 = -4.305; c5 = -1.864
##                c6 = -3.468; c7 = -0.029; c8 = -0.018; c9 = -0.01
##                c10 = 0.208 #divided by 100
##                c11 = -1.5;  c12 = -1.473; c13 = -0.957
##                c14 = -1.439; c15 = -1.69;  c16 = -0.997 #divided by 100
##                bp1 = 0.018591; bp2 = 0.074226; bp3 = 0.118734
##                bp4 = 0.155824; bp5 = 0.200332; bp6 = 0.263385
##                bp7 = 0.348691; bp8 = 0.46367; bp9 = 0.593485
##                op = O.rsf_s3

            #----------------------------------------------------------------------------
            # 9.1-3 Calculation of RSF. RSF creation process goes here, within the while
            #       loop
            #============================================================================
            display("9." + str(season) + " Calculating RSF for SEASON " + str(season) +
                    " and " + popunit + " POPULATION UNIT.", detailsON, log)
            tmpRSF = ((Raster(TR.wtree) * (c1)) + (Raster(TR.rgn) * (c2)) +
                      (Raster(TR.shrub) * (c3)) + (Raster(TR.wherb) * (c4)) +
                      (Raster(TR.uherb) * (c5)) + (Raster(TR.noveg) * (c6)) +
                      (Raster(TR.cc_trd) * (c7)) + (Raster(TR.cc_rgn) * (c8)) +
                      (Raster(TR.sc_x_utree) * (c9)) + (Raster(R.cti150mi) * (c10)) +
                      (Raster(TR.dfor_utree) * (c11)) + (Raster(TR.dfor_wtree) * (c12)) +
                      (Raster(TR.dfor_uhrb) * (c13)) + (Raster(TR.dfor_ergn) * (c14)) +
                      (Raster(TR.dfor_envg) * (c15)) + (Raster(R.d500_strm) * (c16)))
            tmpRSF.save(TR.rsf_raw)

            #----------------------------------------------------------------------------
            # Rescale to between 0 and 1,the rescaled RSF is exp(lp)/(1 + exp(lp))
            #============================================================================
            display("   ...Rescaling RSF between 0 and 1.", detailsON, log)
            rsTmpRSF = ((Exp(Raster(TR.rsf_raw)))/(1 + (Exp(Raster(TR.rsf_raw)))))
            rsTmpRSF.save(TR.rsf_rscld)

            #----------------------------------------------------------------------------
            # Apply masks: Agriculture, water/snow/ice, and nonveg within  alpine
            #               subregion. Then, apply bear range scalar; then removals
            #               (habitat deletion)
            #============================================================================
            display("   ...Applying female range scalar and non-habitat mask.",detailsON,
                    log)
            rngRSF = (Raster(TR.rsf_rscld) * Raster(R.rng_gbz) * Raster(R.mask) * 0.01)
            rngRSF.save(TR.rsf_rng)

            #----------------------------------------------------------------------------
            # Classify seasonal RSF output into 10 classes
            #============================================================================
            display("   ...Reclassifying rescaled RSF to 10 classes.", detailsON, log)
            remap = RemapRange([[0,bp1,1], [bp1,bp2,2], [bp2,bp3,3], [bp3,bp4,4],
                                [bp4,bp5,5], [bp5,bp6,6], [bp6,bp7,7], [bp7,bp8,8],
                                [bp8,bp9,9], [bp9,1,10]])
            sRSF = Reclassify(TR.rsf_rng, "VALUE", remap, "DATA")
            sRSF.save(op)
            display("      -> Saved output: " + op, detailsON, log)

            # increment season and LOOP
            season = season + 1
            
        display("            ...DONE LOOPING...", detailsON, log)
        
        #--------------------------------------------------------------------------------
        # 9.4: Calculate max rsf.
        #================================================================================
        display("9.4 Calculating maximum RSF.", detailsON, log)
        mRSF = CellStatistics([O.rsf_s1, O.rsf_s2, O.rsf_s3], "MAXIMUM")
        mRSF.save(O.rsf_max)

        #--------------------------------------------------------------------------------
        # 9.5: OPTIONAL: Apply habitat deletions (mines, urban, wellsites, etc) to model
        #================================================================================
        display("9.5 Optional Parameter: Habitat Deletions.", detailsON, log)

        if habDel != "#":
            Pdisplay("      ...Habitat deletions will be represented by: \n" +
                     "      ..." + habDel, log)

            # Copy users input so it's not modified
            display("   9.5.0 Copying input deletion layer.", detailsON, log)
            arcpy.FeatureClassToFeatureClass_conversion(habDel, tmpgdb, "removals")

            display("   9.5.1 Adding GRIDCODE field.", detailsON, log)
            # See if there is a GRIDCODE field in AOI, if so delete first
            fields = arcpy.ListFields(tFC.removals)
            for field in fields:
                if field.name == "GRIDCODE":
                    arcpy.DeleteField_management(tFC.removals, "GRIDCODE")

            # Add GRIDCODE field and set value to 1
            arcpy.AddField_management(tFC.removals, "GRIDCODE", "short")
            arcpy.CalculateField_management(tFC.removals, "GRIDCODE", 1)

            display("   9.5.2 Adding new removals to AOI.", detailsON, log)
            # replaces existing one (copy of AOI). Will output a shapefile but not a fc
            arcpy.Union_analysis([S.AOI, tFC.removals], S.update, "ALL", "", "GAPS")

            display("   9.5.3 Converting opening data to raster.", detailsON, log)
            arcpy.PolygonToRaster_conversion(S.update, "GRIDCODE_1", TR.deletion, 
                                             "MAXIMUM_COMBINED_AREA")

            display("   9.5.4 Creating position raster.", detailsON, log)
            pOpen = Reclassify(TR.deletion,"VALUE", RemapValue([[1,2],[0,1]]), "DATA")
            pOpen.save(TR.order)

            display("   9.5.5 Updating RSF raster.", detailsON, log)
            uRSF = Pick(TR.order, [O.rsf_max, 0])
            uRSF.save(TR.rsfupd)

            display("   9.5.6 Replacing existing RSF raster with updated raster.",
                    detailsON, log)
            arcpy.CopyRaster_management(TR.rsfupd, O.rsf_max)

        else:
            Pdisplay("   ...No habitat deletions <section skipped>", log)

        #--------------------------------------------------------------------------------
        # 9.6: Calculate mean rsf_max.
        #================================================================================
        display("9.6 Calculating mean rsf_max.", detailsON, log)
        mRSFm = ZonalStatisticsAsTable(TR.bndgrid, "Value", O.rsf_max, O.rsf_max_tbl,
                                       "DATA", "MEAN")

        #--------------------------------------------------------------------------------
        # 9.7: Reclassify RSF to non-critical/secondary/primary.
        #================================================================================
        display("9.7 Separating RSF (10 classes) into non-critical, secondary & primary",
                detailsON, log)
        remap = RemapRange([[1,4,0],[5,7,1],[8,10,2]])
        rcRSF = Reclassify(O.rsf_max, "VALUE", remap, "DATA")
        rcRSF.save(O.rsf_class)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in calRSF():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 10 Calculate Risk                    
#========================================================================================
def calRisk():
    try:
        #--------------------------------------------------------------------------------
        # 10.0: Apply coefficients and calculate raw scores.
        #================================================================================
        display("10.0 Applying coefficients and calculating raw scores.", detailsON, log)
        c1 = 1.862
        c2 = 0.01352 # divided by 100
        c3 = 3.573
        c4 = -0.00064 # divided by 100
        c6 = -0.02224
        c7 = 0.03909
        c8 = -0.03090
        c9 = -4.865

        tmpRisk = ((Raster(TR.ned_rd) * (c1)) + (Raster(R.dstrm) * (c2)) +
                   (Raster(TR.dfor_rscld) * (c3)) + (Raster(R.tri) * (c4)) +
                   (R.p_uptree) + (Raster(R.protctd6mi) * (c6)) +
                   (Raster(TR.ned_tr500) * Raster(R.protctd6mi) * (c7)) +
                   (Raster(TR.ned_tr100) * Raster(R.whitez_6mi) * (c8)) + c9)
        tmpRisk.save(TR.risk_raw)

        #--------------------------------------------------------------------------------
        # 10.1: Rescale RISK raster between 0 and 1.
        #================================================================================
        display("10.1 Rescaling RISK raster between 0 and 1.", detailsON, log)
        rsRisk = ((Exp(Raster(TR.risk_raw)))/(1 + (Exp(Raster(TR.risk_raw)))))
        rsRisk.save(TR.risk_rscld)

        #--------------------------------------------------------------------------------
        # 10.2: Reclassify RISK raster to 10 classes.
        #================================================================================
        display("10.2 Reclassifying risk to 10 classes.", detailsON, log)
        rcRisk = Reclassify(TR.risk_rscld, "VALUE",
                            RemapRange([[0,0.01,0],[0.01,0.1,1],[0.1,0.2,2],[0.2,0.3,3],
                                        [0.3,0.4,4],[0.4,0.5,5],[0.5,0.6,6],[0.6,0.7,7],
                                        [0.7,0.8,8],[0.8,0.9,9],[0.9,1,10]]), "DATA")
        rcRisk.save(O.risk)

        #--------------------------------------------------------------------------------
        # 10.3: Calculate mean risk.
        #================================================================================
        display("10.3 Calculating mean Risk.", detailsON, log)
        mRisk = ZonalStatisticsAsTable(TR.bndgrid, "Value", O.risk, O.risk_tbl,
                                       "DATA", "MEAN")

        #--------------------------------------------------------------------------------
        # 10.4: Reclassify RISK (10 classes) to sink/source.
        #================================================================================
        display("10.4 Separating risk (10 classes) into sink or source.", detailsON, log)
        rcSSRick = Reclassify(O.risk, "VALUE", RemapRange([[0,5,1],[6,10,-1]]), "DATA")
        rcSSRick.save(O.risk_class)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in calRisk():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise

#----------------------------------------------------------------------------------------
# 11 Calculate Habitat States
#========================================================================================
def calHState():
    try:
        #--------------------------------------------------------------------------------
        # 11.0: Combine RISK and RSF to calculate habitat states.
        #================================================================================
        display("11.0 Combining RISK and RSF to create habitat states.", detailsON, log)
        # Output raster name (based on user input)
        hState = Raster(O.risk_class) * Raster(O.rsf_class)
        hState.save(O.habState)

        #--------------------------------------------------------------------------------
        # 11.1: Calculate mean of habitat states.
        #================================================================================
        display("11.1 Creating output habitat state table.", detailsON, log)
        mHState = ZonalStatisticsAsTable(TR.bndgrid, "Value", O.habState, O.habState_tbl,
                                         "DATA", "MEAN")

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in calHState():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise


#========================================================================================
# Buffer input features, add grid code and export back to raster.
# inLines -> user input; friLines -> base roads or trails; cv -> cell value
#========================================================================================
def prepLines(step, inLines, newBuf, friLines, cv):
    try:
        #--------------------------------------------------------------------------------
        # Buffer line features by 30m. (60m in total)
        #================================================================================
        display(step + ".0 Buffering line features by 30m.", detailsON, log)
        arcpy.Buffer_analysis(inLines, newBuf, "30", "", "", "ALL", "")

        #--------------------------------------------------------------------------------
        # Create new field to populate.
        #================================================================================
        display(step + ".1 Adding GRIDCODE field to layer.", detailsON, log)
        # See if there is a GRIDCODE field in the input layer, if so delete first
        fields = arcpy.ListFields(newBuf)
        for field in fields:
            if field.name == "GRIDCODE":
                arcpy.DeleteField_management(newBuf, "GRIDCODE")

        # Add GRIDCODE field and set value
        arcpy.AddField_management(newBuf, "GRIDCODE", "short")
        arcpy.CalculateField_management(newBuf, "GRIDCODE", cv)

        #--------------------------------------------------------------------------------
        # Append new features to existing trails OR roads 
        #================================================================================
        if process == 1 or process == 3:
            display(step +".2 Appending new features to existing line features.",
                    detailsON, log)
            arcpy.Append_management(inLines, friLines, "NO_TEST")

        #--------------------------------------------------------------------------------
        # Union updated features with boundary.  The field [GRIDCODE] has a value of
        #      3 to indicate it is a Right-of-Way (RoW), a value of 0 everywhere else.
        #================================================================================
        display(step + ".3 Unioning new and existing features with boundary.", detailsON,
                log)
        arcpy.Union_analysis([S.AOI, newBuf] , S.update, "ALL", "", "GAPS")

        #--------------------------------------------------------------------------------
        # Convert updated features to raster.  A value of 3 to indicates it is a
        #       Right-of-Way (RoW), and a value of 0 everywhere else.
        #================================================================================
        display(step + ".4 Converting unioned features to raster.", detailsON, log)
        arcpy.PolygonToRaster_conversion(S.update, "GRIDCODE_1", TR.change,
                                         "MAXIMUM_COMBINED_AREA")

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in prepLines():"
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise


#========================================================================================
# Create distance to forest edge with landcover & regen inputs (within addClearing()
# function)
#========================================================================================
def distFEdge(ipLC, opCont, aBorder, opEuc, opMA):
    try:
        display("             ...creating contours.", detailsON, log)
        # Create contour with Interval = 1, base = 0.5, Z Factor = 1
        Contour(ipLC, opCont, "1", "0.5", "1")

        display("             ...appending to buffered edge.", detailsON, log)
        # Append contour output to the buffered border
        arcpy.Append_management(opCont, aBorder, "NO_TEST")

        display("             ...calculating distance to edge.", detailsON, log)
        # Calculate distance to edge of landcover
        dist = EucDistance(aBorder, "", "30", "")
        dist.save(opEuc)

        display("             ...rescaling.", detailsON, log)
        if ipLC == TR.forest:
            # Scale between 0 and 1
            calScale = (Exp(-0.001 * (Raster(opEuc))))
        else:
            # Scale between 0 and 1
            calScale = ((1 - Exp((-1 * (Raster(opEuc)/500))))) * Raster(ipLC)
        calScale.save(opMA)

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in distFEdge():"
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]), log)
        raise


#========================================================================================
# Create layer files
#========================================================================================
def createLayerFiles():
    try:
        lblLayer = ["Habitat State", "RSF Group", "Risk", "Risk Class"]
        lyrList = ["Habitat State (current).lyr", "RSF Results.lyr", "Risk.lyr",
                   "Risk Class.lyr"]
        rstList = ["hab_state", "rsf", "risk", "risk_class"]
        
        if GBTpath == "X:\\gbp\\Deliverables\\GBTools":
            LYRpath = "X:\\gbp\\Deliverables\\_Common_GBTool_Inputs\\_LayerFiles"
        else:
            LYRpath = fp(GBTpath, "layerfiles")
            
        TARGETpath = fp(GBTpath, "results", "HAB_STATE", sys.argv[2])

        if process == 1: lyrRange = [0,1,2,3]
        elif process == 2: lyrRange = [1]
        elif process == 3: lyrRange = [2,3]
        
        for i in lyrRange:
            display("       creating " + lblLayer[i] + " layer", detailsON, log)
            HSlyr = fp(LYRpath, lyrList[i])
            lyr = arcpy.mapping.Layer(HSlyr)

            if i == 1:
                 for layer in arcpy.mapping.ListLayers(lyr)[0]:
                    oldPath = layer.workspacePath
                    layer.findAndReplaceWorkspacePath(oldPath, TARGETpath)
            else:
                lyr.replaceDataSource(TARGETpath, 'RASTER_WORKSPACE', rstList[i])

            lyr.name = lblLayer[i] + ": " + sys.argv[2]
            lyr.saveACopy(fp(TARGETpath, lyrList[i]))

    except:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        error("\nPYTHON ERRORS in createLayerFiles():"
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
             f.sep25 + "   TOOL: HABITAT STATES" +
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
