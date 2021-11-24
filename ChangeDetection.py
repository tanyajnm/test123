#----------------------------------------------------------------------------------------
#                                      fRI RESEARCH
#                                  GRIZZLY BEAR PROGRAM
#----------------------------------------------------------------------------------------
# SCRIPT NAME:  ChangeDetection.py
#               v.2019.0705
#========================================================================================
# Developed by:  Dan Wismer - GIS Analyst at fRI Research
#
# Date created:  November 9th, 2018
#
# Purpose:       Generates repstive cell value area km2 grouped bar plots of a reference
#               and scenario habitat states rasters and calculates the percent
#               change. Change detection rasters are also generated
#
# Arguments:    [1] Reference habitat states folder containing hab_state, risk,
#               rsf_max, rsf_sa, rsf_s2 and rsf_s3 GRIDS
#
#              [2] Scenario habitat states folder containing hab_state, risk,
#               rsf_max, rsf_s1, rsf_s2 and rsf_s3 GRIDS
#
#              [3] Multi-value check box selection of which rasters to compare
#                  ie. hab_state, risk, rsf_max, rsf_s1, rsf_s2, rsf_s3 
#
# Output:       .png bar plot(s), change detection rasters/layers
#
# Edits:        June 2019 (J.Duval) - Replaced '//' with '\\'
#               July 2019 (J.Duval) - Paths concatenated with os.path.join() as fp
#
#----------------------------------------------------------------------------------------
# Pre-Steps: Get user input and set up tool
#========================================================================================
#Import libraries
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import matplotlib.pyplot as plt
import sys, subprocess
import os
import shutil
import gc
import traceback
from os.path import join as fp
from CheckLicenses import CheckArcInfo, CheckSpatialExt

#----------------------------------------------------------------------------------------
# Only proceed if the required ArcGIS licenses are available
#========================================================================================
if CheckArcInfo() == "yes" and CheckSpatialExt() == "yes":
    pass
else:
    arcpy.AddWarning("!!!ArcGIS Advanced and Spatial Analyst are required")
    sys.exit()

#set-up tool environments
gc.collect()  
arcpy.env.overwriteOutput = True 

#Read-in user input 
refHabStates_folder = arcpy.GetParameterAsText(0)
refUserName = arcpy.GetParameterAsText(1)
scnroHabStates_folder = arcpy.GetParameterAsText(2)
scnroUserName = arcpy.GetParameterAsText(3)
userHabStates = arcpy.GetParameterAsText(4)

#Get script location and set mxd and layer templates
GBTpath = os.path.split(sys.path[0])[0]
MxdPath = fp(sys.path[0], "_MxdTemplates")

mxd_template = fp(MxdPath, "Change_Detection_Template.mxd")
colour_ramp_hab_rsf = fp(MxdPath, "Colour_Ramp_Hab_State_Rsf.jpg")
colour_ramp_risk = fp(MxdPath, "Colour_Ramp_Risk.jpg")

#fRI Research internal paths
if GBTpath == "X:\\gbp\\Deliverables\\GBTools":
    Hab_Lyr = fp(GBTpath, "_Common_GBTool_Inputs", "_LayerFiles",
                           "Change_Detection_Hab_State.lyr")
    Rsf_Lyr = fp(GBTpath, "_Common_GBTool_Inputs", "_LayerFiles",
                           "Change_Detection_Rsf.lyr")
    Risk_Lyr = fp(GBTpath, "_Common_GBTool_Inputs", "_LayerFiles",
                            "Change_Detection_Risk.lyr")

#Get paths when tool is bundled for distribution (Deliv)
else:
    Hab_Lyr =  fp(GBTpath, "LayerFiles", "Change_Detection_Hab_State.lyr")
    Rsf_Lyr =  fp(GBTpath, "LayerFiles", "Change_Detection_Rsf.lyr")
    Risk_Lyr = fp(GBTpath, "LayerFiles", "Change_Detection_Risk.lyr")
    

#Clean nested scratch folder if it exists
try:
    path = fp(scnroHabStates_folder, "Scratch")
    if arcpy.Exists(path):
        shutil.rmtree(path)
except:
    pass

#Create new nested plot output, change detection and scratch folders
arcpy.CreateFolder_management (scnroHabStates_folder, "PctPlots")
plotFolder = fp(scnroHabStates_folder, "PctPlots")
arcpy.CreateFolder_management (scnroHabStates_folder, "ChangeDetection")
ChangeDetectionFolder = fp(scnroHabStates_folder, "ChangeDetection")
#
arcpy.CreateFolder_management (scnroHabStates_folder, "Scratch")
scratch = fp(scnroHabStates_folder, "Scratch")
arcpy.CreateFolder_management (scratch, "Reference")
arcpy.CreateFolder_management (scratch, "Scenario")

#
#----------------------------------------------------------------------------------------
#Step 1: Generate string list with GRID names, referencing user habitat selection
#========================================================================================
arcpy.AddMessage("Step 1:   Verify user selected habitat states GRID names")
userHabStates = userHabStates.split(";")
HabitatRasterNames = []
namesForPlot = []
arcpy.AddMessage("...       User selected the following habitat states:")
for HabState in userHabStates:
    if HabState == "'Habitat State'":
        arcpy.AddMessage("...       Habitat State")
        name = "hab_state"
        HabitatRasterNames.append(name)
        namesForPlot.append("Habitat States")
#
    elif HabState == "Risk":
        arcpy.AddMessage("...       Risk")
        name = "risk"
        HabitatRasterNames.append(name)
        namesForPlot.append("Risk")
#
    elif HabState == "'RSF Max'":
        arcpy.AddMessage("...       RSF Max")
        name = "rsf_max"
        HabitatRasterNames.append(name)
        namesForPlot.append("RSF Max")
#        
    elif HabState == "'RSF Season 1'":
        arcpy.AddMessage("...       RSF Season 1")
        name = "rsf_s1"
        HabitatRasterNames.append(name)
        namesForPlot.append("RSF Season 1")
#        
    elif HabState == "'RSF Season 2'":
        arcpy.AddMessage("...       RSF Season 2")
        name = "rsf_s2"
        HabitatRasterNames.append(name)
        namesForPlot.append("RSF Season 2")
#           
    elif HabState == "'RSF Season 3'":
        arcpy.AddMessage("...       RSF Season 3")
        name = "rsf_s3"
        HabitatRasterNames.append(name)
        namesForPlot.append("RSF Season 3")
#
#----------------------------------------------------------------------------------------
#Step 2: Read-in habitat state reference GRIDS
arcpy.AddMessage(" ")
arcpy.AddMessage("Step 2:   Read-in habitat states reference GRIDS")
arcpy.AddMessage(" ")
arcpy.env.workspace = refHabStates_folder
#========================================================================================
refRasters = []
#raster names: hab_state, risk, rsf_max, rsf_s1, rsf_s2, rsf_s3
for name in HabitatRasterNames:
    scnro_rst = fp(refHabStates_folder, name)
    if arcpy.Exists(scnro_rst):
        arcpy.AddMessage("...       Found '" + name + "' reference")
        refRasters.append(scnro_rst)
    else:
        arcpy.AddError("...       No '" + name + "' raster GRID found in "
                                  + refHabStates_folder + " folder. " +
                                  "Check naming convention to run this tool")
        sys.exit()
            
#----------------------------------------------------------------------------------------
#Step 3: Read-in habitat state scenario GRIDS
arcpy.AddMessage(" ")
arcpy.AddMessage("Step 2:   Read-in habitat states scenario GRIDS")
arcpy.env.workspace = scnroHabStates_folder
#========================================================================================
scnroRasters = []
#raster names: hab_state, risk, rsf_max, rsf_s1, rsf_s2, rsf_s3
for name in HabitatRasterNames:
    scnro_rst = fp(scnroHabStates_folder, name)
    if arcpy.Exists(scnro_rst):
        arcpy.AddMessage("...       Found '" + name + "' scenario")
        scnroRasters.append(scnro_rst)
    else:
        arcpy.AddError("...       No '" + name + "' raster GRID found in " +
                                  scnroHabStates_folder + " folder. " +
                                  "Check naming convention to run this tool")
        sys.exit()

#----------------------------------------------------------------------------------------
#Step 4: Compare reference and scenario raster properties
arcpy.AddMessage(" ")
arcpy.AddMessage("Step 4:   Check GRID properties")
#========================================================================================
maskRasters = []
for ref, scnro in zip(refRasters, scnroRasters):
    #Get/check reference crs
    name = arcpy.Describe(ref).name
    ref_name = arcpy.Describe(ref).name
    ref_path = arcpy.Describe(ref).path
    ref_crs = arcpy.Describe(fp(ref_path, ref_name)).spatialReference
    if str(ref_crs.type) == "Projected":
        arcpy.AddMessage("...       Reference " + ref_name +" is projected:")
        arcpy.AddMessage("...       CRS: " + str(ref_crs.Name))
    else:
        arcpy.AddMessage("...       Reference " + ref_name +" is not projected")
        arcpy.AddMessage("...       Verify that rasters have been projected")
        sys.exit()
#----------------------------------------------------------------------------------------
    #Get/check scenario crs
    scnro_name = arcpy.Describe(scnro).name
    scnro_path = arcpy.Describe(scnro).path
    scnro_crs = arcpy.Describe(fp(scnro_path, scnro_name)).spatialReference
    if str(scnro_crs.type) == "Projected":
        arcpy.AddMessage("...       Scenario " + scnro_name +" is projected:")
        arcpy.AddMessage("...       CRS: " + str(scnro_crs.Name))
    else:
        arcpy.AddMessage("...       Scenario " + scnro_name +" is not projected")
        arcpy.AddMessage("...       Verifiy rasters have been projected")
        sys.exit()
#----------------------------------------------------------------------------------------
    refCRS = arcpy.Describe(ref).spatialReference.alias
    scnroCRS = arcpy.Describe(scnro).spatialReference.alias
    #Check coordinate system
    if refCRS != scnroCRS:
        arcpy.AddMessage("...       " + name + "coordinate reference system failed check")
        arcpy.AddMessage("""...       Check that all reference and scenario rasters are in
                                      the same projection""")
        sys.exit()
    else:    
        arcpy.AddMessage("...       Both " + name + " CRS validated")
#----------------------------------------------------------------------------------------
    #Check cell size
    refCellSizeX = arcpy.GetRasterProperties_management(ref, "CELLSIZEX")    
    refCellSizeY = arcpy.GetRasterProperties_management(ref, "CELLSIZEY")
    refCellSizeTotal = int(refCellSizeX.getOutput(0)) + int(refCellSizeY.getOutput(0))
#
    scnroCellSizeX = arcpy.GetRasterProperties_management(scnro, "CELLSIZEX")
    scnroCellSizeY = arcpy.GetRasterProperties_management(scnro, "CELLSIZEY")
    scnroCellSizeTotal = int(scnroCellSizeX.getOutput(0)) + int(scnroCellSizeY.getOutput(0))
#
    if refCellSizeTotal != scnroCellSizeTotal:
        arcpy.AddMessage("...       " + name + " cell size failed check")
        arcpy.AddMessage("""...       Verify that all reference and scenario rasters have
                                     the same cell size""")
        sys.exit()
    else:
        arcpy.AddMessage("...       " + name + " cell validated: " + str(refCellSizeX) + "m")
#----------------------------------------------------------------------------------------
    #Get reference hab state rows/columns
    refRow = arcpy.GetRasterProperties_management(ref, "ROWCOUNT")
    refCol = arcpy.GetRasterProperties_management(ref, "COLUMNCOUNT")
    refRowCol = int(refRow.getOutput(0)) * int(refCol.getOutput(0))

    scnroRow = arcpy.GetRasterProperties_management(scnro, "ROWCOUNT")
    scnroCol = arcpy.GetRasterProperties_management(scnro, "COLUMNCOUNT")
    scnroRowCol = int(scnroRow.getOutput(0)) * int(scnroCol.getOutput(0))
#
    if scnroRowCol < refRowCol:
        arcpy.AddMessage("...       Scenario " + name + " will be used as the mask.")
        maskRasters.append(scnro)
#
    elif scnroRowCol > refRowCol:
        arcpy.AddMessage("...       Reference " + name + " will be used as the mask.")
        maskRasters.append(ref)
#
    else:
        arcpy.AddMessage("...       Scenario " + name + " will be used as the mask.")
        maskRasters.append(ref)

    arcpy.AddMessage(" ")    
#----------------------------------------------------------------------------------------
#Step 5: Align reference and scenario rasters
arcpy.AddMessage("Step 5:   Align reference and scenario GRIDS:")
arcpy.AddMessage(" ")
#========================================================================================
for mask, ref, scnro in zip(maskRasters, refRasters, scnroRasters):
    try:
        name = arcpy.Describe(ref).name
        arcpy.env.snapRaster = mask
        arcpy.env.mask = mask
        arcpy.AddMessage("...       Aligning reference " + name)
        refMask = ExtractByMask(ref, mask)
        refMask.save(fp(scratch, "Reference", name))

        arcpy.AddMessage("...       Aligning scenario " + name )
        scnroMask = ExtractByMask(scnro, mask)
        scnroMask.save(fp(scratch, "Scenario", name))    
        arcpy.AddMessage(" ")
##    except:
##        arcpy.AddError("...       Reference / scenario rasters do not overlap")
##        arcpy.AddError("...       OR, spatial analyst extension is not available")
##        sys.exit()
    except arcpy.ExecuteError:
        # Get the geoprocessing error messages
        arcpy.AddError("\nGP ERRORS:\n" + 
              arcpy.GetMessage(0) + "\n" +
              arcpy.GetMessages(2))
    except:
        # sys.exc_info() provides a tuple (type, value, traceback)
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError("\nPYTHON ERRORS in main():" 
              "\nTraceback Info:    " + tbinfo + 
              "\nError Type:    " + str(sys.exc_info()[0]) + 
              "\nError Info:    " + str(sys.exc_info()[1]))#
        sys.exit()
#----------------------------------------------------------------------------------------
#Step 6: Generate change raster
arcpy.AddMessage(" ")
arcpy.AddMessage("Step 6:   Generating %Change Report:")
arcpy.AddMessage(" ")
#========================================================================================
#Get list of aligned reference rasters
wksp = fp(scratch, "Reference")
arcpy.env.workspace = wksp
refRasterList = arcpy.ListRasters("*", "GRID")
refRasterList = [fp(wksp, raster) for raster in refRasterList]

#Get list of aligned scenario rasters
wksp = fp(scratch, "Scenario")
arcpy.env.workspace = wksp
scnroRasterList = arcpy.ListRasters("*", "GRID")
scnroRasterList = [fp(wksp, raster) for raster in scnroRasterList]

#Get template map document to populate new raster layers
mxd = arcpy.mapping.MapDocument(mxd_template)  
df = arcpy.mapping.ListDataFrames(mxd, "Layers")[0]
pngPaths = []

for name, ref, scnro in zip(namesForPlot, refRasterList, scnroRasterList):
#.........................Generate change detection rasters..............................
    #Turn off all layers
    layers = arcpy.mapping.ListLayers(mxd, "*", df)
    for layer in layers:
        layer.visible = False
#......................Habitat States....................................................
    arcpy.AddMessage("...       " + name)
    raster_name = arcpy.Describe(ref).name
    if name == "Risk":
        changeRaster = (Raster(ref) - Raster(scnro))
    else:    
        changeRaster = (Raster(ref) - Raster(scnro)) * -1

    #Convert 0 change to NoData
    changeRaster = SetNull(changeRaster, changeRaster, "VALUE = 0")
    cd = fp(ChangeDetectionFolder, raster_name + "_dif")
    changeRaster.save(cd)

    #Get area of change detection
    ChangeDetectionArea = []
    CellSize = arcpy.GetRasterProperties_management(cd, "CELLSIZEX")
    CellSize_x2 = float(CellSize.getOutput(0)) * float(CellSize.getOutput(0))
    try:
        with arcpy.da.SearchCursor(cd, ["VALUE", "COUNT"]) as cursor:
            for row in cursor:
                sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                ChangeDetectionArea.append(sum_area)
        del row
        del cursor
    except:
        ChangeDetectionArea.append(0)
        
    #Create .lyr file and map
    addLayer = arcpy.mapping.Layer(fp(ChangeDetectionFolder, raster_name + "_dif"))
    arcpy.mapping.AddLayer(df,addLayer, "TOP")
    updateLayer = arcpy.mapping.ListLayers(mxd, raster_name + "_dif", df)[0]

    #Use template symbology and map elements for habitat state, rsf and risk
    if name == "Habitat States": sourceLayer = arcpy.mapping.Layer(Hab_Lyr)
    elif name == "Risk": sourceLayer = arcpy.mapping.Layer(Risk_Lyr)
    else: sourceLayer = arcpy.mapping.Layer(Rsf_Lyr)

    arcpy.mapping.UpdateLayer(df, updateLayer, sourceLayer, True)
    updateLayer.name = "Change Detection: " + name
    arcpy.SaveToLayerFile_management(updateLayer,
                        fp(ChangeDetectionFolder, raster_name + ".lyr"))

    for elm in arcpy.mapping.ListLayoutElements(mxd, "PICTURE_ELEMENT"):
        if name == "Habitat States":  elm.sourceImage = colour_ramp_hab_rsf
        elif name == "Risk": elm.sourceImage = colour_ramp_risk
        else: elm.sourceImage = colour_ramp_hab_rsf

    arcpy.RefreshActiveView()
    arcpy.RefreshTOC()
            
#...............................Generate plots...........................................
#Get name and path of reference and scenario raster    
    ref_path = arcpy.Describe(ref).path
    ref_name = arcpy.Describe(ref).name
    scnro_path = arcpy.Describe(scnro).path
    scnro_name = arcpy.Describe(scnro).path
#Get empty lists to populate
    refHabStateAreaValues = [0,0,0,0,0]
    refRiskAreaValues = [0,0,0,0,0,0,0,0,0,0,0]
    refRsfAreaValues = [0,0,0,0,0,0,0,0,0,0]
#
    scnroHabStateAreaValues = [0,0,0,0,0]
    scnroRiskAreaValues = [0,0,0,0,0,0,0,0,0,0,0]
    scnroRsfAreaValues = [0,0,0,0,0,0,0,0,0,0]
#    
    HabStatesValuesCheck = [-2,-1,0,1,2]
    RiskValuesCheck = [0,1,2,3,4,5,6,7,8,9,10]
    RsfValuesCheck = [1,2,3,4,5,6,7,8,9,10]
#      
#Get cell size
    CellSize = arcpy.GetRasterProperties_management(ref, "CELLSIZEX")
    CellSize_x2 = float(CellSize.getOutput(0)) * float(CellSize.getOutput(0))
#Get reference VALUE area from COUNT
    if name == "Habitat States":
        with arcpy.da.SearchCursor(ref, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                index = HabStatesValuesCheck.index(row[0])
                ref_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                refHabStateAreaValues[index] = ref_sum_area
        del row
        del cursor
#
    elif name == "Risk":
        with arcpy.da.SearchCursor(ref, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                index = RiskValuesCheck.index(row[0])
                ref_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                refRiskAreaValues[index] = ref_sum_area
        del row
        del cursor
#
    else:
        with arcpy.da.SearchCursor(ref, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                if row[0] == 0:
                    ref_sum_area_0 = round(((row[1] * CellSize_x2) * 0.000001),2)
#                    
                elif row[0] == 1:
                    try:
                        index = RsfValuesCheck.index(row[0])
                        ref_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                        ref_sum_area += ref_sum_area_0
                        refRsfAreaValues[index] = ref_sum_area
                    except:
                        index = RsfValuesCheck.index(row[0])
                        ref_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                        refRsfAreaValues[index] = ref_sum_area
#               
                else:
                    index = RsfValuesCheck.index(row[0])
                    ref_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                    refRsfAreaValues[index] = ref_sum_area
        del row
        del cursor   
#----------------------------------------------------------------------------------------
#Get scenario VALUE area from COUNT
    if name == "Habitat States":
        with arcpy.da.SearchCursor(scnro, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                index = HabStatesValuesCheck.index(row[0])
                scnro_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                scnroHabStateAreaValues[index] = scnro_sum_area
        del row
        del cursor
#
    elif name == "Risk":
        with arcpy.da.SearchCursor(scnro, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                index = RiskValuesCheck.index(row[0])
                scnro_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                scnroRiskAreaValues[index] = scnro_sum_area
        del row
        del cursor
#
    else:
        with arcpy.da.SearchCursor(scnro, ["VALUE","COUNT"]) as cursor:
            for row in cursor:
                if row[0] == 0:
                    scnro_sum_area_0 = round(((row[1] * CellSize_x2) * 0.000001),2)
#                    
                elif row[0] == 1:
                    try:
                        index = RsfValuesCheck.index(row[0])
                        scnro_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                        scnro_sum_area += scnro_sum_area_0
                        scnroRsfAreaValues[index] = scnro_sum_area
                    except:
                        index = RsfValuesCheck.index(row[0])
                        scnro_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                        scnroRsfAreaValues[index] = scnro_sum_area
#               
                else:
                    index = RsfValuesCheck.index(row[0])
                    scnro_sum_area = round(((row[1] * CellSize_x2) * 0.000001),2)
                    scnroRsfAreaValues[index] = scnro_sum_area
        del row
        del cursor
#----------------------------------------------------------------------------------------
#Get final ref / scenario variables set
    if name == "Habitat States":
        refAreaValues = refHabStateAreaValues
        scnroAreaValues = scnroHabStateAreaValues
        xValueNames = ["Primary\nSink\n-2", "Secondary\nSink\n-1",
                       "Non-Critical\nHabitat\n0",
                       "Secondary\nHabitat\n1", "Primary\nHabitat\n2"]
        txtValueNames = ["-2", "-1", "0","1","2"]
#
    elif name == "Risk":
        refAreaValues = refRiskAreaValues
        scnroAreaValues = scnroRiskAreaValues
        xValueNames = ["Low\n0", "Low\n1", "Low\n2", "Low\n3","Med\n4","Med\n5","Med\n6",
                       "High\n7", "High\n8", "High\n9", "High\n10"]
        txtValueNames = ["L-0","L-1","L-2","L-3","M-4","M-5","M-6",
                         "H-7","H-8","H-9","H-10"]
                
#
    else:
        refAreaValues = refRsfAreaValues
        scnroAreaValues = scnroRsfAreaValues
        xValueNames = ["Low\n1", "Low\n2", "Low\n3", "Med\n4", "Med\n5", "Med\n6",
                       "High\n7", "High\n8", "High\n9", "High\n10"]
        txtValueNames = ["L-1","L-2","L-3","M-4","M-5","M-6","H-7","H-8","H-9","H-10"]  
#----------------------------------------------------------------------------------------
#Calclate percent change between reference and scenario   
    pctChange = []
    for refArea,scnroArea in zip(refAreaValues, scnroAreaValues):
        #Can't divide by 0
        if refArea == 0 and scnroArea == 0:
            pChange = 0
            pctChange.append(pChange)
        #Can't divide by 0
        elif refArea == 0 and scnroArea > 0:
            pChange = 100
            pctChange.append(pChange)
        else:
            pChange = round((((scnroArea - refArea) / refArea) * 100),1)
            pctChange.append(pChange)
#----------------------------------------------------------------------------------------
#Get reference / scneario mean value
    refMean = arcpy.GetRasterProperties_management(ref,"MEAN")
    refMean = round(float(refMean.getOutput(0)),2)
    scnroMean = arcpy.GetRasterProperties_management(scnro,"MEAN")
    scnroMean = round(float(scnroMean.getOutput(0)),2)
#----------------------------------------------------------------------------------------
#Delete scratch rasters, values have been extracted
    arcpy.Delete_management(ref)
    arcpy.Delete_management(scnro)    
#----------------------------------------------------------------------------------------
#Get total area of  scenario/reference raster
    totalArea = sum(scnroAreaValues)
#----------------------------------------------------------------------------------------
#Get total percent area change
    if sum(ChangeDetectionArea) == 0:
        totalPctChange = 0
        ChangeDetectionArea = 0
    else:    
        totalPctChange = round((sum((ChangeDetectionArea))/(totalArea)* 100),1)
        ChangeDetectionArea = round(sum(ChangeDetectionArea),1)
        totalArea = round(sum(scnroAreaValues),1)
#----------------------------------------------------------------------------------------
#Generate Mean value / total sum area string
    meanText = ("Total Area: " + str(totalArea) + " km2" + 
                "\nChanged Area: " + str(ChangeDetectionArea) + " km2" +
                "\nTotal %Changed: " + str(totalPctChange) + "%" +
                "\n\nMean {} Value: ".format(refUserName[:10]) + str(refMean) +
                "\nMean {} Value: ".format(scnroUserName[:10]) + str(scnroMean))   
#----------------------------------------------------------------------------------------
#Generate percent change text string
    valueStr = "Value"
    pctStr = " Change"
    for xValueName, pct in zip(txtValueNames,pctChange):
        xValueName = str(xValueName)
        valueStr += '\n' + xValueName + ":" 
        pct = str(pct)
        pctStr += '\n' +  pct + "%"   
#----------------------------------------------------------------------------------------
#Plot
    N = len(xValueNames)
    ind = np.arange(N)
    width = 0.4
    plt.bar(ind, refAreaValues, width, label=refUserName[:10],
            color = '#fc8d62', edgecolor="black", linewidth=2)
    plt.bar(ind + width, scnroAreaValues, width, label=scnroUserName[:10],
            color = '#66c2a5', edgecolor="black", linewidth=2)
    plt.margins(0.02)
    plt.xticks(ind + (width / 2) + 0.2, xValueNames)
    plt.ylabel('Area, Km2')
    #plt.title("Percent Area Changed")
    plt.text(.99, 0.35, valueStr, transform=plt.gcf().transFigure,
             horizontalalignment='right')
    plt.text(1.1, 0.35, pctStr, transform=plt.gcf().transFigure,
             horizontalalignment='right')
    plt.text(0.95, 0.1, meanText, transform=plt.gcf().transFigure)
    plt.legend(loc='best')
    plt.grid(color='black', which='major', axis='y', linestyle='dashed')
    plt.savefig(fp(plotFolder, name + ".png"), bbox_inches = "tight")
    pngPaths.append(fp(plotFolder, name + ".png"))
    plt.close('all')#Close figure
#----------------------------------------------------------------------------------------
    #Add plot to .mxd template layout
    for element in arcpy.mapping.ListLayoutElements(mxd, "PICTURE_ELEMENT"):
        if element.name == "Change":
            element.sourceImage = fp(plotFolder, name + ".png")
    #Change title
    for element in arcpy.mapping.ListLayoutElements(mxd, "TEXT_ELEMENT"):
        if element.name == "Title":
            element.text = "Change Detection: " + name
    arcpy.mapping.ExportToPNG(mxd, fp(plotFolder, name + ".png"),
                              "PAGE_LAYOUT", resolution=300)
#----------------------------------------------------------------------------------------
    #Display plots
    try:
        path = fp(plotFolder, name + ".png")
        imageViewerFromCommandLine = {'win32':'explorer',
                                  'win64':'explorer',
                                  'linux': 'eog'}[sys.platform]
    
        subprocess.call([imageViewerFromCommandLine, path])
    except:
        continue

del updateLayer        
del mxd
#
#----------------------------------------------------------------------------------------
#Step 7: Clean scratch folder
arcpy.AddMessage(" ")
arcpy.AddMessage("Step 7:   Cleaning scratch folder")
#========================================================================================
#Clean scratch folder if it exists
gc.collect() 
arcpy.ClearWorkspaceCache_management() 
if arcpy.Exists(fp(scnroHabStates_folder, "Scratch")):
    path = fp(scnroHabStates_folder, "Scratch")
    shutil.rmtree(path, ignore_errors=True)
#    
#----------------------------------------------------------------------------------------
#========================================================================================
            
