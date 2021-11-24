#----------------------------------------------------------------------------------------
#                                       fRI RESEARCH
#                                   GRIZZLY BEAR PROGRAM
#----------------------------------------------------------------------------------------
# SCRIPT NAME:  CheckLicenses.py
#               v.2018-1018
#
# PURPOSE:      Checks for license availability.
#               Returns "yes" or "no" value and displays message.
# 
# AUTHOR:       Julie Duval
#               Senior GIS Analyst
#               Foothills Research Institute
#
# EDITORS:      Josh Crough, GIS Analyst
#               Foothills Research Institute
#
# CREATED:      September 2011
#
# LAST UPDATES: January 29, 2015 (J. Crough) - Converted scripts to arcpy
#               October 18, 2018 (J. Duval)  - fixed references to gp. (to .arcpy)
#
# ---------------------------------------------------------------------------------------
# SCRIPT FUNCTIONS:
#       main(args)                 - Main function 
#       CheckArcInfo()             - Checks for ArcInfo availability
#       CheckArcEditor()           - Checks for ArcEditor availability
#       CheckArcView()             - Checks for ArcView availability
#       CheckSpatialExt()          - Checks for Spatial Analyst availability
#       Check3DExt()               - Checks for 3D Analyst availability
#       display(msg)        - Displays message as tool output message, or python print
#========================================================================================
import arcpy

class LicenseError(Exception):
    pass

class UnlicensedError(Exception):
    pass

#****************************************************************************************
#************************************ LICENSE CHECKS ************************************
#****************************************************************************************
def CheckArcInfo(go="yes"):
    #------------------------------------------------------------------------------------
    'Function to checkout an ArcInfo level license'
    #------------------------------------------------------------------------------------
    try:
        if arcpy.CheckProduct("ArcInfo") == "Available" or \
           arcpy.CheckProduct("ArcInfo") == "AlreadyInitialized":
            arcpy.SetProduct("ArcInfo")
            display("*** ArcInfo License is available")
        elif arcpy.CheckProduct("ArcInfo") == "NotLicensed": raise UnlicensedError
        else: raise LicenseError

    except UnlicensedError:    
        display("*** Required ArcInfo license was not found.")
        go = "no"
    except LicenseError:    
        display("*** ArcInfo license is unavailable and is required for this tool to run")
        go = "no"
    except:
        print arcpy.GetMessages(2); raise

    return (go)    


def CheckArcEditor(go="yes"):
    #------------------------------------------------------------------------------------
    'Function to checkout out an ArcEditor level license'
    #------------------------------------------------------------------------------------
    try:
        if arcpy.CheckProduct("ArcEditor") == "Available" or \
           arcpy.CheckProduct("ArcEditor") == "AlreadyInitialized":
            arcpy.SetProduct("ArcEditor")
            display("*** ArcEditor License is available")
        elif gp.CheckProduct("ArcEditor") == "NotLicensed": raise UnlicensedError
        else: raise LicenseError

    except UnlicensedError:    
        display("*** Required ArcEditor license was not found.")
        go = "no"
    except LicenseError:    
        display("*** ArcEditor license is unavailable and is required for this "
                "tool to run")
        go = "no"
    except:
        print arcpy.GetMessages(2); raise

    return (go)

def CheckArcView(go="yes"):
    #------------------------------------------------------------------------------------
    'Fucntion to checkout an ArcView level license'
    #------------------------------------------------------------------------------------
    try:
        if arcpy.CheckProduct("ArcView") == "Available" or \
           arcpy.CheckProduct("ArcView") == "AlreadyInitialized":
            arcpy.SetProduct("ArcView")
            display("*** ArcView License is available")
        elif arcpy.CheckProduct("ArcView") == "NotLicensed": raise UnlicensedError
        else: raise LicenseError

    except UnlicensedError:    
        display("*** Required ArcView license was not found.")
        go = "no"
    except LicenseError:    
        display("*** ArcView license is unavailable and is required for this tool to run")
        go = "no"
    except:
        print arcpy.GetMessages(2); raise

    return (go)


def CheckSpatialExt(go="yes"):
    #------------------------------------------------------------------------------------
    'Function to checkout a spatial analysit extension license'
    #------------------------------------------------------------------------------------
    try:
        if arcpy.CheckExtension("spatial") == "Available":
            arcpy.CheckOutExtension("spatial")
            display("*** Spatial Analyst license is available")
        elif arcpy.CheckOutExtension("spatial") == "NotLicensed": raise UnlicensedError
        else: raise LicenseError

    except UnlicensedError:    
        display("*** Required Spatial Analyst Extension license was not found.")
        go = "no"
    except LicenseError:    
        display("*** Spatial Analyst license is unavailable and is required for this "
                "tool to run")
        go = "no"
    except:
        print arcpy.GetMessages(2); raise

    return (go)


def Check3DExt(go="yes"):
    #------------------------------------------------------------------------------------
    'Function to checkout a 3D analysit extension license'
    #------------------------------------------------------------------------------------
    try:
        if arcpy.CheckExtension("3D") == "Available":
            arcpy.CheckOutExtension("3D")
            display("*** 3D Analyst license is available")
        elif arcpy.CheckOutExtension("3D") == "NotLicensed": raise UnlicensedError
        else: raise LicenseError

    except UnlicensedError:    
        display("*** Required 3D Analyst Extension license was not found.")
        go = "no"
    except LicenseError:    
        display("*** Spatial 3D license is unavailable and is required for this "
                "tool to run")
        go = "no"
    except:
        print arcpy.GetMessages(2); raise

    return (go)


def display(msg):
    arcpy.AddMessage(msg)
    print msg


def main():
    print "This module checks the availability of ArcGIS licenses and extensions."


if __name__ == '__main__':
    main()  
