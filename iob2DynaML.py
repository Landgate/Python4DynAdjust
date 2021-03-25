
# ----------------------------------------------------------------------
#                          iob2DynaML.py
# ----------------------------------------------------------------------
#  Author:  Kent Wheeler
#    Date:  17 March 2020
# Purpose:  Script to create DynaML stn and msr files from Geolab .iob file
# ----------------------------------------------------------------------
#   Usage:  cmd:\> python iob2DynaML.py <*.iob>
# ----------------------------------------------------------------------
    
import os, sys, shutil
from src.Python4Dyna import * 
from src.Dyna2db import * 
from src.geolab import * 
    
if __name__ == "__main__":
    script_path = os.path.abspath(os.path.realpath(__file__))
    script_dir, script_name = os.path.split(script_path)
    os.chdir(script_dir)
    
    files = os.listdir(os.getcwd())
    if len(sys.argv)>1: files=[sys.argv[1]]
    for df in files:
        if df.endswith(".iob"):
            network=df.replace('.iob','').replace(' ','_')
            print('  Converting to DynaML: ' + network)
            datum, geoid = iob2DynaML(df)
            print('  Adjusting: ' + network)
            LSA(network, datum, geoid)
            
            network=df.replace('.iob','.phased-stage').replace(' ','_')
            if os.path.exists(network + '.adj'):
                db='network.db'
                create_DynaML_db(db)
                adj_stats, adj_stns, adj_obs=import_adj (network + '.adj',db)
                if os.path.exists(network + '.apu'):
                    apu_stns=import_apu(network + '.apu',db)
                    export_lst(adj_stats,adj_stns,apu_stns,adj_obs)
                else:
                    print('Lst File not printed...apu not created')
                coords, xyz_bases = adj2Dicts(network + '.adj')
                if len(xyz_bases) > 0: 
                    coords, bases = gnss_xyz2dah(coords, xyz_bases)
                    Create_DynaBAS(network, bases)
                    #DynaML2db(network)
                    #db2Trivials(network+'.db')
            else:
                print('Error in Adjustment...Check Dynadjust reported errors')
                
    files = os.listdir(os.getcwd())
    for f in files:
        if f.endswith('.aml'):os.remove(f)
        if f.endswith('.asl'):os.remove(f)
        if f.endswith('.bms'):os.remove(f)
        if f.endswith('.bst'):os.remove(f)
        if f.endswith('.dbid'):os.remove(f)
        if f.endswith('.dnaproj'):os.remove(f)
        if f.endswith('.imp'):os.remove(f)
        if f.endswith('.map'):os.remove(f)
        if f.endswith('.mtx'):os.remove(f)
        if f.endswith('.seg'):os.remove(f)
        if f.endswith('.xyz'):os.remove(f)
        if f.endswith('.db'):os.remove(f)
        if f.endswith('.lg.xml'):os.remove(f)
        if f.endswith('.adj') or f.endswith('.apu'):
            shutil .move(f, f.replace('.phased-stage',''))
    
    print('\n============== Complete ==============')
