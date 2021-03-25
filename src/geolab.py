from src.Python4Dyna import * 
import os

def iob2DynaML(iob_file):
    #############################################################################
    ##     Input: iobfile.iob                                                  ##
    ##     Output: iobfile.msr.xml (Created in same directory)                 ##
    ##             iobfile.stn.xml (Created in same directory)                 ##
    ##     Return: datum, geoid (These are extraced from the iob and are based ##
    ##                            on the geoid used in the iob)                ##
    #############################################################################
    rel_uncert_sum = 0
    avg_rel_uncert = 0
    stn_cnt = 0
    # Open, run through and close the file
    # for initial information on the adjustment
    print('Reading the Geolab (.iob) File...')
    with open(iob_file, 'r') as f:
        GNSSmarksStr = ';'
        FloatmarksStr = ';'
        w_MksWithObs = ';'
        nme_rec={}
        for ln in f.readlines():
            if ln[0:4] == ' PLO' or ln[0:4] == ' PLH':
                if (ln[72:len(ln)].strip() != ''
                        and ln.find('ppm') != -1):
                    stnRec = ln[72:len(ln)].strip().replace(';', '|').split('|')
                    rel_uncert_sum = (rel_uncert_sum +
                                      float(stnRec[7].replace('ppm', '')))
                    stn_cnt = stn_cnt + 1
            if ln[0:6]==' GFIL ': 
                geoid=ln[6:].strip().replace('P:','\\\\dli\public\\Business Unit\\Operations\\Registrations\\CadastralSubdivisions\\Geodetic\\sipublic1')
                geoid_path='\\'.join(geoid.split('\\')[0:-1])
                gfiles = os.listdir(geoid_path)
                datum='GDA2020'
                if geoid_path.find('09')!=-1:datum='GDA94'
                for g in gfiles:
                    if g.endswith('.gsb'):geoid=geoid_path +'\\' + g
            if (ln[0:5] == ' OHDF' or
                ln[0:5] == ' GAZI' or
                ln[0:5] == ' AZIM' or
                ln[0:5] == ' DIST' or
                ln[0:4] == ' DIR' or
                    ln[0:5] == ' DXYZ'):
                CurrentMsr = DnaMeasurement()
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                if ln[0:5] == ' DXYZ':
                    GNSSmarksStr = (GNSSmarksStr + ';' +
                                    CurrentMsr.first + ';' + CurrentMsr.second)
                if (FloatmarksStr.find(';' + CurrentMsr.first+';') != -1 or
                        FloatmarksStr.find(';'+CurrentMsr.second + ';') != -1):
                    w_MksWithObs = (w_MksWithObs + CurrentMsr.first + ';' +
                                    CurrentMsr.second + ';')

    # bin the guessing of GDA94 relative uncertainty for the float marks
    if stn_cnt != 0:
        avg_rel_uncert = rel_uncert_sum / stn_cnt
        if avg_rel_uncert < 3:
            avg_rel_uncert = 3
        if 3 < avg_rel_uncert <= 7:
            avg_rel_uncert = 7.5
        if 7.5 < avg_rel_uncert <= 10:
            avg_rel_uncert = 10
        if 10 < avg_rel_uncert <= 20:
            avg_rel_uncert = 20
        if 20 < avg_rel_uncert <= 30:
            avg_rel_uncert = 30
        if 30 < avg_rel_uncert <= 50:
            avg_rel_uncert = 50
        if avg_rel_uncert > 50:
            avg_rel_uncert = int(avg_rel_uncert / 10) * 10
    with open(iob_file, 'r') as f:
        stnout = open(iob_file.replace('.iob', '.stn.xml')
                                .replace(' ', '_'), 'w')
        msrout = open(iob_file.replace('.iob', '.msr.xml')
                                .replace(' ', '_'), 'w')
        stnout.write(stn_header(datum))
        msrout.write(msr_header())

        # Run through each line of the input file and extract the relevant lines
        lineCount = 0; msr_cnt=0
        InstHts = DeviceHeight()
        TgtHts = DeviceHeight()
        ControlRec = AdditionalInfoMsr()
        for ln in f.readlines():
            if ln[0:5] == ' TITL':
                jobNumber = find_job_num(ln)
                if jobNumber == '':
                    jobNumber = find_job_num(os.getcwd())
            if ln[0:4] == ' PLO' or ln[0:4] == ' PLH':
                CurrentStn = DnaStation()
                CurrentStn.Name = ln[10:23].strip()
                CurrentStn.Const = ln[6:9].replace('1', 'C')
                CurrentStn.Const = CurrentStn.Const.replace('0', 'F')
                CurrentStn.Const = CurrentStn.Const.strip()
                if ln[0:4] == ' PLO':
                    CurrentStn.Type = 'LLH'
                if ln[0:4] == ' PLH':
                    CurrentStn.Type = 'LLh'
                CurrentStn.XAxis = hms2hp(ln[23:41].strip())
                CurrentStn.YAxis = hms2hp(ln[41:59].strip())
                CurrentStn.Height = ln[59:72].strip()
                CurrentStn.Desc = ln[72:len(ln)].strip()
                stn_addn_rec = AdditionalInfoStn()
                stn_addn_rec.NonGSNumber = 'E' + jobNumber
                if CurrentStn.Desc != '':
                    stnRec = CurrentStn.Desc.replace(';', '|').split('|')
                    stn_addn_rec.Ges_Name = stnRec[0]
                    stn_addn_rec.Ges_Lat = stnRec[1]
                    stn_addn_rec.Ges_Lng = stnRec[2]
                    stn_addn_rec.Ges_Ht = stnRec[3]
                    stn_addn_rec.Ges_Ht_Acc = stnRec[4]
                    stn_addn_rec.Ges_Ht_Mth = stnRec[5]
                    stn_addn_rec.Ges_Hz_Ord = stnRec[6]
                    stn_addn_rec.rel_hz_acc = stnRec[7]
                    stn_addn_rec.Hz_Coord_Mth = stnRec[8]
                    if stnRec[10] != '':
                        CurrentStn.CE = stnRec[9]
                        CurrentStn.aAxis = float(stnRec[10])
                        CurrentStn.bAxis = float(stnRec[11])
                        CurrentStn.ErrAz = float(stnRec[12])
                    stn_addn_rec.SelectPoint = 'true'
                    stn_addn_rec.SelectRL = 'true'

                if CurrentStn.Const[0:2] == 'FF' and avg_rel_uncert != 0:
                    stn_addn_rec.rel_hz_acc = str(avg_rel_uncert) + 'ppm'
                    if GNSSmarksStr.find(';'+CurrentStn.Name+';'):
                        stn_addn_rec.Hz_Coord_Mth = 'GNSS'

                stnout.write(stn_xml_str(CurrentStn, stn_addn_rec))

            if ln[0:5] == ' HI  ':
                add_device_height(InstHts, ln[10:23].strip(), ln[23:33].strip())

            if ln[0:5] == ' HT  ':
                add_device_height(TgtHts, ln[10:23].strip(), ln[23:33].strip())

            if (ln[0:5] == ' OHDF' or
                ln[0:5] == ' OHGT' or
                ln[0:5] == ' GAZI' or
                ln[0:5] == ' AZIM' or
                ln[0:5] == ' DIST' or
                ln[0:5] == ' DSET' or
                ln[0:5] == ' GRP '):
                    msr_cnt+=1
                    CurrentMsr = DnaMeasurement()
                    CurrentMsr.id = msr_cnt
                    
            if ln[0:5] == ' OHGT':
                CurrentMsr.type = 'H'
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.add_values(ln[36:65].strip())
                CurrentMsr.stddev = ln[65:76].strip()
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))

            if ln[0:5] == ' OHDF':
                CurrentMsr.type = 'L'
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.add_values(ln[50:65].strip())
                CurrentMsr.stddev = ln[65:76].strip()
                ControlRec.LevelDistance = ln[36:50].strip()
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))

            if ln[0:5] == ' GAZI':
                CurrentMsr.type = 'B'
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                CurrentMsr.stddev = ln[65:76].strip()
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))

            if ln[0:5] == ' AZIM':
                CurrentMsr.type = 'K'
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                CurrentMsr.stddev = ln[65:76].strip()
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))

            if ln[0:5] == ' DIST':
                CurrentMsr.type = 'S'
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.add_values(ln[36:65].strip())
                rw = 0
                for Stn in InstHts.StnName:
                    if Stn == CurrentMsr.first:
                        CurrentMsr.instheight = InstHts.RefHeight[rw]
                    rw = rw+1
                rw = 0
                for Stn in TgtHts.StnName:
                    if Stn == CurrentMsr.first:
                        CurrentMsr.targheight = TgtHts.RefHeight[rw]
                    rw = rw+1
                CurrentMsr.stddev = ln[65:76].strip()
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                
            if ln[0:5] == ' DSET':
                CurrentMsr.type = 'D'
                lineCount = 0
            if ln[0:4] == ' DIR' and lineCount == 1 and CurrentMsr.type == 'D':
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.add_targets(ln[23:35].strip())
                CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                CurrentMsr.add_targetstddevs(ln[65:76].strip())
                CurrentMsr.total = lineCount
            if ln[0:4] == ' DIR' and lineCount > 1 and CurrentMsr.type == 'D':
                CurrentMsr.add_targets(ln[23:35].strip())
                CurrentMsr.add_values(hms2hp(ln[36:65].strip()))
                CurrentMsr.add_targetstddevs(ln[65:76].strip())
                CurrentMsr.total = lineCount
            if CurrentMsr.type == 'D' and ln[0:4] != ' DIR' and lineCount > 1:
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                CurrentMsr = DnaMeasurement()

        # Scrape information from Landgate GESMAR control Records
        # eg.*CONTROL;GPS;201810010007;012057;E;B;TRIM;TRIM;D;ST  ;FX;015;N
        # eg.*CONTROL;OHDF;20181001;SLEV;LC;2.71;MEAS
        # eg.*CONTROL;DIS;20181213;TS 16;C;MEAS
        # eg.*CONTROL;ANG;20181213;TS16;C;MEAS
            if ln.startswith('*CONTROL;') and ln.startswith('*CONTROL;NAM')==False:
                ControlRec = AdditionalInfoMsr()
                alinestr = ln.split(';')
                stdatetimestr = alinestr[2] + '0000'
                yr = int(stdatetimestr[0:4])
                mth = int(stdatetimestr[4:6])
                ddy = int(stdatetimestr[6:8])
                hr = int(stdatetimestr[8:10])
                mn = int(stdatetimestr[10:12])
                ControlRec.StartDateTime = dt.datetime(yr, mth, ddy, hr, mn, 00)
                ControlRec.NonGSNumber = 'E' + jobNumber
                if ln[0:13] == '*CONTROL;DIS;':
                    ControlRec.InstrumentModel = alinestr[3].strip()
                    ControlRec.Class = alinestr[4].strip()
                    ControlRec.Derivation = alinestr[5].strip()
                if ln[0:13] == '*CONTROL;ANG;':
                    ControlRec.InstrumentModel = alinestr[3].strip()
                    ControlRec.Class = alinestr[4].strip()
                    ControlRec.Derivation = alinestr[5].strip()
                if ln[0:14] == '*CONTROL;OHDF;':
                    ControlRec.SurveyTechnique = alinestr[3].strip()
                    ControlRec.Class = alinestr[4].strip()
                    ControlRec.LevelDistance = alinestr[5].strip()
                    ControlRec.Derivation = alinestr[6].strip()
                if ln[0:13] == '*CONTROL;GPS;':
                    durationstr = alinestr[3]
                    hr = int(durationstr[0:2])
                    mn = int(durationstr[2:4])
                    sec = int(durationstr[4:6])
                    ControlRec.Duration = (dt.datetime(1900, 1, 1, 0, 0, 0) +
                                           dt.timedelta(hours=hr,
                                                        minutes=mn,
                                                        seconds=sec))
                    ControlRec.TimeStatus = alinestr[4].strip()
                    ControlRec.EphemerisType = alinestr[5].strip()
                    ControlRec.AtReceiver = alinestr[6].strip()
                    ControlRec.ToReceiver = alinestr[7].strip()
                    ControlRec.FrequencyMode = alinestr[8].strip()
                    ControlRec.SurveyTechnique = alinestr[9].strip()
                    ControlRec.Solution = alinestr[10].strip()
                    if alinestr[11].strip().isnumeric():
                        ControlRec.EpochInterval = int(alinestr[11])
                    ControlRec.Class = alinestr[12].strip()

            if ln[0:4] == ' GRP':
                CurrentMsr.type = 'G'
            if ln[0:5] == ' DXYZ':
                CurrentMsr.first = ln[10:23].strip()
                CurrentMsr.second = ln[23:35].strip()
                CurrentMsr.dx = ln[36:50].strip()
                CurrentMsr.dy = ln[50:64].strip()
                CurrentMsr.dz = ln[64:78].strip()
            if ln[0:6] == ' COV  ' or ln[0:6] == ' CORR ':
                m_type = ln[0:13]
                lineCount = 0
                if ln[25:36].strip()!='': CurrentMsr.vscale = ln[25:36].strip()
                diagscale=1
                if ln[48:58].strip()!='': diagscale = float(ln[48:58].strip())
                if diagscale == 0: diagscale=1
            if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 1:
                CurrentMsr.Vs[0, 0] = float(ln[7:30].strip())*diagscale
                if m_type== ' COV  LG DIAG':
                    CurrentMsr.Vs[1, 1] = ln[30:54].strip()
                    CurrentMsr.Vs[2, 2] = ln[54:78].strip()
                    msrout.write(msr_xml_str(CurrentMsr, ControlRec))
                else:    
                    CurrentMsr.Vs[0, 1] = ln[30:54].strip()
                    CurrentMsr.Vs[0, 2] = ln[54:78].strip()
            if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 2:
                CurrentMsr.Vs[1, 1] = float(ln[7:30].strip())*diagscale
                CurrentMsr.Vs[1, 2] = ln[30:54].strip()
            if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 3:
                CurrentMsr.Vs[2, 2] = float(ln[7:30].strip())*diagscale
                if m_type== ' COV  CT UPPR':
                    msrout.write(msr_xml_str(CurrentMsr, ControlRec))
            if ln[0:5] == ' ELEM' and CurrentMsr.type == 'G' and lineCount == 4:
                Ds=np.zeros([3,3])
                Ds[0,0] = ln[7:30].strip()
                Ds[1,1] = ln[30:54].strip()
                Ds[2,2] = ln[54:78].strip()
                CurrentMsr.Vs[1, 0] =CurrentMsr.Vs[0, 1]
                CurrentMsr.Vs[2, 0] =CurrentMsr.Vs[0, 2]
                CurrentMsr.Vs[2, 1] =CurrentMsr.Vs[1, 2]
                CurrentMsr.Vs=np.matmul(np.matmul(Ds,CurrentMsr.Vs),Ds)
                msrout.write(msr_xml_str(CurrentMsr, ControlRec))

            lineCount = lineCount+1

        # Write footers, close
        stnout.write(dml_footer())
        msrout.write(dml_footer())
        stnout.close()
        msrout.close()
        
        return datum, geoid

def export_lst(adj_stats,adj_stns,apu_stns,adj_obs):
    lst_file=adj_stats['FILE_NAME'][:-3]+'lst'
    with open(lst_file, 'w') as f:
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write(adj_stats['FILE_CREATED']+'\n')
        f.write('Input file:   '+adj_stats['FILE_NAME']+'\n')
        f.write('Output file:  '+lst_file+'\n')
        f.write('Options file: '+adj_stats['COMMAND_LINE_ARGUMENTS']+'\n')
        f.write('Geoid File:   '+adj_stats['GEOID_MODEL']+'\n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('|             PARAMETERS               |             OBSERVATIONS             |\n')
        f.write('|-----------------------------------------------------------------------------|\n')
        f.write('|   Description    |      Number       |   Description     |      Number      |\n')
        f.write('|-----------------------------------------------------------------------------|\n')
        f.write('| No. of Stations  |{:12d}       | Directions        |{:12d}      |\n'.format(len(adj_stns)
                                                                                        ,sum(1 for i in adj_obs.values() if i['M']=='D')))
        f.write('| Coord Parameters |{:12d}       | Distances         |{:12d}      |\n'.format(len(adj_stns)*2,
                                                                                        sum(1 for i in adj_obs.values() if (i['M']=='M' or i['M']=='S'))))
        f.write('| Free Latitudes   |{:12d}       | Azimuths          |{:12d}      |\n'.format(sum(1 for i in adj_stns.values() if i['CONST'][:1]=='F'),
                                                                                        sum(1 for i in adj_obs.values() if (i['M']=='B' or i['M']=='K'))))
        f.write('| Free Longitudes  |{:12d}       | Vertical Angles   |           0      |\n'.format(sum(1 for i in adj_stns.values() if i['CONST'][1:2]=='F')))
        f.write('| Free Heights     |{:12d}       | Zenithal Angles   |           0      |\n'.format(sum(1 for i in adj_stns.values() if i['CONST'][-1:]=='F')))
        f.write('| Fixed Coordinates|{:12d}       | Angles            |           0      |\n'.format(sum(1 for i in adj_stns.values() if i['CONST'][:2]=='CC')))
        f.write('| Astro. Latitudes |           0       | Heights           |{:12d}      |\n'.format(sum(1 for i in adj_obs.values() if i['M']=='H')))
        f.write('| Astro. Longitudes|           0       | Height Differences|           0      |\n')
        f.write('| Geoid Records    |           0       | Auxiliary Params. |           0      |\n')
        f.write('| All Aux. Pars.   |           0       | 2-D Coords.       |           0      |\n')
        f.write('| Direction Pars.  |           0       | 2-D Coord. Diffs. |           0      |\n')
        f.write('| Scale Parameters |           0       | 3-D Coords.       |           0      |\n')
        f.write('| Constant Pars.   |           0       | 3-D Coord. Diffs. |{:12d}      |\n'.format(sum(1 for i in adj_obs.values() if i['M']=='G' or i['M']=='X')))
        f.write('| Rotation Pars.   |           0       |                   |                  |\n')
        f.write('| Translation Pars.|           0       |                   |                  |\n')
        f.write('|                  |                   |                   |                  |\n')
        f.write('|                  |    --------       |                   |    --------      |\n')
        f.write('| Total Parameters |{:>12}       | Total Observations|{:12d}      |\n'.format(adj_stats['NUMBER_OF_UNKNOWN_PARAMETERS'],len(adj_obs)))
        f.write('|-----------------------------------------------------------------------------|\n')
        f.write('|                      Degrees of Freedom ={:>10}                         |\n'.format(adj_stats['DEGREES_OF_FREEDOM']))
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('-------------------------------------------------------------------------------\n')
        f.write('                        SUMMARY OF SELECTED OPTIONS                            \n')
        f.write('-------------------------------------------------------------------------------\n')
        f.write('    OPTION                             |   SELECTION                           \n')
        f.write('-------------------------------------------------------------------------------\n')
        f.write('    Computation Mode                   |   Adjustment                          \n')
        f.write('    Maximum Iterations                 |   {:36}\n'.format(adj_stats['MAXIMUM_ITERATIONS']))
        f.write('    Convergence Criterion              |   {:36}\n'.format(adj_stats['ITERATION_THRESHOLD']))
        f.write('    Residual Rejection Criterion       |   Tau Max                             \n')
        f.write('    Confidence Region Types            |   1D 2D Station Relative              \n')
        f.write('    Relative Confidence Regions        |   Connected Only                      \n')
        f.write('    Variance Factor (VF) Known         |   Yes                                 \n')
        f.write('    Scale Covariance Matrix With VF    |   Yes                                 \n')
        f.write('    Scale Residual Variances With VF   |   No                                  \n')
        f.write('    Force Convergence in Max Iters     |   No                                  \n')
        f.write('    Distances Contribute To Heights    |   No                                  \n')
        f.write('    Compute Full Inverse               |   Yes                                 \n')
        f.write('    Optimize Band Width                |   Yes                                 \n')
        f.write('    Generate Initial Coordinates       |   No                                  \n')
        f.write('    Re-Transform Obs After 1st Pass    |   Yes                                 \n')
        f.write('    Geoid Interpolation Method         |   Bi-Cubic                            \n')
        f.write('-------------------------------------------------------------------------------\n')
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write('Adjusted PLH Coordinates:                                                                            \n')
        f.write('                                                                         LATITUDE LONGITUDE   HEIGHT \n')
        f.write(' CODE FFF STATION               LATITUDE         LONGITUDE  ELIP-HEIGHT   STD DEV   STD DEV   STD DEV\n')    
        f.write(' ---- --- ------------ ----------------- ----------------- ------------ --------- --------- ---------\n')
        for s in adj_stns.values():
            f.write(' PLH  {:3} {:12} S {:>15} E{:>16} {:>12} {:>9} {:>9} {:>9}\n'.format(
                    s['CONST'].replace('C','1').replace('F','0'),
                    s['STATION'],
                    dd2hms(-1*float(s['LATITUDE']),5),
                    dd2hms(float(s['LONGITUDE']),5),
                    s['H_ELLIPSE'],
                    s['SD_E'],
                    s['SD_N'],
                    s['SD_UP']))
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write('Adjusted PLO Coordinates:                                                                            \n')
        f.write('                                                                         LATITUDE LONGITUDE   HEIGHT \n')
        f.write(' CODE FFF STATION               LATITUDE         LONGITUDE  ORTH-HEIGHT   STD DEV   STD DEV   STD DEV\n')      
        f.write(' ---- --- ------------ ----------------- ----------------- ------------ --------- --------- ---------\n')      
        for s in adj_stns.values():
            f.write(' PLO  {:3} {:12} S {:>15} E{:>16} {:>12} {:>9} {:>9} {:>9}\n'.format(
                    s['CONST'].replace('C','1').replace('F','0'),
                    s['STATION'],
                    dd2hms(-1*float(s['LATITUDE']),5),
                    dd2hms(float(s['LONGITUDE']),5),
                    s['H_ORTHO'],
                    s['SD_E'],
                    s['SD_N'],
                    s['SD_UP']))
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write('Geoid Values: \n')
        f.write(' CODE     STATION         N/S DEFLECTION    E/W DEFLECTION   UNDULATION  \n')
        f.write(' ----     ------------ ----------------- ----------------- ------------  \n')
        for s in adj_stns.values():
            f.write(' GEOI     {:12}    0  0   0.00000 -  0  0   0.00000 {:12.3f} m\n'.format(
                    s['STATION'],
                    float(s['H_ELLIPSE'])-float(s['H_ORTHO'])))
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')                
        f.write('\n')
        f.write('                                                                        OBSERV   RESIDUAL                  \n')
        f.write('TYPE AT           FROM         TO                 OBSERVATION RESIDUAL  STD DEV  STD DEV  STD RES      PPM \n')
        f.write('---- ------------ ------------ ------------ ----------------- -------- -------- -------- -------- -------- \n')                
        PrevTyp=''
        for o in adj_obs.values():
            Grp=''; sAt=''; ppm=''
            sFrm=o['STATION_1']
            sTo=o['STATION_2']
            if o['M']=='H': Grp='OHGT'; sAt=o['STATION_1']; sFrm=''
            if o['M']=='L': Grp='OHDF'
            if o['M']=='B': Grp='GAZI'
            if o['M']=='K': Grp='AZIM'
            if o['M']=='S': Grp='DIST'
            if o['M']=='D': Grp='DIR'
            if o['M']=='G': Grp='D'+o['C']+'CT';
            if PrevTyp!=o['M']:f.write('\n')
            if o['C']=='X':f.write('GROUP:\n')
            f.write('{:4} {:12} {:12} {:12} {:>17} {:>8} {:>8} {:>8} {:>8} {:>8}\n'.format(
                    Grp,
                    sAt,
                    sFrm,
                    sTo,
                    o['MEASURED'],
                    o['CORRECTION'],
                    o['MEAS_SD'],
                    o['CORR_SD'],
                    o['N_STAT'],
                    ppm))
            if o['OUTLIER']=='*':f.write('                                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n')
            if o['C']=='Z':f.write('\n')
            PrevTyp=o['M']
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write('Residuals (critical value = 3.042):                                            \n')
        f.write('|                                                                             |\n')
        f.write('|                   S T A T I S T I C S     S U M M A R Y                     |\n')
        f.write('|                                                                             |\n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('|                                     |                                       |\n')
        f.write('|     Residual Critical Value Type    |                Tau Max                |\n')
        f.write('|     Residual Critical Value         |                 0.0000                |\n')
        f.write('|     Number of Flagged Residuals     |{:23d}                |\n'.format(sum(1 for i in adj_obs.values() if i['OUTLIER']=='*')))
        f.write('|     Convergence Criterion           |{:>23}                |\n'.format(adj_stats['ITERATION_THRESHOLD']))
        f.write('|     Final Iteration Counter Value   |                      0                |\n')
        f.write('|     Confidence Level Used           |{:>23}                |\n'.format(adj_stats['TEST_CONFIDENCE_INTERVAL']))
        f.write('|     Estimated Variance Factor       |{:>23}                |\n'.format(adj_stats['RIGOROUS_SIGMA_ZERO']))
        f.write('|     Number of Degrees of Freedom    |{:>23}                |\n'.format(adj_stats['DEGREES_OF_FREEDOM']))
        f.write('|                                     |                                       |\n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('|                                                                             |\n')
        f.write('|                  Chi-Square Test on the Variance Factor:                    |\n')
        f.write('|                                                                             |\n')
        f.write('|'+adj_stats['CHI_SQUARE_TEST__950P'][:30].strip().center(77)+'|\n')
        f.write('|                                                                             |\n')
        if adj_stats['CHI_SQUARE_TEST__950P'].find('FAILED')!=-1:
            f.write('|                    ********   THE TEST FAILS   ********                     |\n')
        else:
            f.write('|                   ********   THE TEST PASSES   ********                     |\n')
        f.write('|                                                                             |\n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('|                                                                             |\n')
        f.write('| NOTE:  All confidence regions were computed using the following factors:    |\n')
        f.write('|        ----------------------------------------------------------------     |\n')
        f.write('|        Variance factor used      ={:>16}                          |\n'.format(adj_stats['RIGOROUS_SIGMA_ZERO']))
        f.write('|        1-D expansion factor      =          1.9600                          |\n')
        f.write('|        2-D expansion factor      =          2.4477                          |\n')
        f.write('|                                                                             |\n')
        f.write('|        Note that, for relative confidence regions, precisions are           |\n')
        f.write('|        computed from the ratio of the major semi-axis and the spatial       |\n')
        f.write('|        distance between the two stations.                                   |\n')
        f.write('|                                                                             |\n')
        f.write(' ----------------------------------------------------------------------------- \n')
        f.write('\n')
        f.write('================================================================================\n')
        f.write(adj_stats['SEGMENTATION_FILE'][:-4].center(80)+'\n')
        f.write('Dynadjust, '+adj_stats['BUILD'].ljust(40)+adj_stats['REFERENCE_FRAME']+'\n')
        f.write('================================================================================\n')
        f.write('2-D and 1-D Station Confidence Regions (95.000 and 95.000 percent):             \n')
        f.write('STATION            MAJOR SEMI-AXIS  AZ     MINOR SEMI-AXIS             VERTICAL \n')
        f.write('------------ --------------------- --- ------------------- -------------------- \n')
        for s in apu_stns.values():
            f.write('{:12} {:>21} {:3d} {:>19} {:>20}\n'.format(
                    s['STATION'],
                    s['SEMI_MAJOR'],
                    int(float(s['ORIENTATION'])),
                    s['SEMI_MINOR'],
                    s['VT_POSU']))        
       
        f.write(adj_stats['FILE_CREATED']+'\n\n')
        