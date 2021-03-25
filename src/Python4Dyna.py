import os, subprocess
import shutil
import datetime as dt
import numpy as np
from math import sqrt, sin, cos, tan, radians, degrees
from geodepy.convert import hp2dec, llh2xyz, xyz2llh, dec2hp
from geodepy.geodesy import vincinv
from geodepy.constants import grs80
from datetime import datetime

class DnaStation:
    def __init__(self):
        self.name = ''
        self.Const = ''
        self.Type = ''
        self.XAxis = ''
        self.YAxis = ''
        self.Height = ''
        self.Desc = ''
        self.aAxis = 0
        self.bAxis = 0
        self.ErrAz = 0


class AdditionalInfoStn:
    def __init__(self):
        self.Ges_Name = ''
        self.Ges_Lat = ''
        self.Ges_Lng = ''
        self.Ges_Ht = ''
        self.Ges_Ht_Acc = ''
        self.Ges_Ht_Mth = ''
        self.Ges_Hz_Ord = ''
        self.Hz_Coord_Mth = ''
        self.rel_hz_acc = ''
        self.CE = ''
        self.NonGSNumber = ''
        self.SelectPoint = 'true'
        self.SelectRL = 'true'


class AdditionalInfoMsr:
    def __init__(self):
        # Additional information is included as a comment in the DynaML file.
        # This can be used for database import
        self.StartDateTime = dt.datetime(1994, 1, 1, 00, 00, 00)
        self.Duration = dt.datetime(1966, 1, 1, 00, 00, 00)
        self.TimeStatus = ''
        self.EphemerisType = ''
        self.AtReceiver = ''
        self.ToReceiver = ''
        self.FrequencyMode = ''
        self.SurveyTechnique = 'SLEV'
        self.Solution = ''
        self.EpochInterval = ''
        self.Class = 'LC'
        self.LevelDistance = '0.01'
        self.InstrumentModel = ''
        self.Derivation = 'MEAS'
        self.NonGSNumber = ''


class DnaMeasurement:
    def __init__(self):
        self.type = ''
        self.vscale = '1'
        self.pscale = '1'
        self.lscale = '1'
        self.hscale = '1'
        self.epoch = ''
        self.first = ''
        self.second = ''
        self.third = ''
        self.stddev = ''
        self.total = ''
        self.instheight = 0
        self.targheight = 0
        self.targets = []
        self.values = []
        self.targetstddevs = []
        self.dx = ''
        self.dy = ''
        self.dz = ''
        self.id = 0
        self.Vs = np.zeros([3, 3])

    def add_targets(self, target):
        self.targets.append(target)

    def add_values(self, value):
        self.values.append(value)

    def add_targetstddevs(self, targetstddev):
        self.targetstddevs.append(targetstddev)


class DeviceHeight:
    def __init__(self):
        # Device Height might be the height of instrument at a point
        # or Height of target
        self.StnName = []
        self.RefHeight = []


def add_device_height(self, stn, hgt):
    self.StnName.append(stn)
    self.RefHeight.append(hgt)


def hms2hp(hms_ang):
    # Input: HH MM SS.ssss used by Geolab
    # Output: HH.MMSSSsssss used by DynAdjust
    sign = 1
    if hms_ang.find('S') != -1 or hms_ang.find('-') != -1:
        sign = -1
    while hms_ang.find('  ') != -1:
        hms_ang = hms_ang.replace('  ', ' ')
    hms_ang = hms_ang.replace('S', '')
    hms_ang = hms_ang.replace('E', '')
    hms_ang = hms_ang.replace('.', ' ')
    ang_part = hms_ang.split()
    if len(ang_part)>1:
        ang_part[0] = str(sign * abs(int(ang_part[0])))
        ang_part[1] = "%02d" % float(ang_part[1])
        ang_part[2] = "%02d" % float(ang_part[2])
        return ang_part[0] + '.' + ang_part[1] + ang_part[2] + ang_part[3]


def dd2hms(dd_ang,ds):
    #Input: DD.ddddd used by Geolab
    #Output: HH MM SS.sssss used by DynAdjust
    hp_ang = dec2hp(dd_ang)
    h=int(hp_ang)
    dec_m = (hp_ang-h)*100
    m=int(dec_m)
    s=(dec_m-m)*100
    return '{:} {:>2} {:>3}'.format(h, m, round(s,ds))


def find_job_num(strg):
    # search a string for 8 consecutive numbers, this is probably the Job Number
    job_num = ''
    i = 0
    while i+7 != len(strg):
        if strg[i:i + 8].isnumeric():
            job_num = strg[i:i + 8]
        i = i+1
    return job_num


def stn_xml_str(stn, stn_rec):
    # Output: String for printing to xml that is one complete station
    xml_str = '''<DnaStation>
    <Name>''' + stn.Name + '''</Name>
    <Constraints>''' + stn.Const + '''</Constraints>
    <Type>''' + stn.Type + '''</Type>
    <StationCoord>
    <Name>''' + stn.Name + '''</Name>
    <XAxis>''' + str(stn.XAxis) + '''</XAxis>
    <YAxis>''' + str(stn.YAxis) + '''</YAxis>
    <Height>''' + stn.Height + '''</Height>
    </StationCoord>
    <Description>''' + stn.Desc + '''</Description>
    <!--AdditionalInfoStn>
    <HorizCoordMethod>''' + stn_rec.Hz_Coord_Mth + '''</HorizCoordMethod>
    <RelativeHorizAccuracy>''' + stn_rec.rel_hz_acc + '''</RelativeHorizAccuracy>
    <NonGSNumber>''' + stn_rec.NonGSNumber + '''</NonGSNumber>
    <SelectPoint>''' + stn_rec.SelectPoint + '''</SelectPoint>
    <SelectRL>''' + stn_rec.SelectRL + '''</SelectRL>
    </AdditionalInfoStn-->
    </DnaStation>\n'''
    return xml_str


def msr_xml_str(msr, cntrl):
    # Output: xml string for printing to file.
    # Caters for type G, D, S, B, D, L, H
    xml_str = '<DnaMeasurement>\n'+ \
    '<Type>' + msr.type + '</Type>\n'+ \
    '<Ignore/>\n'
    if msr.type == 'G':
        xml_str = (xml_str + '<ReferenceFrame>' + 
                   gnss_date_2_ref(cntrl.StartDateTime) +
                   '</ReferenceFrame>\n'+ 
                   '<Epoch>' + 
                   cntrl.StartDateTime.strftime('%d.%m.%Y') + 
                   '</Epoch>\n' + 
                   '<Vscale>{0:.1f}</Vscale>\n'.format(float(msr.vscale)) 
                   + '<Pscale>{0:.1f}</Pscale>\n'.format(float(msr.pscale)) 
                   + '<Lscale>{0:.1f}</Lscale>\n'.format(float(msr.lscale)) 
                   + '<Hscale>{0:.1f}</Hscale>\n'.format(float(msr.hscale)))
    xml_str = xml_str + '<First>' + msr.first + '</First>\n'
    if msr.second != '':
        xml_str = xml_str + '<Second>' + msr.second + '</Second>\n'
    if msr.type != 'G' and msr.type != 'D':
        xml_str = (xml_str + '<Value>' + msr.values[0] + '</Value>\n' + 
             '<StdDev>' + msr.stddev + '</StdDev>\n')
    if msr.type == 'G':
        xml_str = (xml_str + '<GPSBaseline>\n'+
        '<X>' + str(msr.dx) + '</X>\n'+
        '<Y>' + str(msr.dy) + '</Y>\n'+
        '<Z>' + str(msr.dz) + '</Z>\n'+
        '<!--MeasurementID>' + str(msr.id) + '</MeasurementID-->\n'+
        '<SigmaXX>' + str(msr.Vs[0, 0]) + '</SigmaXX>\n'+
        '<SigmaXY>' + str(msr.Vs[0, 1]) + '</SigmaXY>\n'+
        '<SigmaXZ>' + str(msr.Vs[0, 2]) + '</SigmaXZ>\n'+
        '<SigmaYY>' + str(msr.Vs[1, 1]) + '</SigmaYY>\n'+
        '<SigmaYZ>' + str(msr.Vs[1, 2]) + '</SigmaYZ>\n'+
        '<SigmaZZ>' + str(msr.Vs[2,2]) + '</SigmaZZ>\n'+
        '</GPSBaseline>\n '+
        '<!--AdditionalInfoMsrG>\n'+ 
        '<StartDateTime>' + 
        cntrl.StartDateTime.strftime('%Y-%m-%dT%H:%M:%S') + 
        '</StartDateTime>\n'+ 
        '<Duration>' + 
        'P' + str(cntrl.Duration.year - 1900) + 
        'Y' + str(cntrl.Duration.month - 1) + 
        'M' + str(cntrl.Duration.day - 1) + 
        'DT' + cntrl.Duration.strftime('%HH%MM%SS') + 
        '</Duration>\n'+ 
        '<TimeStatus>' + cntrl.TimeStatus + '</TimeStatus>\n'+ 
        '<EphemerisType>' + cntrl.EphemerisType + '</EphemerisType>\n'+ 
        '<AtReceiver>' + cntrl.AtReceiver + '</AtReceiver>\n'+ 
        '<ToReceiver>' + cntrl.ToReceiver + '</ToReceiver>\n'+ 
        '<FrequencyMode>' + cntrl.FrequencyMode + '</FrequencyMode>\n'+ 
        '<SurveyTechnique>' + cntrl.SurveyTechnique + '</SurveyTechnique>\n'+ 
        '<Solution>' + cntrl.Solution + '</Solution>\n'+ 
        '<EpochInterval>' + str(cntrl.EpochInterval) + '</EpochInterval>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + str(cntrl.NonGSNumber) + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrG-->\n')
    
    if msr.type == 'L':
        xml_str = (xml_str + '<!--AdditionalInfoMsrL>\n'+ 
        '<SurveyTechnique>' + cntrl.SurveyTechnique + '</SurveyTechnique>\n'+ 
        '<LevelDistance>' + cntrl.LevelDistance + '</LevelDistance>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrL-->\n')
    
    if msr.type == 'S':
        xml_str = (xml_str + '<InstHeight>' + str(msr.instheight) + '</InstHeight>\n'+ 
        '<TargHeight>' + str(msr.targheight) + '</TargHeight>\n'+ 
        '<!--AdditionalInfoMsrS>\n'+ 
        '<InstrumentModel>' + cntrl.InstrumentModel + '</InstrumentModel>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrS-->\n')

    if msr.type == 'D':
        xml_str = (xml_str + '<Value>' + msr.values[0] + '</Value>\n'+ 
        '<StdDev>' + msr.targetstddevs[0] + '</StdDev>\n'+ 
        '<Total>' + str(msr.total - 1) + '</Total>\n')
        obs_num = 1
        while obs_num < msr.total:
            xml_str = (xml_str + '<Directions>\n'+ 
            '<Ignore/>\n'+ 
            '<Target>' + msr.targets[obs_num] + '</Target>\n'+ 
            '<Value>' + msr.values[obs_num] + '</Value>\n'+ 
            '<StdDev>' + msr.targetstddevs[obs_num] + '</StdDev>\n'+ 
            '</Directions>\n')
            obs_num = obs_num + 1
        xml_str = (xml_str + '<!--AdditionalInfoMsrD>\n'+ 
        '<InstrumentModel>' + cntrl.InstrumentModel + '</InstrumentModel>\n'+ 
        '<ObsDate>' + cntrl.StartDateTime.strftime('%Y-%m-%d') + '</ObsDate>\n'+ 
        '<Derivation>' + cntrl.Derivation + '</Derivation>\n'+ 
        '<Class>' + cntrl.Class + '</Class>\n'+ 
        '<NonGSNumber>' + cntrl.NonGSNumber + '</NonGSNumber>\n'+ 
        '</AdditionalInfoMsrD-->\n')
    xml_str = xml_str+'<Source></Source>\n'
    if msr.type != 'G':
        xml_str = xml_str+'<!--MeasurementID>' + str(msr.id) + '</MeasurementID-->\n'
    xml_str = xml_str+'</DnaMeasurement>\n'
    
    return xml_str


def get_rl(p,l):
    rlat = radians(p)
    rlng = radians(l)
    rl = np.zeros([3, 3])
    rl[0, 0] = -sin(rlng)
    rl[0, 1] = -sin(rlat)*cos(rlng)
    rl[0, 2] = cos(rlat)*cos(rlng)
    rl[1, 0] = cos(rlng)
    rl[1, 1] = -sin(rlat)*sin(rlng)
    rl[1, 2] = cos(rlat)*sin(rlng)
    rl[2, 1] = cos(rlat)
    rl[2, 1] = sin(rlat)
    return rl


def err_ellip_2_ycluster(stn):
    # Input: Supply a station with coordinates
    #        and error ellipse for coordinate uncertainty
    # Output: xml string for  point cluster (Y-type observation)
    x, y, z = llh2xyz(float(stn.XAxis), float(stn.YAxis), float(stn.Height))
    
    a = stn.aAxis / 2.44774683068
    b = stn.bAxis / 2.44774683068
    az = 90 - stn.ErrAz
    
    r_az = radians(az)
    rl = get_rl(float(stn.XAxis),float(stn.YAxis))
    
    i_a = np.zeros([3, 3])
    i_a[0, 0] = (cos(r_az)*cos(r_az)*a*a)+(b*b*sin(r_az)*sin(r_az))
    i_a[0, 1] = (a*a-b*b)*cos(r_az)*sin(r_az)
    i_a[1, 0] = i_a[0, 1]
    i_a[1, 1] = (a*a*sin(r_az)*sin(r_az))+(b*b*cos(r_az)*cos(r_az))
    i_a[2, 2] = 0.000001
    
    wt = np.matmul(np.matmul(rl, i_a), rl.transpose())
    
    xml_str = ('<DnaMeasurement>\n'+ 
        '<Type>Y</Type>\n'+ 
        '<Ignore/>\n'+ 
        '<ReferenceFrame>GDA2020</ReferenceFrame>\n'+ 
        '<Epoch>01.01.2020</Epoch>\n'+ 
        '<Vscale>1.000</Vscale>\n'+ 
        '<Pscale>1.000</Pscale>\n'+ 
        '<Lscale>1.000</Lscale>\n'+ 
        '<Hscale>1.000</Hscale>\n'+ 
        '<Coords>XYZ</Coords>\n'+ 
        '<Total>1</Total>\n'+ 
        '<First>' + stn.Name + '</First>\n'+ 
        '<Clusterpoint>\n'+ 
        '<X>'+str(x)+'</X>\n'+ 
        '<Y>'+str(y)+'</Y>\n'+ 
        '<Z>'+str(z)+'</Z>\n'+ 
        '<SigmaXX>'+str(wt[0, 0])+'</SigmaXX>\n'+ 
        '<SigmaXY>'+str(wt[0, 1])+'</SigmaXY>\n'+ 
        '<SigmaXZ>'+str(wt[0, 2])+'</SigmaXZ>\n'+ 
        '<SigmaYY>'+str(wt[1, 1])+'</SigmaYY>\n'+ 
        '<SigmaYZ>'+str(wt[1, 2])+'</SigmaYZ>\n'+ 
        '<SigmaZZ>'+str(wt[2, 2])+'</SigmaZZ>\n'+ 
        '</Clusterpoint>\n'+ 
        '</DnaMeasurement>\n')
    
    return xml_str


def stn_header(datum):
    xml_str='<?xml version="1.0" encoding="utf-8"?>\n'
    epoch='01.01.2020'
    if datum.find('09')!=-1: epoch='01.01.1994'
    xml_str=(xml_str+'<DnaXmlFormat type="Station File" referenceframe="'+datum+
             '" epoch="'+epoch+
             '" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'+
             ' xsi:noNamespaceSchemaLocation="DynaML.xsd">\n')
    return xml_str


def msr_header():
    xml_str = '<?xml version="1.0" encoding="utf-8"?>\n'+\
               '<DnaXmlFormat type="Measurement File" '+\
               'referenceframe="GDA94" epoch="01.01.1994" '+\
               'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '+\
               'xsi:noNamespaceSchemaLocation="DynaML.xsd">\n'
    return xml_str


def dml_footer():
    xml_str = '</DnaXmlFormat>\n'
    return xml_str


def gnss_date_2_ref(obs_date):
    # Use the date of GNSS baseline observation
    # to determine the reference frame used by broadcast ephemeris
    if dt.datetime(1900, 1, 1) <= obs_date < dt.datetime(1994, 1, 2):
        ref = 'ITRF1991'
    if dt.datetime(1994, 1, 2) <= obs_date < dt.datetime(1995, 1, 1):
        ref = 'ITRF1992'
    if dt.datetime(1995, 1, 1) <= obs_date < dt.datetime(1996, 6, 30):
        ref = 'ITRF1993'
    if dt.datetime(1996, 6, 30) <= obs_date < dt.datetime(1998, 3, 1):
        ref = 'ITRF1994'
    if dt.datetime(1998, 3, 1) <= obs_date < dt.datetime(1999, 8, 1):
        ref = 'ITRF1996'
    if dt.datetime(1999, 8, 1) <= obs_date < dt.datetime(2001, 12, 2):
        ref = 'ITRF1997'
    if dt.datetime(2001, 12, 2) <= obs_date < dt.datetime(2006, 11, 5):
        ref = 'ITRF2000'
    if dt.datetime(2006, 11, 5) <= obs_date < dt.datetime(2011, 4, 17):
        ref = 'ITRF2005'
    if dt.datetime(2011, 4, 17) <= obs_date < dt.datetime(2017, 1, 29):
        ref = 'ITRF2008'
    if obs_date >= dt.datetime(2017, 1, 29):
        ref = 'ITRF2014'
    return ref


def adj2Dicts(adj_file):
    #############################################################################
    ##     Input: adj_file.adj                                                 ##
    ##     Output: adj_file.adj (The original file is reformated to include)   ##
    ##                            a blank line after each GNSS basline)        ##
    ##     Return: coords, xyz_bases (Station Coordintes (dd) and baselines    ##
    ##                            returned as a dictionary)                    ##
    #############################################################################
    with open(adj_file) as f:
        first_line = f.readline()
    if first_line.startswith('This DynAdjust file has been re-formatted by program DynaGEOL2020.'):
        print(adj_file +' has already been re-formatted by program DynaGEOL2020.\n')
    else:
        try: os.remove('$'+adj_file)
        except: ''
        os.rename( adj_file, '$' + adj_file)
        f_in = open('$'+adj_file, 'r')
        f_out = open(adj_file, 'w')
        part_G_cnt=1; base_cnt=0; coords = {}; xyz_bases = {}
        E=''; N=''; z=''; P=''; L=''; H=''; h=''
        read_msr=False; read_coord=-1; angular_msr_type=''
        f_out.write('This GEOLAB file has been re-formatted by program DynaGEOL2020.\n')
        f_out.write('DYNADJUST ADJUSTMENT OUTPUT FILE\n')
        for linestr in f_in.readlines():
            f_out.write(linestr)
            if (linestr.find('Command line arguments:')!=-1 and 
                linestr.find('--angular-msr-type 1')!=-1):angular_msr_type='dd'
            if linestr.find('Station coordinate types:')!=-1: coord_types = linestr[35:].strip()
            if linestr.find('Adjusted Measurements')!=-1: read_msr=True
            if read_msr and linestr[0:2]=='G ':
                a_obs=linestr[67:].split()
                if part_G_cnt==1: dX=float(a_obs[0]); sdev_x=float(a_obs[3]); nstat_x=float(a_obs[6])
                if part_G_cnt==2: dY=float(a_obs[0]); sdev_y=float(a_obs[3]); nstat_y=float(a_obs[6])
                if part_G_cnt==3: 
                    dZ=float(a_obs[0]); sdev_z=float(a_obs[3]); nstat_z=float(a_obs[6])
                    xyz_bases[base_cnt] = {'first':linestr[2:22].strip(),'second':linestr[22:42].strip(),
                             'dX':dX, 'dY':dY, 'dZ':dZ,
                             'sdev_x':sdev_x, 'sdev_y':sdev_y, 'sdev_z':sdev_z,
                             'nstat_x':nstat_x, 'nstat_y':nstat_y, 'nstat_z':nstat_z}
                    f_out.write('\n')
                    base_cnt+=1
                    part_G_cnt=0
                part_G_cnt+=1
                
            if linestr.find('Adjusted Coordinates')!=-1:
                read_msr=False; read_coord=0
            if read_coord>4 and linestr!='\n':
                stn = linestr[0:20].strip()
                results = linestr[25:].split()
                r_count = 0
                for ct in coord_types:
                    if ct == 'E': E = float(results[r_count])
                    if ct == 'N': N = float(results[r_count])
                    if ct == 'z': z = int(results[r_count])
                    if ct == 'P':
                        if angular_msr_type!='dd':
                            P = hp2dec(float(results[r_count]))
                        else:
                            P = float(results[r_count])
                    if ct == 'L':
                        if angular_msr_type!='dd':
                            L = hp2dec(float(results[r_count]))
                        else:
                            L = float(results[r_count])
                    if ct == 'H': H = float(results[r_count])
                    if ct == 'h': h = float(results[r_count])
                    r_count += 1
                        
                coords[str(stn)] = {'E': E,'N': N,'Z': z,'LAT': P,'LON': L,'OHGT': H,'EHGT': h }
            if read_coord>=0:read_coord+=1
        f_in.close()
        f_out.close()
        os.remove('$'+adj_file)
        return coords, xyz_bases
    

def gnss_xyz2dah(coords, xyz_bases):
    dah_bases = {}
    print('   Calculating baseline transformations ')
    for b in xyz_bases:
        f_x, f_y, f_z = llh2xyz(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],coords[xyz_bases[b]['first']]['EHGT'])
        Calc_point_x = f_x+xyz_bases[b]['dX']
        Calc_point_y = f_y+xyz_bases[b]['dY']
        Calc_point_z = f_z+xyz_bases[b]['dZ']
        p, l, h = xyz2llh(Calc_point_x, Calc_point_y, Calc_point_z)
        e_dist, b_az, b_rev_az = vincinv(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],p, l)
        adj_dist, adj_az, b_rev_az = vincinv(coords[xyz_bases[b]['first']]['LAT'],coords[xyz_bases[b]['first']]['LON'],coords[xyz_bases[b]['second']]['LAT'],coords[xyz_bases[b]['second']]['LON'])
        e_dist_res = adj_dist - e_dist
        b_az_res = adj_az - b_az
        e_ht_diff = h - coords[xyz_bases[b]['first']]['EHGT']
        n1=coords[xyz_bases[b]['first']]['OHGT']-coords[xyz_bases[b]['first']]['EHGT']
        n2=coords[xyz_bases[b]['second']]['OHGT']-coords[xyz_bases[b]['second']]['EHGT']
        o_ht_diff = -1*n1 + e_ht_diff + n2
        o_ht_diff_res =(coords[xyz_bases[b]['second']]['OHGT']-coords[xyz_bases[b]['first']]['OHGT']) - o_ht_diff
        
        #transform the baseline to distance, azimuth, up
        #translated code from http://peas13.dli.wa.gov.au/svn-intsvc/Delphi/GesmarComponents/trunk/Geodetic.pas
        m_lat = (coords[xyz_bases[b]['first']]['LAT']+p)/2
        m_x =(f_x+Calc_point_x)/2
        m_y =(f_y+Calc_point_y)/2
        
        dR       = sqrt(m_x*m_x + m_y*m_y)
        dE2      = grs80.ecc1sq
        dQ       = radians(m_lat)
        sinQ     = sin(dQ)
        dTemp    = sqrt(1.0 - dE2 * sinQ * sinQ)
        dN       = grs80.semimaj / dTemp
        dM       = dN * ((1.0 - dE2) / (dTemp * dTemp))
        sinQ     = sin(dQ)
        cosQ     = cos(dQ)
        secQ     = 1.0/cosQ
        dF       = grs80.f
        dE2      = 2 * dF -(dF * dF)
        dA       = m_x * tan(dQ)/(dR*dR*(dE2 - secQ*secQ))
        dB       = 1.0/(dR*secQ*secQ - dE2*dN*cosQ)
        dC       = dR*sinQ/(cosQ*cosQ)-dN*dE2*sinQ*cosQ/(1.0-dE2*sinQ*sinQ)
        
        JMatrix      = np.zeros([3,3])
        JMatrix[0,0] = dM*dA
        JMatrix[0,1] = dM*(m_y*dA/m_x)
        JMatrix[0,2] = dM*dB
        JMatrix[1,0] = -(dN*cosQ)*m_y/(dR*dR)
        JMatrix[1,1] = (dN*cosQ)*m_x/(dR*dR)
        JMatrix[1,2] = 0.0
        JMatrix[2,0] = m_x/(dR*cosQ) + dA*dC
        JMatrix[2,1] = m_y/(dR*cosQ) + m_y*dA*dC/m_x
        JMatrix[2,2] = dB*dC
        
        b_var_matrix = np.zeros([3,3])
        b_var_matrix[0,0] = xyz_bases[b]['sdev_x']**2
        b_var_matrix[1,1] = xyz_bases[b]['sdev_y']**2
        b_var_matrix[2,2] = xyz_bases[b]['sdev_z']**2
                
        b_nst_matrix = np.zeros([3,3])
        b_nst_matrix[0,0] = xyz_bases[b]['nstat_x']**2
        b_nst_matrix[1,1] = xyz_bases[b]['nstat_y']**2
        b_nst_matrix[2,2] = xyz_bases[b]['nstat_z']**2

        cosAz        = cos(radians(b_az))
        sinAz        = sin(radians(b_az))
        AMatrix      = np.zeros([3,3])
        AMatrix[0,0] = cosAz
        AMatrix[0,1] = sinAz
        if e_dist==0:
            AMatrix[1,0] = degrees(sinAz/0.00001)
            AMatrix[1,1] = -degrees(cosAz/0.00001)
        else:
            AMatrix[1,0] = degrees(sinAz/e_dist)
            AMatrix[1,1] = -degrees(cosAz/e_dist)

        GMatrix = np.matmul(np.matmul(JMatrix,b_var_matrix),JMatrix.transpose())
        dah_Matrix = np.matmul(np.matmul(AMatrix,GMatrix),AMatrix.transpose())

        n_GMatrix = np.matmul(np.matmul(JMatrix,b_nst_matrix),JMatrix.transpose())
        nstat_Matrix = np.matmul(np.matmul(AMatrix,n_GMatrix),AMatrix.transpose())

        dah_bases[b] = {'first':xyz_bases[b]['first'],'second':xyz_bases[b]['second'],
                        'dX':xyz_bases[b]['dX'], 'dY':xyz_bases[b]['dY'], 'dZ':xyz_bases[b]['dZ'],
                        'e_dist':e_dist, 'e_dist_res':e_dist_res,
                        'e_dist_sdev':sqrt(abs(dah_Matrix[0,0])),
                        'e_dist_nstat':sqrt(abs(nstat_Matrix[0,0])),
                        'b_az':b_az, 'b_az_res':b_az_res, 
                        'b_az_sdev':sqrt(abs(dah_Matrix[1,1])),
                        'b_az_nstat':sqrt(abs(sin(radians(nstat_Matrix[1,1]))*e_dist)),
                        'o_ht_diff':o_ht_diff, 'e_ht_diff':e_ht_diff, 'o_ht_diff_res':o_ht_diff_res,
                        'e_ht_diff_sdev':sqrt(abs(GMatrix[2,2])),
                        'e_ht_diff_nstat':sqrt(abs(n_GMatrix[2,2]))}
    return coords, dah_bases


def Create_DynaBAS(out_file,bases):
    print('   writing: '+out_file+'.BAS')
    f_bas = open(out_file+'.BAS', 'w')
    f_bas.write( DynaBAS_header(out_file))
    f_bas.write('                                       Ellipsoidal\n')
    f_bas.write('BASE                                   Distance                Standard\n')
    f_bas.write('LINE        From         To            Observation   Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  -----------  ---------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>12} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'],
                '{:.3f}'.format(bases[b]['e_dist']), '{:.3f}'.format(bases[b]['e_dist_res']),
                '{:.3f}'.format(bases[b]['e_dist_sdev']), '{:.2f}'.format(bases[b]['e_dist_nstat'])) + '\n')
    
    f_bas.write('\n')
    f_bas.write( DynaBAS_header(out_file))
    f_bas.write('\n')
    f_bas.write('BASE                                   Azimuth                  Standard\n')
    f_bas.write('LINE        From         To            Observation    Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  ------------  ---------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>13s} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'], 
                dd2hms(bases[b]['b_az'],2), '{:.2f}'.format(bases[b]['b_az_res']*3600),
                '{:.2f}'.format(bases[b]['b_az_sdev']*3600),'{:.2f}'.format(bases[b]['b_az_nstat'])) + '\n')
    
    f_bas.write('\n')
    f_bas.write( DynaBAS_header(out_file))
    f_bas.write('                                       Height Difference\n')
    f_bas.write('                                       Observation\n')
    f_bas.write('BASE                                   -----------------------             Standard\n')
    f_bas.write('LINE        From         To            Orthometric  Spheroidal   Residual  Deviation     N-Stat\n')
    f_bas.write('----------- ------------ ------------  ----------- -----------   --------  ---------   --------\n')
    for b in bases:
        f_bas.write('{:11s} {:12s} {:12s} {:>12} {:>11} {:>10} {:>10} {:>10}'.format(
                str(b), bases[b]['first'], bases[b]['second'], '{:.3f}'.format(
                        bases[b]['o_ht_diff']), '{:.3f}'.format(bases[b]['e_ht_diff']), '{:.3f}'.format(bases[b]['o_ht_diff_res']),
                        '{:.3f}'.format(bases[b]['e_ht_diff_sdev']),'{:.2f}'.format(bases[b]['e_ht_diff_nstat'])) + '\n')
    
    f_bas.close()


def DynaBAS_header(out_file):
    header_str='================================================================================\n'
    header_str=header_str + out_file.center(80)+'\n'
    header_str=header_str + 'GRS80\n'
    header_str=header_str + '================================================================================\n'
    return header_str


def LSA(network, datum, geoid):
    xsd_file='\\\\dli\public\\Business Unit\\Operations\\Registrations\\CadastralSubdivisions\\Geodetic\\AssetsPCsAndSoftware\\Software\\DynAdjust\\DynaML.xsd'
    if os.path.exists(xsd_file): shutil.copy(xsd_file, os.getcwd())

    ####################################
    ### Run the DynAdjust Adjustment ###
    ####################################
    print('  Adjusting: ' + network)
    subprocess.run("import -n " + network + " " + 
                   network + ".msr.xml "+ 
                   network + ".stn.xml "+
                   "--flag-unused-stations "+
                   "--remove-ignored-msr -r " + datum)
    subprocess.run("reftran " + network + " -r " + datum)
    subprocess.run("geoid " + network + 
                   " -g \""+geoid+"\" --convert-stn-hts")
    subprocess.run("segment " + network + 
                   " --min-inner-stns 500 "+
                   "--max-block-stns 500")
    subprocess.run("adjust " + network + 
                   " --staged "+
                   "--create-stage-files "+
                   "--output-adj-msr "+
                   "--output-pos-uncertainty "+
                   "--output-adj-gnss-units 0 "+
                   "--max-iterations 20 "+
                   "--free-stn-sd 5 "+
                   "--iteration-threshold 0.0005 "+
                   "--stn-coord-types ENzPLHhXYZ")

        
def lg_xml(file):
    f_in = open(file, 'r')
    f_out = open(file+'.lg.xml', 'w')
    for ln in f_in.readlines():
        ln=ln.replace('<!--AdditionalInfo','<AdditionalInfo')
        if ln.find('-->')>0 and ln.find('/AdditionalInfo')>0:
            ln=ln.replace('-->','>')
        f_out.write(ln)
    return file+'.lg.xml'


def kmlHeader(nme):
    strg = ('<?xml version="1.0" encoding="UTF-8"?>\n'
    + '<kml xmlns="http://earth.google.com/kml/2.0">\n'
    + ' <Document>\n'
    + '   <name>' + nme + '</name>\n'
    + '    <Style id="trivlineStyle">\n'
    + '        <LineStyle>\n'
    + '            <color>660000ff</color>\n'
    + '            <width>6</width>\n'
    + '        </LineStyle>\n')
    return strg + '    </Style>'


def fmt_int_date(v):
    return datetime.strftime(datetime.strptime(str(v),'%Y%m%d%H%M%S'),'%d-%m-%Y %H:%M:%S')


def MkLine(ln,s):
    strg = ('   <Placemark>\n'
    +'      <name>' + str(ln[1])+' - - ' + str(ln[2])+ '</name>\n'
    + '    <description>'+str(s[0])+', '+str(s[1])+'\nStart: ' + fmt_int_date(ln[7]) + '\nFinish: ' + fmt_int_date(ln[8]) + '</description>\n'
    + '    <styleUrl>#trivlineStyle</styleUrl>\n'
    + '    <LineString>\n'
    + '      <extrude>1</extrude>\n'
    + '      <tessellate>1</tessellate>\n'
    + '      <altitudeMode>ClampToGround</altitudeMode>\n'
    + '      <coordinates>\n'
    + '       ' +  str(hp2dec(ln[4])) + ',' +  str(hp2dec(ln[3])) + ',0 ' +  str(hp2dec(ln[6])) + ',' +  str(hp2dec(ln[5])) + ',0 \n'
    + '      </coordinates>\n'
    + '    </LineString>\n')
    return strg +'  </Placemark>'


def kmlFooter():
    return '</Document>\n</kml>'
