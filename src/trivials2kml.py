import isodate
import Python4Dyna

def create_DynaDb(conn):
    cursor = conn.cursor()
        
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJUSTMENT (
                  ID integer PRIMARY KEY, 
                  ADJUSTMENT_NAME short text NOT NULL UNIQUE);""")
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJ (
                  ID integer PRIMARY KEY, 
                  ADJUSTMENT_ID integer);""")
    cursor.execute("""CREATE TABLE IF NOT EXISTS APU (
                  ID integer PRIMARY KEY, 
                  ADJUSTMENT_ID integer);""")
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJ_MEASUREMENTS (
                  ID integer PRIMARY KEY, 
                  ADJUSTMENT_ID integer);""")
    cursor.execute("""CREATE TABLE IF NOT EXISTS X_G_Y_COMPONENTS (
                  ID integer PRIMARY KEY, 
                  OBSERVATION_ID integer);""")
    cursor.execute("""CREATE TABLE IF NOT EXISTS ADJUSTMENT (
                  ID integer PRIMARY KEY, 
                  ADJUSTMENT_NAME short text);""")
    conn.commit()


def add_tbl_clm(tbl,clms,conn):
    cursor = conn.cursor()
    e_clms=cursor.execute('PRAGMA table_info('+tbl+')').fetchall()
    for c in clms:
        add_c=False
        for e in e_clms:
            if c in e: add_c=True
        if add_c==False:
            sql='ALTER TABLE ' + tbl + ' ADD '
            sql= sql + c
            if (c =='M' 
                or c.startswith('Station') 
                or c =='C' or c =='Outlier?' 
                or c =='Const' 
                or c == 'Description'):
                sql = sql + ' short text;'
                cursor.execute(sql)
            else:
                sql = sql + ' double;'
                cursor.execute(sql)
        conn.commit()
 
       
def chg4sql(stg):
    stg=(stg.replace('*',' ')
            .replace('?',' ')
            .replace('.','')
            .replace('-','_')
            .replace('(','_')
            .replace(')',' ')
            .replace('M Station 1','M  Station 1')
            .upper())
    clms=[s.strip().replace(' ','_') for s in stg.split('  ') if s]
    return clms


def list2sql (clms):
    sql = ' ('
    v=''
    for c in clms:
        sql=sql + c + ', '
        v = v + '?,'
    sql = sql[:-2] + ') VALUES ('+v[:-1] +')'     
    return sql


def stn2db(stn,conn):
    cursor = conn.cursor()
    if stn['Description']==None:
        cntrl='|||||||||||||'.split('|')
    else:
        cntrl=stn['Description'].replace(';', '|').split('|')
    if len(cntrl)<12:cntrl='|||||||||||||'.split('|')
    
    clms =['COORD_TYPE', 'STATION_NAME', 'CONSTRAIN', 
            'LATITUDE', 'LONGITUDE', 
            'HEIGHT', 'DESC', 
            'GES_NAME', 'GES94_LATITUDE', 'GES94_LONGITUDE', 'GES_HEIGHT', 
            'HT_ACCURACY', 'HT_METHOD', 
            'HZ_ORDER', 'HZ_ACCURACY', 'HZ_METHOD', 
            'CE', 'A_AXIS', 'B_AXIS', 'ERR_AZ']
    sql =list2sql(clms)
    cursor.execute('INSERT INTO STATIONS '+sql,
            [stn['Type'], stn['Name'],stn['Constraints'], 
             stn['StationCoord']['XAxis'],stn['StationCoord']['YAxis'], 
             stn['StationCoord']['Height'], stn['Description'], 
             cntrl[0], cntrl[1], cntrl[2], cntrl[3],
             cntrl[4], cntrl[5],
             cntrl[6], cntrl[7], cntrl[8],
             cntrl[9], cntrl[10],cntrl[11], cntrl[12]])
    conn.commit()       


def msr2db(msr,conn):
    cursor = conn.cursor()
    cntrl={}
    if 'AdditionalInfoMsrS' in msr: cntrl=msr['AdditionalInfoMsrS']
    if 'AdditionalInfoMsrM' in msr: cntrl=msr['AdditionalInfoMsrM']
    if 'AdditionalInfoMsrL' in msr: cntrl=msr['AdditionalInfoMsrL']
    if 'AdditionalInfoMsrK' in msr: cntrl=msr['AdditionalInfoMsrK']
    if 'AdditionalInfoMsrE' in msr: cntrl=msr['AdditionalInfoMsrE']
    if 'AdditionalInfoMsrD' in msr: cntrl=msr['AdditionalInfoMsrD']
    if 'AdditionalInfoMsrB' in msr: cntrl=msr['AdditionalInfoMsrB']
    if 'AdditionalInfoMsrG' in msr: cntrl=msr['AdditionalInfoMsrG']
    
    m = DnaMeasurement()
    for k in msr.keys():
        if k.lower() in vars(m): vars(m)[k.lower()]=msr[k]
    if 'Value' in msr:m.add_values(msr['Value'])
    if 'StdDev' in msr: m.stddev=msr['StdDev']
    if 'MeasurementID' in msr: m.id=msr['MeasurementID']
    if 'Directions' in msr:
        m.add_targets(m.second)
        m.add_targetstddevs(m.stddev)
        for d in msr['Directions']:
            m.add_targets(d[0])
            m.add_values(d[1])
            m.add_targetstddevs(d[2])
    if 'GPSBaseline' in msr:
        m.dx=msr['GPSBaseline']['X']
        m.dy=msr['GPSBaseline']['Y']
        m.dz=msr['GPSBaseline']['Z']
        m.Vs[0,0]=msr['GPSBaseline']['SigmaXX']
        m.Vs[0,1]=msr['GPSBaseline']['SigmaXY']
        m.Vs[0,2]=msr['GPSBaseline']['SigmaXZ']
        m.Vs[1,1]=msr['GPSBaseline']['SigmaYY']
        m.Vs[1,2]=msr['GPSBaseline']['SigmaYZ']
        m.Vs[2,2]=msr['GPSBaseline']['SigmaZZ']
        #m.id=msr['GPSBaseline']['MeasurementID']
    
    c = AdditionalInfoMsr()
    for k in cntrl.keys():
        if k in vars(c): vars(c)[k]=cntrl[k]
    if 'ObsDate' in cntrl:
        c.StartDateTime=isodate.parse_date(cntrl['ObsDate'])
    if 'StartDateTime' in cntrl: 
        c.StartDateTime=isodate.parse_date(cntrl['StartDateTime'])
    
    fnsh=0; drtn=0
    if msr['Type']=='G':  
        t_zero = datetime.strptime('00:00:00', '%H:%M:%S')
        drtn = isodate.parse_duration(c.Duration)
        fnsh = c.StartDateTime + drtn
        fnsh = fnsh.strftime('%Y%m%d%H%M%S')
        drtn = (t_zero + drtn).strftime('%Y%m%d%H%M%S')
        m.add_values('')
    strt = c.StartDateTime.strftime('%Y%m%d%H%M%S')
        
    clms =['TYPE', 
        'VSCALE', ' PSCALE', ' LSCALE', ' HSCALE', 
        'FIRST', ' SECOND', ' VALUE', ' SDEV', 
        'DX', ' DY', ' DZ', 
        'VS_1_1', ' VS_1_2', ' VS_1_3', 
        'VS_2_1', ' VS_2_2', 
        'VS_3_1', 
        'StartDateTime', ' FinishDateTime', ' Duration', ' TimeStatus', 
        'EphemerisType', 
        'AtReceiver', ' ToReceiver', ' FrequencyMode', 
        'SurveyTechnique', ' LevelDistance', ' Derivation', 
        'Solution', ' EpochInterval', ' Class', ' NON_GS']
    sql =list2sql(clms)
    cursor.execute('INSERT INTO OBSERVATIONS ' + sql, 
       [m.type, 
        m.vscale, m.pscale, m.lscale, m.hscale, 
        m.first, m.second, m.values[0], m.stddev,
        m.dx, m.dy, m.dz,
        m.Vs[0,0], m.Vs[0,1], m.Vs[0,2],
        m.Vs[1,1], m.Vs[1,2],
        m.Vs[2,2],
        int(strt), int(fnsh), int(drtn), c.TimeStatus,
        c.EphemerisType,
        c.AtReceiver,c.ToReceiver,c.FrequencyMode, 
        c.SurveyTechnique, c.LevelDistance, c.Derivation,
        c.Solution,c.EpochInterval, c.Class, c.NonGSNumber])

    if m.type=='D':
        obs_id = cursor.execute("SELECT MAX(ID) FROM OBSERVATIONS").fetchall()
        d = 0
        while d !=int(msr['Total']):
            cursor.execute('''INSERT INTO DIR_TARGETS (OBSERVATIONS_ID, 
            TARGETS, VALUE, TARGETS_SDEV)  
            VALUES (?,?,?,?)''', 
               [obs_id[0][0], msr['Directions']['Target'],
               msr['Directions']['Value'], 
               msr['Directions']['StdDev']])
            d+=1
    conn.commit() 
    
def create_DynaML_db(conn):
    # Connect to Sqlite and Open a new database
    cursor = conn.cursor()
        
    cursor.execute('''CREATE TABLE IF NOT EXISTS STATIONS (
                  ID integer PRIMARY KEY, 
                  STATION_NAME short text, 
                  COORD_TYPE text, CONSTRAIN text, W_CONSTRAIN text, 
                  LATITUDE double, LONGITUDE double, 
                  HEIGHT double, E_HEIGHT double,
                  DESC text, GES_NAME text, 
                  GES94_LATITUDE double, GES94_LONGITUDE double, 
                  GES_HEIGHT double, HT_ACCURACY text, HT_METHOD text,
                  HZ_ORDER text, HZ_ACCURACY text, HZ_METHOD text, 
                  CE double, A_AXIS double, B_AXIS double, ERR_AZ double,
                  GES2020_LATITUDE double, GES2020_LONGITUDE double);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS OBSERVATIONS (
                  ID integer PRIMARY KEY, TYPE text, 
                  VSCALE double, PSCALE double, LSCALE double, HSCALE double, 
                  FIRST short text, SECOND short text, THIRD short text,
                  VALUE text, SDEV text, TOTAL integer, 
                  INST_HEIGHT double, TARG_HEIGHT double, 
                  DX double, DY double, DZ double, 
                  VS_1_1 double, VS_1_2 double, VS_1_3 double, 
                  VS_2_1 double, VS_2_2 double, 
                  VS_3_1 double,
                  StartDateTime integer,FinishDateTime integer, Duration time, 
                  TimeStatus text, EphemerisType text, 
                  AtReceiver text, ToReceiver text, FrequencyMode text,
                  SurveyTechnique text, Solution text, EpochInterval text, 
                  Class text, LevelDistance text, InstrumentModel text,
                  Derivation text, NON_GS text);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS DIR_TARGETS (
                  ID integer PRIMARY KEY, OBSERVATIONS_ID integer, 
                  TARGETS short text, VALUE text, TARGETS_SDEV text);''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS GNSS_SESSIONS (
                  ID integer PRIMARY KEY, NETWORK integer, 
                  SESSION integer, STATION_NAME text);''')
    conn.commit()
 
    
def DynaML2db(network):    
    conn = sqlite3.connect(network+'.db')
    cursor = conn.cursor()
    cursor.execute('DROP TABLE IF EXISTS STATIONS')
    cursor.execute('DROP TABLE IF EXISTS OBSERVATIONS')
    cursor.execute('DROP TABLE IF EXISTS DIR_TARGETS')
    cursor.execute('DROP TABLE IF EXISTS GNSS_SESSIONS')
    conn.commit()

    create_DynaML_db(conn)
    
    f = lg_xml(network +'.stn.xml')
    with open(f, 'r') as x:
        xml = xmltodict.parse(x.read())

    for s in xml['DnaXmlFormat']['DnaStation']:
        stn2db(s,conn)
        
    f = lg_xml(network+'.msr.xml')
    with open(f, 'r') as x:
        xml = xmltodict.parse(x.read())
    
    if type(xml['DnaXmlFormat']['DnaMeasurement'])==list:
        for m in xml['DnaXmlFormat']['DnaMeasurement']:
            msr2db(m,conn)
    else:
        msr2db(xml['DnaXmlFormat']['DnaMeasurement'],conn)
        
    conn.close()
    os.remove(f)
    
def db2Trivials(f):
    adjustment_name=f[:-3]
    conn = sqlite3.connect(f)
    cursor = conn.cursor()
    ####################################
    ### Search for Trivial Baselines ###
    ####################################
    del_kml=True
    tvl_kml = open(adjustment_name+'_tvl.kml', 'w')
    tvl_kml.write(kmlHeader(adjustment_name))
    
    #Break into individual networks for each session (event)
    sqlstring="SELECT DISTINCT OBSERVATIONS.StartDateTime AS EventTime \
        FROM OBSERVATIONS WHERE ((OBSERVATIONS.TYPE)='G') \
        ORDER BY EventTime;"
    Events=cursor.execute(sqlstring).fetchall()
    for e in Events:
        network_num=0; prevNtCnt=-1
        while prevNtCnt!=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall():
            prevNtCnt=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall()
            network_num=network_num+1
            # Seed the loop with a new network number
            sqlstring="INSERT INTO GNSS_SESSIONS (NETWORK, SESSION, STATION_NAME) \
            SELECT "+ str(network_num) + " AS NETWORK, "+ str(e[0]) + " AS SESSION, OBSERVATIONS.FIRST AS STATION_NAME \
            FROM OBSERVATIONS LEFT JOIN (SELECT GNSS_SESSIONS.NETWORK, GNSS_SESSIONS.STATION_NAME \
            FROM GNSS_SESSIONS \
            WHERE GNSS_SESSIONS.SESSION="+ str(e[0]) + ") AS M \
            ON OBSERVATIONS.FIRST = M.STATION_NAME \
            WHERE (((OBSERVATIONS.StartDateTime)<="+ str(e[0]) + ") AND ((OBSERVATIONS.FinishDateTime)>="+ str(e[0]) + ") AND ((OBSERVATIONS.TYPE)='G') AND ((M.STATION_NAME) Is Null)) \
            LIMIT 1"
            cursor.execute(sqlstring)
            conn.commit()
            prevCnt=-1
            while prevCnt!=cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall():
                prevCnt = cursor.execute("SELECT Count(GNSS_SESSIONS.STATION_NAME) AS CountOfSTATION_NAME FROM GNSS_SESSIONS;").fetchall()
                sqlstring="INSERT INTO GNSS_SESSIONS (NETWORK, SESSION, STATION_NAME) \
                SELECT C.NETWORK, C.SESSION, C.STATION_NAME \
                FROM (SELECT C_FIRST.NETWORK, C_FIRST.SESSION, C_FIRST.STATION_NAME \
                FROM( \
                SELECT MARKS_IN_SESS.NETWORK, MARKS_IN_SESS.SESSION, BASE_IN_SESS.FIRST AS STATION_NAME \
                FROM (SELECT GNSS_SESSIONS.NETWORK, GNSS_SESSIONS.SESSION, GNSS_SESSIONS.STATION_NAME \
                FROM GNSS_SESSIONS \
                WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)='"+ str(e[0]) + "')) \
                )  AS MARKS_IN_SESS  \
                INNER JOIN (SELECT OBSERVATIONS.FIRST, OBSERVATIONS.SECOND \
                FROM OBSERVATIONS \
                WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<="+ str(e[0]) + ") AND ((OBSERVATIONS.FinishDateTime)>="+ str(e[0]) + ")) \
                )  AS BASE_IN_SESS ON MARKS_IN_SESS.STATION_NAME = BASE_IN_SESS.SECOND) AS C_FIRST \
                UNION  \
                SELECT C_SECOND.NETWORK, C_SECOND.SESSION, C_SECOND.STATION_NAME \
                FROM( \
                SELECT MARKS_IN_SESS.NETWORK, MARKS_IN_SESS.SESSION, BASE_IN_SESS.SECOND AS STATION_NAME \
                FROM (SELECT GNSS_SESSIONS.NETWORK, GNSS_SESSIONS.SESSION, GNSS_SESSIONS.STATION_NAME \
                FROM GNSS_SESSIONS \
                WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)="+str(e[0])  + ")) \
                )  AS MARKS_IN_SESS INNER JOIN (SELECT OBSERVATIONS.FIRST, OBSERVATIONS.SECOND \
                FROM OBSERVATIONS \
                WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<="+ str(e[0]) + ") AND ((OBSERVATIONS.FinishDateTime)>="+ str(e[0]) + ")) \
                )  AS BASE_IN_SESS ON MARKS_IN_SESS.STATION_NAME = BASE_IN_SESS.FIRST) AS C_SECOND)  AS C LEFT JOIN (SELECT GNSS_SESSIONS.STATION_NAME \
                FROM GNSS_SESSIONS \
                WHERE (((GNSS_SESSIONS.NETWORK)="+ str(network_num) + ") AND ((GNSS_SESSIONS.SESSION)="+ str(e[0]) + ")) \
                )  AS MARKS_IN_SESS ON C.STATION_NAME = MARKS_IN_SESS.STATION_NAME \
                WHERE (((MARKS_IN_SESS.STATION_NAME) Is Null))"
                cursor.execute(sqlstring)
                conn.commit()
    
    ## Count the number of baselines and stations
    ## used in a session and test if #_baselines > #stations - 1
    sqlstring="SELECT DISTINCT NETWORK, SESSION \
        FROM GNSS_SESSIONS \
        ORDER BY SESSION, NETWORK"
    session=cursor.execute(sqlstring).fetchall()
    for s in session:
        baseCntSql = "SELECT COUNT(OBSERVATIONS.ID) AS BASE_CNT \
            FROM (GNSS_SESSIONS INNER JOIN OBSERVATIONS ON GNSS_SESSIONS.STATION_NAME = OBSERVATIONS.FIRST) \
            INNER JOIN GNSS_SESSIONS AS GNSS_SESSIONS_1 ON OBSERVATIONS.SECOND = GNSS_SESSIONS_1.STATION_NAME \
            WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<="+ str(e[0]) + ") AND ((OBSERVATIONS.FinishDateTime)>"+ str(e[0]) + ") AND ((GNSS_SESSIONS.SESSION)="+ str(s[1]) + ") AND ((GNSS_SESSIONS_1.SESSION)="+ str(s[1]) + ") AND ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS_1.NETWORK)="+ str(s[0]) + "))"
        stationCntSql = "SELECT COUNT(ID) AS STATION_CNT \
            FROM GNSS_SESSIONS \
            WHERE ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS.SESSION)="+ str(s[1]) + ")"
        stationCnt = cursor.execute(stationCntSql).fetchall()
        baseCnt = cursor.execute(baseCntSql).fetchall()
        if stationCnt[0][0]-1<baseCnt[0][0]:
            del_kml=False
            kmlsql="SELECT OBSERVATIONS.ID, OBSERVATIONS.FIRST, OBSERVATIONS.SECOND , STATIONS.LATITUDE, STATIONS.LONGITUDE, STATIONS_1.LATITUDE, STATIONS_1.LONGITUDE, OBSERVATIONS.StartDateTime, OBSERVATIONS.FinishDateTime \
                    FROM (((OBSERVATIONS INNER JOIN GNSS_SESSIONS ON OBSERVATIONS.FIRST = GNSS_SESSIONS.STATION_NAME) INNER JOIN GNSS_SESSIONS AS GNSS_SESSIONS_1 ON OBSERVATIONS.SECOND = GNSS_SESSIONS_1.STATION_NAME) LEFT JOIN STATIONS AS STATIONS_1 ON GNSS_SESSIONS_1.STATION_NAME = STATIONS_1.STATION_NAME) LEFT JOIN STATIONS ON GNSS_SESSIONS.STATION_NAME = STATIONS.STATION_NAME \
                    WHERE (((OBSERVATIONS.TYPE)='G') AND ((OBSERVATIONS.StartDateTime)<="+ str(e[0]) + ") AND ((OBSERVATIONS.FinishDateTime)>"+ str(e[0]) + ") AND ((GNSS_SESSIONS.SESSION)="+ str(s[1]) + ") AND ((GNSS_SESSIONS_1.SESSION)="+ str(s[1]) + ") AND ((GNSS_SESSIONS.NETWORK)="+ str(s[0]) + ") AND ((GNSS_SESSIONS_1.NETWORK)="+ str(s[0]) + "))"
            kml_Trivial=cursor.execute(kmlsql).fetchall()  
            for tv in kml_Trivial:
                tvl_kml.write(MkLine(tv,s))
    tvl_kml.write(kmlFooter())
    tvl_kml.close()
    if del_kml: os.remove(adjustment_name+'_tvl.kml')
    conn.close()

