def download(today = True):
    """
    1. Description
    : Download IMSNG imaging data using SCP and crontab.

    * Main commands ;
    > ssh -p 222 admin@***.***.***.** # PWD : *********
    > ls /volume1/1mobs/STX16803/
    > scp -p port id@ip:path_original upper_directory_for_destination

    On the terminal at qso server, type this command.
    > scp -P 222 -r admin@147.46.139.34:/volume1/1mobs/STX16803/20190516 /data3/SAO1m/obsdata/STX16803

    2. Usage
    
    if today is True, the code uses today's date.
    if today is False, You should put specific date like '20190606'

    3. History
    2019.06.05 : 

    """
    import glob
    import os, sys
    from datetime import datetime

    nas_data_path = '/volume1/1mobs/STX16803/'
    qso_data_path = '/data3/SAO1m/obsdata/STX16803/'
    qso_red_path  = '/data3/SAO1m/red/STX16803/'

    # Are there data yesterday?
    if today == True :
        now = datetime.now()
        curdate = str(now.year)+str(format(now.month,'02'))+str(format((now.day-1), '02'))
    elif today == False :
        curdate = input('Specify date : ')
    try:
        print("Copy today's data on obsdata path...")
        sshpass = os.system('sshpass -p*********** scp -o StrictHostKeyChecking=no -P 222 -r admin@***.***.***.**:'+ nas_data_path+curdate+' '+qso_data_path)
        if sshpass == 0 :
            print('Data is downloaded. Copy to working directory.')
            #os.system('cp -r '+qso_data_path+curdate+' '+qso_red_path)
            #print('Data copy is done. See you tomorrow!')
        elif sshpass == 256 :
            print('No data in '+curdate+'.')
    except :
        print('No Recent data.')
    exit()

download(today=True)
#download(today=False)
