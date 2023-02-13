# ------------------------------------------------------------------------------
# -------------------------------------------------------------- IMPORT PACKAGES
import os
import shutil
import sys
import math

# ------------------------------------------------------------------------------
# --------------------------------------------------------------------- ADD PATH
sys.path.insert(0, os.getcwd() + '/htpp_run')

# ------------------------------------------------------------------------------
# ----------------------------------------------------------------- IMPORT FILES
from htpp_run_simulation import htpp_run_simulation

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------- MAIN
# ------------------------------------------------------------------------------
def htpp_run():
    # ------------------------------------------------------ FIND SETTINGS FILES
    path = os.getcwd()[:-4]
    files = []
    for file in os.listdir(path):
        if file.endswith('.htpp'):
            files.append(file)

    for file in files:
        # ------------------------------------------ OPEN AND READ SETTINGS FILE
        data_settings = []
        try:
            path = os.getcwd()[:-4]
            with open(path + file, 'r') as f:
                for line in f:
                    if len(line.split()) > 0:
                        if line.split()[0][0:2][0:2] == '\\*':
                            continue
                        else:
                            data_settings.append(line.split(", "))
        except Exception as e:
            print('Error: ' + file)
            return

        i = 1
        for line in data_settings:
            if '*ODB' in line[0].upper():
                odbFile = line[-1][6:-2]
                name = odbFile
            elif '*RUN' in line[0].upper():
                data_begin = i
            elif '*ENDRUN' in line[0].upper():
                data_end = i
            i += 1

        # ------------------------------------------------- EXTRACT DATA OPTIONS
        data_options = data_settings[data_begin:data_end]
        all_data = []
        for line in data_options:
            if '**JOB' in line[0].upper():
                all_data.append(line)


        # ------------------------------------------------ CALL DATA SUBROUTINES
        if len(all_data) > 0:
            tlen = 20.0
            st = int(round(tlen / len(all_data)))
            stri = "> HTPP RUN  "
            sys.stdout.write(stri+"[%s]    %d%s  (%s)"%(" "*int(tlen), 0,'%',name))
            sys.stdout.flush()
            n = 1
            for option in all_data:
                opt = option[0]
                if opt == '**JOB':
                    htpp_run_simulation(option,name)

                per = 100.0 * float(st * n) / float(tlen)
                bars = st * n
                form = "[%s%s] ~ %d%s  (%s)"
                if per >= 100.0:
                    per = 100.0
                    bars = int(tlen)
                    form = "[%s%s]  %d%s  (%s)"
                sys.stdout.write("\r"+stri+form %("-"*bars," "*(int(tlen)-bars),per,'%',name))
                sys.stdout.flush()
                n += 1
            sys.stdout.write("\r"+stri+"[%s]  %d%s  (%s)\n"%("-"*int(tlen),100,'%',name))

        # --------------------------------------------------- SAVE ODB AND CLOSE

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    htpp_run()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
