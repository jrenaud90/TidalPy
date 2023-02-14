
import subprocess
import time
changes = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
changes = [0.1, 0.3, 0.6, 0.9]

with open('test.py', 'r') as pyfile:
    lines = pyfile.readlines()

def run():
    
    for change in changes:
        
        print('Preping alpha', change)
        lines[19] = "venus_solar_day_freq = days2rads(" + f'{change}' + ")\n"
        # lines[67] = "complex_shear = andrade(venus_solar_day_freq, shear_array**(-1), viscosity_array, " + f'{change}' + ", 1.0)**(-1)\n"
        with open('test.py', 'w') as pyfile:
            pyfile.writelines(lines)
            
        print('Calling Script')
        subprocess.Popen('python -W ignore test.py', shell=False,
             stdin=None, stdout=None, stderr=None, close_fds=True)
             
        print('sleeping')
        time.sleep(2)
            
run()
            
            
            