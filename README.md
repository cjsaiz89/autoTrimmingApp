# autoTrimmingApp
Hurricane Seagliders auto trimmer app: reads the NetCDF files containing dive data and returns the recommended values to adjust the command file.

- Execute the autoTrimApp.exe (Windows only). It needs [Matlab runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html#:~:text=The%20MATLAB%20Runtime%20is%20a%20standalone%20set%20of,numerical%20applications%20or%20software%20components%20quickly%20and%20securely).
- Select the directory named as "sgGGG" where all the NetCDF files are located. *When you select to process diveN and J dives back, all the .nc files from diveN to diveN-J must be in the directory otherwise it will not process.*
- Use the switch to toggle the type of average:
  - Simple --> for J dives --> (diveN + diveN-1 + diveN-2 + ... + diveN-J ) / J
  - Weighed --> for J dives --> (diveN * J + diveN-1 * (J-1) + diveN-2 * (J-2) + ... + diveN-J * 1 ) / (J + (J-1) + (J-2) + ... + 1) 

The *sgGGG_cmdfile* is created in the same directory as the program, and it will have the implied parameters, like:
>$MAX_BUOY,120  |
>$SM_CC,380   | 
>$C_VBD,2835  |
>$C_PITCH,2654    | 
>$PITCH_GAIN,27.1   |
>$C_ROLL_DIVE,2154  |
>$C_ROLL_CLIMB,2080 

MAXBUOY is not implied but can be increased/decreased by checking the last field 'DISTKM' (distance between GPS fixes in km), for the last dives and change accordingly.

![autoTrimmingApp](https://user-images.githubusercontent.com/89260258/187809712-9e9be9c2-cfa3-4699-91dc-3a2e92d68117.png)
