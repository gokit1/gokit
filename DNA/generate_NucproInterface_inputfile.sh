#! /usr.bin/bash
echo ";"  > nucpro_interface.input
echo ";Input paramters for nucpro_interface."	>> nucpro_interface.input
echo ";usage $ python nucproSBM.py --interface" >> nucpro_interface.input
echo ";teh terms in square brackets are not readed by the program. Donot remove any row as the input is based on the line index" >> nucpro_interface.input
echo "; Add C10/C6 attreaction term between DNA/RNA and Protein" >> nucpro_interface.input
echo ";Example:" >> nucpro_interface.input
echo ";[ LJ attraction (False) ] True" >> nucpro_interface.input
echo "[ LJ attraction (False) ] " >> nucpro_interface.input
echo ";Include LJ term for aromatic stacking" >> nucpro_interface.input
echo "[ Aromacti Stacking (True)]" >> nucpro_interface.input
echo "; Radius of beads N bases (BA,BT,BG,BC,BU), Sugar (Deoxy=DS,Oxy=S), Phosphate (XP), CA, CB" >> nucpro_interface.input
echo ";Example:" >> nucpro_interface.input
echo ";[ BA radius (Same as bulk) ] 2.0 " >> nucpro_interface.input
echo "[ BA radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ BG radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ BT radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ BC radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ BU radius (Same as bulk) ]" >> nucpro_interface.input
echo "[  S radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ XP radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ CA radius (Same as bulk) ]" >> nucpro_interface.input
echo "[ CB radius (Same as bulk) ]" >> nucpro_interface.input
echo ";Use Native (True) or non-specific electostatic interaction. Only applicable for native DNA" >> nucpro_interface.input
echo "[Native Charge (False) ]" >> nucpro_interface.input
