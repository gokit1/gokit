sed -i s/"BA"/"NA"/g $1
sed -i s/"BG"/"NG"/g $1
sed -i s/"BC"/"NC"/g $1
sed -i s/"BT"/"NT"/g $1
sed -i s/"BU"/"NU"/g $1
sed -i s/"XP"/" P"/g $1
sed -i s/"RS"/"CR"/g $1
sed -i s/"DS"/"CD"/g $1
######################################
#changing back to CBA (it as changed to CNA in above sed commands)
sed -i s/"CNA"/"CBA"/g $1
sed -i s/"CNC"/"CBC"/g $1
sed -i s/"CNT"/"CBT"/g $1
