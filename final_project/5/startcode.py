# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
line=[]

import csv
with open('hep-a-mi.csv', 'rb') as csvfile:
     reader = csv.reader(csvfile)
     for row in reader:
         
         line.append(row)
del line[0]
import numpy as np
num=np.zeros((1,269))
num=num[0]
final=[]
for i in range(269):
    num[i]=int(line[i*4][2])+int(line[i*4+1][2])+int(line[i*4+2][2])+int(line[i*4+3][2])
    final.append([i+1,num[i]])

with open("hep-a-mi2.csv","wb") as output:
    writer= csv.writer(output)
    writer.writerows(final)
    

         
         
    

