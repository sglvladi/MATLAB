from openpyxl import load_workbook
import ecef
import pickle
import json

wb = load_workbook("C:\Users\sglvladi\Documents\WebSync\University of Liverpool\PhD\resources\geolife.xlsx")
gps_data = wb['Sheet1']
in_file = open("/home/sglvladi/Mini_Challenge_2/cart_table.txt","r");
i=1;
for line in in_file:
	tokens = line.split(",")
	taxi_id = tokens[0]
	timestamp = tokens[1]
	lat = tokens[2]
	long = tokens[3]
	gps_data["A%s"%i]=taxi_id
	gps_data["B%s"%i]=timestamp
	gps_data["C%s"%i]=lat
	gps_data["D%s"%i]=long
	
