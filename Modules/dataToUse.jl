### Contain the loading of all the data used in the thesis.

using  CSV, DataFrames, DelimitedFiles ;

#### Import data from CVS file
nb_days2load = 365 # the total number of day 
Pv_max = 4;
firstRow = 2; 
lastRow = 48*nb_days2load; # Number of line to read in the file


data2010_11 = DataFrame(CSV.File("/home/supelec/My_Jupiter/GitJupLab/SolarHomeData/2010-2011/2010_2011_Cstm61_90.csv";
    skipto=firstRow, limit=lastRow));
data2012_13 = DataFrame(CSV.File("/home/supelec/My_Jupiter/GitJupLab/SolarHomeData/2012-2013/2012_2013_Cstm61_90.csv";
    skipto=firstRow, limit=lastRow));


# 2010 2011----------------------------------------------

Pv_star = Pv_max*data2010_11.GG61/data2010_11.Generator_Capacity61[1]; 
Pl_star = data2010_11.GC61;

dt=.5;

# reshape in order to have each row representing a day
PvData2010_11 = reshape(Pv_star,(48,nb_days2load))'
LoadData2010_11 = reshape(Pl_star,(48,nb_days2load))';


month2010_11 = [
            "Jul 2010"
            "Aug 2010"
            "Sep 2010"
            "Oct 2010"
            "Nov 2010"
            "Dec 2010"
            "Jan 2011"
            "Fev 2011"
            "Mar 2011"
            "Apr 2011"
            "May 2011"
            "Jun 2011" ];


# reshape data so that each row represent a day 
PvData = reshape(Pv_star,(48,nb_days2load))'
LoadData = reshape(Pl_star,(48,nb_days2load))';


Jul_10_Ld = LoadData2010_11[1:31,:]
Aug_10_Ld = LoadData2010_11[32:62,:]
Sep_10_Ld = LoadData2010_11[63:92,:]
Oct_10_Ld = LoadData2010_11[93:123,:]
Nov_10_Ld = LoadData2010_11[124:153,:]
Dec_10_Ld = LoadData2010_11[154:184,:]
Jan_11_Ld = LoadData2010_11[185:215,:]
Fev_11_Ld = LoadData2010_11[216:243,:]
Mar_11_Ld = LoadData2010_11[244:274,:] 
Apr_11_Ld = LoadData2010_11[275:304,:] 
May_11_Ld = LoadData2010_11[305:335,:] 
Jun_11_Ld = LoadData2010_11[336:365,:];

Jul_10_Pv = PvData2010_11[1:31,:]
Aug_10_Pv = PvData2010_11[32:62,:]
Sep_10_Pv = PvData2010_11[63:92,:]
Oct_10_Pv = PvData2010_11[93:123,:]
Nov_10_Pv = PvData2010_11[124:153,:]
Dec_10_Pv = PvData2010_11[154:184,:]
Jan_11_Pv = PvData2010_11[185:215,:]
Fev_11_Pv = PvData2010_11[216:243,:]
Mar_11_Pv = PvData2010_11[244:274,:] 
Apr_11_Pv = PvData2010_11[275:304,:] 
May_11_Pv = PvData2010_11[305:335,:] 
Jun_11_Pv = PvData2010_11[336:365,:];

PlStarYear2010_11=[ 
                [Jul_10_Ld]
                [Aug_10_Ld]
                [Sep_10_Ld]
                [Oct_10_Ld]
                [Nov_10_Ld]
                [Dec_10_Ld]
                [Jan_11_Ld]
                [Fev_11_Ld]
                [Mar_11_Ld]
                [Apr_11_Ld]
                [May_11_Ld]
                [Jun_11_Ld]   ];

PvStarYear2010_11=[ 
                [Jul_10_Pv]
                [Aug_10_Pv]
                [Sep_10_Pv]
                [Oct_10_Pv]
                [Nov_10_Pv]
                [Dec_10_Pv]
                [Jan_11_Pv]
                [Fev_11_Pv]
                [Mar_11_Pv]
                [Apr_11_Pv]
                [May_11_Pv]
                [Jun_11_Pv]   ];



# 2012 2013----------------------------------------------

Pv_star = Pv_max*data2012_13.GG61/data2012_13.Generator_Capacity61[1]; 
Pl_star = data2012_13.GC61;



PvData2012_13 = reshape(Pv_star,(48,nb_days2load))'
LoadData2012_13 = reshape(Pl_star,(48,nb_days2load))';


month2012_13 = [
            "Jul 2012"
            "Aug 2012"
            "Sep 2012"
            "Oct 2012"
            "Nov 2012"
            "Dec 2012"
            "Jan 2013"
            "Fev 2013"
            "Mar 2013"
            "Apr 2013"
            "May 2013"
            "Jun 2013" ];


Jul_12_Ld = LoadData2012_13[1:31,:]
Aug_12_Ld = LoadData2012_13[32:62,:]
Sep_12_Ld = LoadData2012_13[63:92,:]
Oct_12_Ld = LoadData2012_13[93:123,:]
Nov_12_Ld = LoadData2012_13[124:153,:]
Dec_12_Ld = LoadData2012_13[154:184,:]
Jan_13_Ld = LoadData2012_13[185:215,:]
Fev_13_Ld = LoadData2012_13[216:243,:]
Mar_13_Ld = LoadData2012_13[244:274,:] 
Apr_13_Ld = LoadData2012_13[275:304,:] 
May_13_Ld = LoadData2012_13[305:335,:] 
Jun_13_Ld = LoadData2012_13[336:365,:];

Jul_12_Pv = PvData2012_13[1:31,:]
Aug_12_Pv = PvData2012_13[32:62,:]
Sep_12_Pv = PvData2012_13[63:92,:]
Oct_12_Pv = PvData2012_13[93:123,:]
Nov_12_Pv = PvData2012_13[124:153,:]
Dec_12_Pv = PvData2012_13[154:184,:]
Jan_13_Pv = PvData2012_13[185:215,:]
Fev_13_Pv = PvData2012_13[216:243,:]
Mar_13_Pv = PvData2012_13[244:274,:] 
Apr_13_Pv = PvData2012_13[275:304,:] 
May_13_Pv = PvData2012_13[305:335,:] 
Jun_13_Pv = PvData2012_13[336:365,:];

PlStarYear2012_13=[ 
                [Jul_12_Ld]
                [Aug_12_Ld]
                [Sep_12_Ld]
                [Oct_12_Ld]
                [Nov_12_Ld]
                [Dec_12_Ld]
                [Jan_13_Ld]
                [Fev_13_Ld]
                [Mar_13_Ld]
                [Apr_13_Ld]
                [May_13_Ld]
                [Jun_13_Ld]   ];

PvStarYear2012_13=[ 
                [Jul_12_Pv]
                [Aug_12_Pv]
                [Sep_12_Pv]
                [Oct_12_Pv]
                [Nov_12_Pv]
                [Dec_12_Pv]
                [Jan_13_Pv]
                [Fev_13_Pv]
                [Mar_13_Pv]
                [Apr_13_Pv]
                [May_13_Pv]
                [Jun_13_Pv]   ];
