
import time
from scipy.optimize import curve_fit
import numpy as np
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pltlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import Tkinter as Tk
from matplotlib.figure import Figure
import Tkinter
from Tkinter import *
import scipy as sc

class variables:
    def __init__(self):
        self.side=0.9
        self.airDensity=1.225
        self.noOfSides=[5,1]
        self.cAir=1.01 #J/gram-celcius
        self.diameter=12.0/100.0
        self.transCoeff=[0.1/.02,4]
        self.TRoom=32
        self.TIn=32
        self.inletSpeed=10
        self.stefanBoltzmann=5.670373*(10**(-8))
        self.emissivity=0.4
        self.tungstenArea=3*(10**(-5))
        self.wattageRatio=400/60
        self.tungstenTemp=3000
        self.sample=1
        self.bulbWattage=100
        self.percentOpen=100
        self.wallRoomCoeff=10
        self.measuredTemp=self.TRoom
        self.numberOfHumans=0
        self.numberOfBulbs=2
        self.recalculate()

    def recalculate(self):
        self.boxAirMass=self.airDensity*(self.side**3)
        self.area=self.side**2
        self.specMassAir=self.boxAirMass*self.cAir*1000
        self.massPerSecond=self.sample*(math.pi*self.diameter*self.diameter/4)*self.inletSpeed*self.airDensity*(self.percentOpen/100)
        self.plywoodCMass=(self.area)*(2/100.0)*545*1215*5.2*0.2
        self.wattage=self.numberOfHumans*120+self.numberOfBulbs*self.bulbWattage


##########   Setters   ###########        
    def setSide(self,side):
        self.side=side
        self.recalculate()
    def setAirDensity(self,airDensity):
        self.airDensity=airDensity
        self.recalculate()
    def setNoOfSides(self,noOfSides):
        self.noOfSides=noOfSides
        self.recalculate()
    def setCAir(self,cAir):
        self.cAir=cAir
        self.recalculate()
    def setDiameter(self,diameter):
        self.diameter=diameter
        self.recalculate()
    def setTransCoeff(self,transCoeff):
        self.transCoeff=transCoeff
        self.recalculate()
    def setTRoom(self,TRoom):
        self.TRoom=TRoom
        self.recalculate()
    def setTIn(self,TIn):
        self.TIn=TIn
        self.recalculate()
    def setInletSpeed(self,inletSpeed):
        self.inletSpeed=inletSpeed
        self.recalculate()    
    def setEmissivity(self,emissivity):
        self.emissivity=emissivity
    def setPercentOpen(self,percentOpen):
        self.percentOpen=percentOpen
        self.recalculate()
    def setTungstenArea(self,tungstenArea):
        self.tungstenArea=tungstenArea
        
    def setWattageRatio(self,wattageRatio):
        self.wattageRatio=wattageRatio
       
    def setSample(self,sample):
        self.sample=sample
        self.recalculate() 
    def setTungstenTemp(self,tungstenTemp):
        self.tungstenTemp=tungstenTemp
    def setBulbWattage(self,bulbWattage):
        self.bulbWattage=bulbWattage
    def setNumberOfBulbs(self,n):
        self.numberOfBulbs=n
    def setNumberOfHumans(self,n):
        self.numberOfHumans=n
    def getWattage(self):
        return self.wattage
##################################
class model:
    def __init__(self):
        self.var=variables()
        self.allTemps=[]
        self.i=0
        
        self.setPoint=0
        self.wallTemp=self.var.TRoom
    
    def getConduction(self):
        conduction=0
        for i in range(len(self.var.noOfSides)):
            conduction+=(self.var.transCoeff[i]*self.var.noOfSides[i])
     
        conduction=conduction*(self.wallTemp-self.var.TIn)*self.var.area
        return conduction
     

    def getRadiation(self):
        return (self.var.emissivity*self.var.stefanBoltzmann*self.var.tungstenArea*(self.var.tungstenTemp**4))

    def getRadiation2(self):
        return (self.var.getWattage()*.95/6)

    def getNewMixedHeat(self):
        massLeft=self.var.boxAirMass-self.var.massPerSecond
        t1=(massLeft*self.var.TRoom)
        t2=(self.var.massPerSecond*self.var.TIn)
        t=(t1+t2)*(self.var.cAir*1000)
        #print t1+t2
        return t
    def wallTempFunc(self):
        self.wallTemp+=(.95*self.var.getWattage()*self.var.sample-self.wallsAirTransfer()-self.getConduction()-self.getRadiation2())/self.var.plywoodCMass

    
    def wallsAirTransfer(self):
        return self.var.wallRoomCoeff*self.var.area*5.2*(self.wallTemp-self.var.TRoom)
      
    def getTRoom(self):
        return self.var.TRoom
    def getTIn(self):
        return self.var.TIn

    def setTIn(self,TIn):
        self.var.setTIn(TIn)
    def setTRoom(self,TRoom):
        self.var.setTRoom(TRoom)
    def getBulbWattage(self):
        return self.var.bulbWattage
    def setbulbWattage(self, bulbWattage):
        self.var.bulbWattage=bulbWattage
    def getBoxAirMass(self):
        return self.var.boxAirMass
    def setPercentOpen(self,percentOpen):
        self.var.setPercentOpen(percentOpen)
    def setI(self,i):
        self.i=i
    def setSetPoint(self,setPoint):
        self.setPoint=setPoint
    def getSetPoint(self):
        return self.setPoint



########## MAIN PART OF MODEL####################
 
    def runModel(self):

        newTRoom=self.getTRoom()
        
        self.wallTempFunc()
       
        """
        self.pid.setTRoom(self.var.measuredTemp) #interface with the PID sensor is sending value to PID
        self.pid.setSetPoint(self.setPoint)# PID is being told the setPoint from the GUI
        if self.tuningStatus='Autotune':
        
        self.setPercentOpen(50)#self.pid.getPid()) #PID is retruning the control variable.
        """
        #do emissivity of walls -> which takes the radiation and heats up and then the walls supply the air with heat. Conduction is with the walls temps
        deltaNewTRoom=(self.getNewMixedHeat()+self.wallsAirTransfer())/(self.getBoxAirMass()*self.var.cAir*1000)
        newTRoom=deltaNewTRoom
        self.setTRoom(newTRoom)
       
    def getAllTemps(self):
        self.runModel()
        return self.allTemps




####   PID   ######

class PID:
    def __init__(self,model):
        
        self.model=model
        self.heur=heuristic(self.model)
        self.TRoom=self.model.getTRoom()
        self.setPoint=0
        self.Ku=3.14*(100)/(2*.47)
        self.Pu=8
        self.Kp=12.9595829#.75*self.Ku*20
        self.Ki=0.044444#0.5/(0.625*self.Pu)
        self.Kd=5.625*100#self.Pu/10
        self.prevError=0
        self.sumError=2000
        self.prop_min=0
        self.prop_max=100
        self.prop=0
        self.measured=self.TRoom
        self.alpha=0.95
        self.prevProp=self.prop
        self.it=0
        self.currentIt=0
        self.riseStatus=True
        self.trend=0
        self.TRoomArray=[]
        
    def getProp(self):
        prop=self.Kp*(self.TRoom-self.setPoint)
        
        return prop  
    def getInt(self):
        if self.prop!=self.prop_max and self.prop!=self.prop_min:
            self.sumError+=(self.TRoom-self.setPoint)
        print self.sumError
        return self.Ki*self.sumError
    def getDiff(self):
        error=(self.TRoom-self.setPoint)
        diff=self.Kd*(error-self.prevError)
        self.prevError=error
        return diff
 
                  #  main #  method  #  for  #  PID #
    def Pid(self):
        self.measured=np.random.normal(loc=self.model.getTRoom(),scale=0.5)
        self.filtered= self.alpha*self.TRoom+(1-self.alpha)*self.measured
        self.TRoom=self.filtered
        self.TRoomArray.append(self.TRoom)
        self.prop= self.getProp()+self.getDiff()+self.getInt()
        if (self.prop<self.prop_min):
            self.prop= self.prop_min
        if (self.prop>self.prop_max):    
            self.prop= self.prop_max
        #print self.prop
        #self.prop=0.95*self.prevProp+0.05*self.prop
        self.prevProp=self.prop
        self.model.setPercentOpen(self.prop)
        self.heur.setPercentOpen(self.prop)
        self.heur.calcCostFunction()
        #self.Kp+=self.heur.delKp()
        #print self.heur.delKp()
        #print self.Kp
        #print 'KP', self.Kp, 'Ki ', self.Ki, 'Kd ', self.Kd
        #print '_____________'
        self.it+=1
        
        self.trend+=(self.setPoint-self.TRoom)


        if abs(self.setPoint-np.mean([temp for temp in self.TRoomArray[-10:]]))<0.1 and self.riseStatus==True and self.it-self.currentIt>20:
            
            self.Kp+=self.heur.timeToReach((self.it-self.currentIt),(self.trend/abs(self.trend)))
            self.riseStatus=False
            self.heur.reset()
            #b=curveFit(self.TRoomArray[-self.currentIt+100:-self.currentIt+150]).curvefit()
            #print 'Time Contant=  ', 1.0/b[1]
        self.model.runModel()

                    #       #         #        #      #
  
    def setKu(self,Ku):
        self.Ku=Ku
        self.autoSet()
    def setPu(self,Pu):
        self.Pu=Pu
        self.autoSet()
    def autoSet(self):
        self.Kp=0.6*self.Ku
        self.Ki=1/(0.5*self.Pu)
        self.Kd=0.125*self.Pu
        #print 'KP', self.Kp, 'Ki ', self.Ki, 'Kd ', self.Kd
    def setKp(self,Kp):
        self.Kp=Kp
    def setKi(self,Ki):
        self.Ki=Ki
    def setKd(self,Kd):
        self.Kd=Kd
    def getKp(self):
        return self.Kp
    def setKi(self):
        return self.Ki
    def setKd(self):
        return self.Kd
    def setTRoom(self,TRoom):
        self.TRoom=TRoom
    def setSetPoint(self,setPoint):
        self.setPoint=setPoint
    def getTRoom(self):
        return self.TRoom
    
#########################################
class autoTuneRelay:
    def __init__(self,model):
        self.model=model
        self.open=0
        self.status=True
        self.setPoint=0
        self.alpha=0.6
        self.TRoom=self.model.getTRoom()
        self.openTime=[]
        self.closedTime=[]
        self.previ=self.model.i
        self.peakStatus=False
        self.peaknumber=0
        self.cycles=5
        self.closedPercent=0.0
        self.openPercent=100.0
        self.allTemps=[]
    def autoTune(self):
        measured=np.random.normal(loc=self.model.getTRoom(),scale=.5)
        filtered= self.alpha*self.TRoom+(1-self.alpha)*measured
        self.TRoom=filtered#self.model.getTRoom()
        if(self.TRoom<self.setPoint*0.98 and self.status):
            self.status=False
            self.open=self.closedPercent
            self.closedTime.append(self.model.i-self.previ)
            self.previ=self.model.i
            self.model.setPercentOpen(self.open)
            self.peakStatus=True
            self.peaknumber+=1
        elif (self.TRoom>=self.setPoint*1.02 and not self.status):
            self.status=True
            self.open=self.openPercent
            self.openTime.append(self.model.i-self.previ)
            self.model.setPercentOpen(self.open)
            self.peakStatus=True
            self.previ=self.model.i
        self.model.runModel()
        if self.peaknumber>=2:
            self.allTemps.append(self.TRoom)
        #print self.status
        if self.peaknumber>=self.cycles  and self.peakStatus==True:
            openTimeAvg=(sum([time for time in self.openTime[1:]])/(len(self.openTime)-1))
            closedTimeAvg=(sum([time for time in self.closedTime[1:]])/(len(self.closedTime)-1))
            a=max(self.allTemps)-min(self.allTemps)
            e=.02*self.setPoint
            pi=3.14
            delta=self.openPercent-self.closedPercent
            self.ku=4*delta/(pi*(((a**2)-(e**2))**.5))
            self.pu=openTimeAvg+closedTimeAvg
            self.peakStatus=False
            
         
   
    def setSetPoint(self,setPoint):
        self.setPoint=setPoint
    def getTRoom(self):
        return self.TRoom
        
    
 
################ Auto Tune Heuristic#########################    
class heuristic:
    def __init__(self,model):
        self.model=model
        self.costFunction=0
        self.air=0
        self.percentOpen=0
        self.i=1
        self.prevAvgAir=0
        #self.timeToReach=0
    def calcCostFunction(self):
        if self.percentOpen!= 0.0:
            self.air+=self.percentOpen
            self.i+=1
        
    def delKp(self):
        if self.percentOpen!=0:
            return (-self.avgAirDiff/100)
        else: return 0.0
    def timeToReach(self, riseTime, trend):
        if trend==-1:
            self.avgAir=self.air/self.i
            self.avgAirDiff=self.prevAvgAir-self.avgAir
            self.prevAvgAir=self.avgAir
            print 'Difference in air volume: ',self.avgAirDiff 
            print riseTime
            return riseTime/50.0+self.avgAirDiff/5
        return 0
    def reset(self):
        self.i=1
        self.air=0
        self.avgAir=0
    def setPercentOpen(self,percent):
        self.percentOpen=percent
#############################################################
class curveFit:
    def __init__(self, y):
         
         self.x=np.array(range(len(y)))
         self.y=np.array(y)
         self.y=self.y-(min(self.y)*np.ones(len(self.y)))
    def func(self,x,a,b,c):
        return a * np.exp(-b * x)+c 
    
    def curvefit(self):
        self.popt,self.pcov=curve_fit(self.func,self.x,self.y)

        return self.popt
      
################  Graphical User Interface     ############

class App_Window(Tkinter.Tk):
    def __init__(self,parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent=parent
        self.initialize()
    def initialize(self):
        self.grid()
        self.it=0
        self.model=model()
        self.pid=PID(self.model)
        self.ATR=autoTuneRelay(self.model)
        self.model.setTRoom(32)
        self.model.setTIn(32)
        self.title("HVAC Model with AutoTune PID")
#--------------------------------------------------------
        
        box1=Frame(master=self, height=2, bd=1, padx=5, pady=5, relief=SUNKEN) 
        box2= Frame(master=self, height=10, bd=1, padx=5, pady=5, relief=SUNKEN)
        box3=Frame(master=self,height=10,width=10,bd=1, padx=5,pady=5,relief=SUNKEN)
        box1.grid(column=0,row=1,columnspan=3)
        box2.grid(column=0,row=3)
        box3.grid(column=0,row=4)
        self.secLabel=Label(box2, text="Time= "+str(self.it)+" secs")
        self.tempLabel=Label(box2, text="Room Temperature= %.2f" %self.model.getTRoom()+" degree C")
        self.dampLabel=Label(box2,text="Damper Position= %.2f" %self.model.var.percentOpen+"%")
        self.secLabel.pack(anchor='w')
        self.tempLabel.pack(anchor='w')
        self.dampLabel.pack(anchor='w')       
        self.config(padx=20, pady=20)

        button=Tkinter.Button(self, text="Start Process \n Submit SetPoint", command=self.OnStartClick).grid(column=0,row=2)

        button2=Tkinter.Button(self, text="Start Relay Tune",command=self.startRelayTune).grid(column=0,row=5)

        closeButton=Tkinter.Button(self, text="Close Window", command=self.close_window).grid(column=0,row=6)

        Label(box1,text="Enter Set Point").pack()
        #self.model.setbulbWattage(200)
        buttonBulbRise=Tkinter.Button(master=box3, text="Increase bulb", command=self.incrementBulb)
        buttonBulbRise.grid(column=0, row=0,padx=2,pady=2)
        buttonBulbFall=Tkinter.Button(master=box3, text="Decrease bulb", command=self.decrementBulb)
        buttonBulbFall.grid(column=1, row=0,padx=2,pady=2)
        self.bulbNumberLabel=Label(master=box3, text='No of Bulbs= %d' %self.model.var.numberOfBulbs)
        self.bulbNumberLabel.grid(column=2,row=0,padx=2,pady=2)
        buttonHumanRise=Tkinter.Button(master=box3, text="Increase human", command=self.incrementHuman)
        buttonHumanRise.grid(column=0, row=1,padx=2,pady=2)
        buttonHumanFall=Tkinter.Button(master=box3, text="Decrease human", command=self.decrementHuman)
        buttonHumanFall.grid(column=1, row=1,padx=2,pady=2)
        self.humanNumberLabel=Label(master=box3, text='No of Humans= %d' %self.model.var.numberOfHumans)
        self.humanNumberLabel.grid(column=2,row=1,padx=2,pady=2)
        self.entryVariable=Tkinter.DoubleVar()

        self.entry=Tkinter.Entry(box1,textvariable=self.entryVariable)

        self.entry.pack()
  
        self.entry.bind("<Return>", self.OnSetpointEnter)

        self.entryVariable.set(50.0)
#------------------------------------------------------------
        self.autoTuningStatus=False
        self.canvasFig=pltlib.figure(1)
        Fig=matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        FigSubPlot=Fig.add_subplot(111)
        x=[]
        y=[]
        y1=[]
        y2=[]
        self.job=[]
        self.X=[]
        self.allTemps=[]
        self.setPointArray=[]
        self.dampPos=[]
        self.line1,self.line2,self.line3,=FigSubPlot.plot(x,y,'r-',x,y1,'b--',x,y2,'g-')
        self.canvas=matplotlib.backends.backend_tkagg.FigureCanvasTkAgg(Fig,master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=5,row=0, sticky=E )
        self.canvas._tkcanvas.grid(column=3,row=0, sticky=E ,rowspan=9, padx=50)
        self.resizable(True,True)
        self.grid_columnconfigure(0,weight=1)
        self.update()
    

    def refreshFigure(self,x,y,setPoint,z):
        self.line1.set_data(x,y)
        self.line2.set_data(x,setPoint)
        self.line3.set_data(x,z)
        ax = self.canvas.figure.axes[0]
        ax.set_xlim(x.min(), x.max())
        k=max(y.max(),setPoint.max())
        ax.set_ylim(y.min()-10, k+10)
        self.bulbNumberLabel.config(text='No of Bulbs= %d' %self.model.var.numberOfBulbs)
        self.humanNumberLabel.config(text='No of Humans= %d' %self.model.var.numberOfHumans)        
        self.canvas.draw()
    
#--------------All the buttons Methods goes here---------------
    def OnStartClick(self):
        self.setPoint=self.entryVariable.get()
        self.pid.setSetPoint(self.entryVariable.get())
        self.pid.it=self.it
        self.pid.currentIt=self.it
        self.pid.riseStatus=True
        self.pid.trend=0
        self.compute()
    
    def OnSetpointEnter(self,event):
        self.OnStartClick()
       
 
    def startRelayTune(self):
        self.setPoint=self.entryVariable.get()
        self.ATR.setSetPoint(self.entryVariable.get())
        self.autoTuningStatus=True
        self.startIter=self.it
        self.autoTune()

    def incrementBulb(self):
        self.model.var.numberOfBulbs+=1
    def decrementBulb(self):
        if self.model.var.numberOfBulbs>0:
            self.model.var.numberOfBulbs-=1
    def incrementHuman(self):
        self.model.var.numberOfHumans+=1
    def decrementHuman(self):
        if self.model.var.numberOfHumans>0:
            self.model.var.numberOfHumans-=1
    def close_window(self):
        self.quit()
#---------------------------------------------------------------
    def autoTune(self):
        self.allTemps.append(self.ATR.getTRoom())
        self.TRoom=self.ATR.getTRoom()
        self.X.append(self.it) 
        self.setPointArray.append(self.setPoint)
        self.dampPos.append(self.model.var.percentOpen)
        X = np.array(self.X)
        Y = np.array(self.allTemps)
        Z= np.array(self.dampPos)
        setPointArray=np.array(self.setPointArray)
        self.refreshFigure(X,Y,setPointArray,Z)
        self.model.setI(self.it)
        self.cleanup()
        self.it+=1
        self.ATR.autoTune()
        if self.ATR.peaknumber<self.ATR.cycles:
            self.pid.setTRoom(self.ATR.getTRoom())
            self.after(1,self.autoTune)
        else: 
            self.autoTuningStatus=False
            self.pid.setSetPoint(self.entryVariable.get())
            print "Ku= ", self.ATR.ku," |  Pu= ", self.ATR.pu
            self.pid.setKu(self.ATR.ku)
            self.pid.setPu(self.ATR.pu) 
            self.pid.sumError=0
            self.compute()

    def compute(self):
        self.pid.Pid()
        self.allTemps.append(self.pid.getTRoom())
        self.TRoom=self.model.getTRoom()
        self.X.append(self.it)
        self.setPointArray.append(self.setPoint)
        self.dampPos.append(self.model.var.percentOpen)
        X = np.array(self.X)
        Y = np.array(self.allTemps)
        Z= np.array(0)#self.dampPos)
        setPointArray=np.array(self.setPointArray)
        self.refreshFigure(X,Y,setPointArray,Z)
        self.model.setI(self.it)
        self.it+=1
        self.cleanup()
        if self.autoTuningStatus==False:
            self.after(1,self.compute)
  
    def cleanup(self):
        self.secLabel.config(text="Time= "+str(self.it)+" secs")
        self.tempLabel.config(text="Room Temperature= %.2f" %     self.TRoom+" degree C")
        self.dampLabel.config(text="Damper Position= %.2f" %self.model.var.percentOpen+"%")
        
        

if __name__ == "__main__":
    MainWindow = App_Window(None)
    MainWindow.mainloop()


###GUI tkinter to put set point and do PID
#see reinforcement learning








