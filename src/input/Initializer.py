class Initializer:


#    spatial

    def getXleft(self):
        return -500.0e-4
    
    def getXright(self):
        return 500.0e-4

    def getYleft(self):
        return  0.0e-4
        
    def getYright(self):
        return 30.0e-4
    
    def getZleft(self):
        return 0.0
    
    def getZright(self):
        return 0.0
    
    def getXresolution(self):
        return 500
    
    def getYresolution(self):
        return 30
    
    def getZresolution(self):
        return 1

    def getXmpiDomainNum(self):
        return 2
    
    def getYmpiDomainNum(self):
        return 1
    
    def getZmpiDomainNum(self):
        return 1
    
#   time

    def getTimestep(self):
        return 0.5e-10 # 0.1 ns
    
    def getMaxTimestepsNum(self):
        return 10
    
#   output
#   output folder must be created manually
    def getOutputDir(self):
        return "output/"
    
    def getOutputFilenameTemplate(self):
        return "test_template_"
    
    def getEachTimestepsNum2WriteFile(self):
        return 1
    
    
#   physics
#   Al
    def getDensity(self, x, y, z):
        if x > 50.0e-4 and x < 70.0e-4 and  y > -20.0e-4 and y < 20.0e-4:
            return 1.2
        else:
            return 1.2


    def getVelocityX(self, x, y, z):
        return 0.0

    def getVelocityY(self, x, y, z):
        return 0.0

    def getVelocityZ(self, x, y, z):
        return 0.0
