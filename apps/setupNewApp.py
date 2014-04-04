import sys, os

from os.path import join


def main():
    
    if len(sys.argv) == 1:
        print "Usage: python %s [app name]" % sys.argv[0]
        return
        
    appName = sys.argv[1]
    
    cwd = os.getcwd()    
    
    absAppPath = join(cwd, appName)
    
    if os.path.exists(absAppPath):
        print "App %s already exists." % appName
        return

    os.mkdir(absAppPath)
    
    infilesPath = join(absAppPath, "infiles")

    os.mkdir(infilesPath)


    configFileAbsPath = join(infilesPath, "%s.cfg" % appName)    
    
    appProFileAbsPath = join(absAppPath, "%s.pro" % appName)
    
    appMainFileAbsPath = join(absAppPath, "%smain.cpp" % appName)


    with open(appProFileAbsPath, 'w') as appProFile:
        
        with open(join(cwd, "defaults", "default.pro.bones"), 'r') as defaultFile:
            
            appProFile.write(defaultFile.read().replace("__name__", appName))
            
    with open(appMainFileAbsPath, 'w') as appMainFile:
            
        with open(join(cwd, "defaults", "defaultmain.cpp.bones"), 'r') as defaultFile:
            
            appMainFile.write(defaultFile.read().replace("__name__", appName))
           
    with open(configFileAbsPath, 'w') as appConfigFile:
            
        with open(join(cwd, "defaults", "defaultconfig.cfg.bones"), 'r') as defaultFile:
            
            appConfigFile.write(defaultFile.read().replace("__name__", appName))
          
    
    appsProFileAbspath = join(cwd, "apps.pro")
    
    with open(appsProFileAbspath, 'r') as appsProFile:
        
        newAppsProFile = appsProFile.read().replace("#__next_app__", "\\\n           %s #__next_app__" % appName)
        
    with open(appsProFileAbspath, 'w') as appsProFile:
        
        appsProFile.write(newAppsProFile)
        
                
            
        
if __name__ == "__main__":
    main()                
    