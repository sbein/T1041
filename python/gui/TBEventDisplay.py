#!/usr/bin/env python
#-----------------------------------------------------------------------------
# File:        TBEventDisplay.py
# Description: TB 2014 Event Display
# Created:     09-May-2014 Harrison B. Prosper & Sam Bein
#              18-May-2014 HBP - simplify
#              19-May-2014 HBP - fix rewind/forward
#              09-Jun-2014 HBP - streamline calling of displays
#              12-Jun-2014 HBP - fix bugs in heatmap and surface plots
#-----------------------------------------------------------------------------
import sys, os, re
from ROOT import *
from time import ctime, sleep
from string import lower, replace, strip, split, joinfields, find
from glob import glob
from array import array
from gui.utils import *
from gui.TBFileReader import TBFileReader
from gui.TBWaveFormPlots import TracePlot, SurfacePlot
from gui.TBShashlikFaces import ShashlikHeatmap
from gui.TDCtiming import TDCtiming
from gui.TBFiberADC import FiberADC
from gui.TBDisplay3D import Display3D
#------------------------------------------------------------------------------
WIDTH        = 1000            # Width of GUI in pixels
HEIGHT       =  500            # Height of GUI in pixels
VERSION      = \
"""
TBEventDisplay.py %s
Python            %s
Root              %s
""" % ('v1.0',
       platform.python_version(),
       gROOT.GetVersion())
#-----------------------------------------------------------------------------
MINDELAY = 1.0 # (seconds) minimum delay between event displays

# Help
HELP = \
"""
The time has come the walrus said
"""

ABOUT = \
"""
%s
\tTB 2014
\te-mail: harry@hep.fsu.edu
\te-mail: samuel.bein@gmail.com
""" % VERSION

# read modes
R_REWIND  =-1
R_ONESHOT = 0
R_FORWARD = 1

DEBUG = 1
#-----------------------------------------------------------------------------
# (A) Root Graphical User Interfaces (GUI)
#
#   A Root GUI is a double tree of widgets (that is, graphical objects) with
#   which the user can interact in order to orchestrate a series of actions.
#   One tree describes the child-to-parent relationships, while the other
#   describes the parent-to-child relationships. The latter describes the
#   graphical layout of widgets. In the Root GUI system the two trees are not
#   isomorphic. For example, the child-to-parent relationship of a TGPopupMenu
#   is TGPopupMenu -> TGWindow, however, TGMenuBar -> TGPopupMenu is 
#   a typical parent-to-child relationship.
#
#   o A child-to-parent relationship is defined by the child widget when the
#     latter is created.
#
#   o A parent-to-child relationship, that is, a specific layout of a widget
#     within another, is defined by the parent after it has been created
#     using its AddFrame method in which the child is specified.
#
#   Each widget can emit one or more signals, usually triggered by some user
#   manipulation of, or action on, the widget. For example, clicking on a
#   tab of a  TGTab, that is, a notebook, will cause the notebook to emit
#   the signal "Selected(Int_t)" along with the identity of the selected tab.
#   Signals are connected to "Slots", that is, actions. This is how a user
#   GUI interaction can be made to yield one or more actions. Any signal can be
#   connected to any slot. Indeed, the relationship between signals and slots
#   can be many-many. In practice, a slot is modeled as a procedure such as
#   a method of the GUI class.
#  
#   In summary, a Root GUI is a (double) tree of widgets with which the user
#   can interact, whose signals---usually generated by user interactions---are
#   connected to slots, that is, actions modeled as methods of the GUI class.
#
# (B) This GUI
#
#   window                   (TGWindow)
#
#     main                      (TGMainFrame)
#
#       menuBar                    (TGMenuBar)
#         menuFile                    (TGPopupMenu)
#         menuEdit                    (TGPopupMenu)
#         menuEvent                   (TGPopupMenu)
#         menuHelp                    (TGPopupMenu)
#
#       vframe                     (TGVerticalFrame)
#         toolBar                    (TGToolBar)
#           nextButton                  (TGPictureButton)
#           previosButton               (TGPictureButton)
#
#         hframe                     (TGHorizontalFrame)
#           noteBook                    (TGTab)
#
#         statusBar                  (TGSTatusBar)
#-----------------------------------------------------------------------------
class TBEventDisplay:
    """
    gui = TBEventDisplay(title)
    """

    def __init__(self, title, filename=None, width=WIDTH, height=HEIGHT):

        # Initial directory for open file dialog
        self.openDir  = os.environ['PWD']
        self.filename = filename

        #-------------------------------------------------------------------
        # Create main frame
        #-------------------------------------------------------------------
        # Establish a connection between the main frame's "CloseWindow()"
        # signal and the GUI's "close" slot, modeled as a method.
        # When the main frame issues the signal CloseWindow() this
        # triggers a call to the close method of this class.
        #-------------------------------------------------------------------
        self.root = gClient.GetRoot()
        self.main = TGMainFrame(self.root, width, height)
        self.main.SetWindowName(title)
        self.main.SetCleanup(kDeepCleanup)
        self.connection = Connection(self.main, "CloseWindow()",
                         self,      "close")

        self.util = Util()
        #-------------------------------------------------------------------
        # Create menu bar
        #-------------------------------------------------------------------
        self.menuBar = MenuBar(self, self.main)

        self.menuBar.Add('File',
                 [('&Open',  'openFile'),
                  ('&Close', 'closeFile'),
                  0,
                  ('E&xit',  'exit')])

        self.menuBar.Add('Edit',
                 [('&Undo',  'undo')])


        self.menuBar.Add('Event',
                 [('&Next',     'nextEvent'),
                  ('&Previous', 'previousEvent'),
                  ('&Goto',     'gotoEvent'),
                0,
                  ('Set ADC cut', 'setADCcut'),
                  ('Set delay',   'setDelay')])

        self.menuBar.Add('Help',
                 [('About', 'about'),
                  ('Usage', 'usage')])

        #-------------------------------------------------------------------
        # Add vertical frame to the main frame to contain toolbar, notebook
        # and status window
        #-------------------------------------------------------------------
        self.vframe = TGVerticalFrame(self.main, 1, 1)
        self.main.AddFrame(self.vframe, TOP_X_Y)

        #-------------------------------------------------------------------
        # Add horizontal frame to contain toolbar
        #-------------------------------------------------------------------
        self.toolBar = TGHorizontalFrame(self.vframe)
        self.vframe.AddFrame(self.toolBar, TOP_X)

        # Add picture buttons        
        #self.refreshButton = PictureButton(self, self.toolBar,
        #                picture='Button-Refresh.png',
        #                method='refreshFile',
        #                text='refresh file for new events')

        
        self.enchiladaButton = PictureButton(self, self.toolBar,
                        picture='Enchilada.jpg',
                        method='wholeEnchilada',
                        text='process the whole enchilada')

        self.nextButton = PictureButton(self, self.toolBar,
                        picture='GoForward.gif',
                        method='nextEvent',
                        text='go to next event')

        self.forwardButton = PictureButton(self, self.toolBar,
                           picture='Forward.png',
                           method='forwardPlayer',
                           text='foward event player')

        self.stopButton = PictureButton(self, self.toolBar,
                        picture='Stop.png',
                        method='stopPlayer',
                        text='stop event player')

        self.rewindButton = PictureButton(self, self.toolBar,
                          picture='Rewind.png',
                          method='rewindPlayer',
                          text='rewind event player')

        self.previousButton = PictureButton(self, self.toolBar,
                            picture='GoBack.gif',
                            method='previousEvent',
                            text='go to previous event')
        
        self.snapCanvasButton = PictureButton(self, self.toolBar,
                        picture='Camera.png',
                        method='snapCanvas',
                        text='snap this canvas')

        self.cycleButton = PictureButton(self, self.toolBar,
                        picture='Bicycle.jpg',
                        method='cycleTabs',
                        text='snap this canvas')

        self.accumulateButton = CheckButton(self, self.toolBar,
                        hotstring='Accumulate',
                        method='toggleAccumulate',
                        text='Accumulate')


        #-------------------------------------------------------------------
        # Add a notebook with multiple pages
        #-------------------------------------------------------------------

    
        
        #general util struct members:
        self.util.accumulate = False   
        self.util.stealthmode = False  
        self.util.filename = self.filename
        #WC trace bools:
        self.util.showBoard112 = True
        self.util.showBoard113 = True
        self.util.showBoard115 = True
        self.util.showBoard116 = True
        #WC bools:
        self.util.WC_showIThits = True
        self.util.WC_showQhits = True
        #track struct members:
        self.util.x1hit = 10
        self.util.x2hit = 1
        self.util.y1hit = 10
        self.util.y2hit = 4
        #3D bools
        self.util._3D_showWC1 = True
        self.util._3D_showWC2 = True
        self.util._3D_isolateClusters = False

        WFbuttons = [('Board 112','toggleShowBoard112','toggleShowBoard112',True),('Board 113','toggleShowBoard113','toggleShowBoard113',True),('Board 115','toggleShowBoard115','toggleShowBoard115',True),('Board 116','toggleShowBoard116','toggleShowBoard116',True)]
        WCbuttons = [('In-time hits','toggleShowIThits','toggleShowIThits', True),('Quality hits','toggleShowQhits','toggleShowQhits', True)]
        ThreeDbuttons = [('show WC1','toggleShowWC1','toggleShowWC1', True),('show WC2','toggleShowWC2','toggleShowWC2', True),('isolate cluster','toggleIsolateClusters','toggleIsolateCluster', False)]
        
        
        
        self.reraw = True
        self.shutterOpen = False
        self.pageName = 'default'
                 
        self.noteBook = NoteBook(self, self.vframe, 'setPage', width, height)
        # Add pages 
        self.display = {}
        for pageName, constructor, buttons in [('WF traces', 'TracePlot(canvas)', WFbuttons),
                          ('ADC heatmap',   'ShashlikHeatmap(canvas)', None),
                          ('ADC fibers', 'FiberADC(canvas)', None), 
                          ('Wire chambers', 'WCPlanes(canvas)', WCbuttons),
                          ('TDC timing','TDCtiming(canvas)', None),
                          ('3D display',    'Display3D(page)', ThreeDbuttons)]:
            self.noteBook.Add(pageName, buttons)
            self.noteBook.SetPage(pageName)
            page = self.noteBook.page
            canvas = page.canvas
            print '==> building display: %s' % pageName
            self.display[pageName] = eval(constructor)
            
        #-------------------------------------------------------------------
        # Create a status bar, divided into two parts
        #-------------------------------------------------------------------
        self.statusBar = TGStatusBar(self.vframe, 1, 1, kDoubleBorder)
        self.statusBar.SetHeight(22)
        status_parts = array('i')
        status_parts.append(18)
        status_parts.append(18)
        status_parts.append(24)
        status_parts.append(20)
        status_parts.append(20)
        self.statusBar.SetParts(status_parts, len(status_parts))
        self.progressBar = ProgressBar(self, self.statusBar)
        self.vframe.AddFrame(self.statusBar, TOP_X)

    
        # Initial state
        self.wcplanes = self.display['Wire chambers']	
        self.ADCcut  = 500
        self.nevents = 0
        self.eventNumber = -1
        self.DELAY  = int(1000*MINDELAY)
        self.mutex  = TMutex(kFALSE)
        self.timer  = TTimer()
        self.timerConnection = Connection(self.timer, 'Timeout()',
                                          self, 'managePlayer')
        
        self.DEBUG  = DEBUG
        self.DEBUG_COUNT = 0

        # Initialize layout        
        self.main.MapSubwindows()
        self.main.Resize()
        self.main.MapWindow()
	
        
        #DEBUG
        # to debug a display uncomment next line
        #filename = "data/test.root"

        if filename != None: self.__openFile(filename)
	else: self.__openFile('data/test.root')	

        #DEBUG
        # to debug a display uncomment next two lines
        self.noteBook.SetPage('ADC heatmap')
        self.displayEvent()
        
    def __del__(self):
        pass


    def debug(self, message):
        if self.DEBUG < 1: return
        self.DEBUG_COUNT += 1
        print "%10d> %s" % (self.DEBUG_COUNT, message)

    #-----------------------------------------------------------------------
    #	M E T H O D S
    #-----------------------------------------------------------------------

    #	S L O T S    (that is, callbacks)

    def openFile(self):
        dialog = Dialog(gClient.GetRoot(), self.main)
        self.filename = dialog.SelectFile(kFDOpen, self.openDir)
        self.util.filename = self.filename
        self.openDir = dialog.IniDir()
        if self.filename[-5:] != '.root':
            dialog.ShowText("Oops!",
                    "Please select a root file",
                    230, 30)
            return
        self.__openFile(self.filename)


    def __openFile(self, filename):
        self.filename = filename
        self.util.filename = self.filename
        self.closeFile()		
        self.reader = TBFileReader(filename)
        self.nevents= self.reader.entries()
        self.statusBar.SetText('events: %d' % self.nevents, 0)
        self.statusBar.SetText(self.filename[self.filename.rfind('/')+1:][-32:], 2)
        self.eventNumber = -1
        self.util.eventNumber = -1
        self.nextEvent()
        self.progressBar.SetRange(0, self.nevents)
        self.progressBar.SetPosition(self.eventNumber)
        self.wcplanes.CacheWCMeans("meanfile.txt", filename)
        self.util.tableX = self.reader.tableX
        self.util.tableY = self.reader.tableY
        self.statusBar.SetText('table(x, y) = (%d, %d) ' % (self.util.tableX, self.util.tableY), 1)
        self.filetime = time.ctime(os.path.getctime(filename))

        
    def refreshFile(self):
        try:
            t = self.filetime
        except:
            return
        if self.filetime == time.ctime(os.path.getctime(self.filename)):
            print "file up to date"
            return                 
        else:      
            print "refreshing"
            eventNumber = self.eventNumber
            self.closeFile()	            
            self.reader = TBFileReader(self.filename)
            self.nevents= self.reader.entries()
            self.statusBar.SetText('event: %d / %d' % (self.eventNumber, self.nevents-1), 0)
            self.filetime = time.ctime(os.path.getctime(self.filename))
            self.wcplanes.CacheWCMeans("meanfile.txt", self.filename)
            print "in refresh fashion:"
            self.eventNumber = eventNumber
            self.util.eventNumber = eventNumber
        
    def closeFile(self):
        try:
            if self.reader.file().IsOpen():
                #gSystem.Exit()
                self.reader.file().Close()
                del self.reader
        except:
            pass

    def setPage(self, id):
        if id==0:
            self.pageName = 'Traces'
        if id==1:
            self.pageName = 'Heatmap'
        if id==2:
            self.pageName = 'Fibers'
        if id==3:
            self.pageName = 'WC'
        if id==4:
            self.pageName = 'TDC'
        if id==5:
            self.pageName = '3D_'
        self.noteBook.SetPage(id)
        if self.eventNumber >= 0:
            self.displayEvent()

    def nextEvent(self):
        self.debug("begin:nextEvent")
        if self.eventNumber > self.nevents-2:
            self.eventNumber = 0
        self.readEvent(R_FORWARD)
        self.displayEvent()
        self.debug("end:nextEvent")			

    def previousEvent(self):
        self.debug('begin:previousEvent')
        if self.eventNumber < 1:
            self.eventNumber = self.nevents 
        self.readEvent(R_REWIND)
        self.displayEvent()
        self.debug('end:previousEvent')

    def gotoEvent(self):
        self.debug('begin:gotoEvent')
        from string import atoi
        dialog = Dialog(gClient.GetRoot(), self.main)
        self.eventNumber = atoi(dialog.GetInput('Goto event %d - %d' % \
                            (0, self.nevents-1),
                            '0'))
        self.eventNumber = max(self.eventNumber, 0)
        self.eventNumber = min(self.eventNumber, self.nevents-1)

        # do a one-shot read
        self.readEvent(R_ONESHOT)
        self.displayEvent()
        self.debug('end:gotoEvent')

    def forwardPlayer(self):
        self.debug('begin:forwardPlayer')
        self.mutex.Lock()
        self.forward = True
        self.timer.Start(self.DELAY, kFALSE)
        self.mutex.UnLock()
        self.debug('end:forwardPlayer')

    def cycleTabs(self):
        self.debug('begin:cycleTabs')
        self.mutex.Lock()
        if self.noteBook.pageNumber == 5:#len(self.noteBook.pages)-1:
            self.noteBook.pageNumber = 0
        self.timer.Start(self.DELAY*5, kFALSE)
        self.noteBook.pageNumber+=1
        self.noteBook.SetPage(self.noteBook.pageNumber)
        self.displayEvent()
        self.forward = True
        self.cycle = True
        self.mutex.UnLock()
        self.debug('end:cycleTabs')
        
        
    def rewindPlayer(self):
        self.debug('begin:rewindPlayer')
        self.mutex.Lock()
        self.forward = False
        self.timer.Start(self.DELAY, kFALSE)
        self.mutex.UnLock()
        self.debug('end:rewindPlayer')

    def stopPlayer(self):
        self.debug('begin:stopPlayer')
        self.mutex.Lock()
        self.timer.Stop()
        self.cycle = False
        self.mutex.UnLock()
        self.debug('end:stopPlayer - STOP REQUESTED')

    def managePlayer(self):
        if self.cycle:
            self.cycleTabs()
        if self.forward:
            self.nextEvent()
        else:
            self.previousEvent()

    def toggleAccumulate(self):
        self.debug("toggle:accumulate")
        self.util.accumulate = not self.util.accumulate

    def wholeEnchilada(self):
        self.debug("begin:wholeEnchilada")
        self.util.accumulate = True
        self.util.accumulateButton.SetState(True) 
        self.util.stealthmode = True
        self.eventNumber = 0
        while self.eventNumber<self.nevents-2:
            self.readEvent(R_FORWARD)
            self.displayEvent()
            self.statusBar.Draw()
        self.util.stealthmode = False
        self.readEvent(R_FORWARD)
        self.displayEvent()     
        self.eventNumber = 0
        self.debug("end:wholeEnchilada")
        
    def snapCanvas(self):
        self.debug("begin:snapCanvas")
        self.shutterOpen = True
        self.displayEvent()
        self.shutterOpen = False
        self.debug("end:snapCanvas")

    def toggleShowBoard112(self):
            self.debug("begin:toggleShowBoard112")
            self.util.showBoard112 = not self.util.showBoard112
            self.redraw = True
            self.displayEvent()
            self.debug("end:toggleShowBoard112")

    def toggleShowBoard113(self):
            self.debug("begin:toggleShowBoard113")
            self.util.showBoard113 = not self.util.showBoard113
            self.redraw = True
            self.displayEvent()
            self.debug("end:toggleShowBoard113")
            
    def toggleShowBoard115(self):
            self.debug("begin:toggleShowBoard115")
            self.util.showBoard115 = not self.util.showBoard115
            self.redraw = True
            self.displayEvent()
            self.debug("end:toggleShowBoard115")

    def toggleShowBoard116(self):
            self.debug("begin:toggleShowBoard116")
            self.util.showBoard116 = not self.util.showBoard116
            self.redraw = True
            self.displayEvent()
            self.debug("end:toggleShowBoard116")

    def toggleShowIThits(self):
        self.debug("begin:ShowIThits")
        self.util.WC_showIThits = not self.util.WC_showIThits
        self.redraw = True
        self.displayEvent()
        self.debug("end:ShowIThits")

    def toggleShowQhits(self):
        self.debug("begin:ShowQhits")
        self.util.WC_showQhits = not self.util.WC_showQhits
        self.redraw = True
        self.displayEvent()
        self.debug("end:ShowQhits")
        

    def toggleShowWC2(self):
        self.debug("begin:toggleShowWC2")
        self.util._3D_showWC2 = not self.util._3D_showWC2
        self.redraw = True
        self.displayEvent()
        self.debug("end:toggleShowWC2")

        
    def toggleShowWC1(self):
        self.debug("begin:toggleShowWC2")
        self.util._3D_showWC1 = not self.util._3D_showWC1
        self.redraw = True
        self.displayEvent()
        self.debug("end:toggleShowWC1")

        
    def toggleIsolateClusters(self):
            self.debug("begin:toggleIsolateClusters")
            self.util._3D_isolateClusters = not self.util._3D_isolateClusters
            self.redraw = True
            self.displayEvent()
            self.debug("end:toggleIsolateClusters")

        
 
        
    def setDelay(self):
        from string import atof
        dialog = Dialog(gClient.GetRoot(), self.main)
        seconds= atof(dialog.GetInput('Enter delay in seconds', '2.0'))
        self.playerDelay = max(MINDELAY, int(1000*seconds))
        self.statusBar.SetText('delay set to: %8.1f s' % seconds, 1)

    def setADCcut(self):
        from string import atoi
        dialog = Dialog(gClient.GetRoot(), self.main)
        self.ADCcut = atoi(dialog.GetInput('Enter ADC cut', '50'))
        self.statusBar.SetText('table(x, y) = (%d, %d) ' % (self.util.tableX, self.util.tableY), 1)

    def close(self):
        gSystem.Abort()
        gApplication.Terminate(0)
        #gApplication.Close()

    def usage(self):
        dialog = Dialog(gClient.GetRoot(), self.main)
        dialog.SetText('Not done', 'Sorry!', 230, 30)
        #dialog.SetText('Help', HELP)

    def about(self):
        dialog = Dialog(gClient.GetRoot(), self.main)
        dialog.SetText('Not done', 'Sorry!', 230, 30)
        #dialog.SetText('About', ABOUT)

    def exit(self):
        gEve.CloseEveWindow()
        self.close()

    def notdone(self):
        dialog = Dialog(gClient.GetRoot(), self.main)
        dialog.SetText('Not done', 'Sorry!', 230, 30)

    def run(self):
        gApplication.Run()


    #-----------------------------------------------------------------------
    # O T H E R   M E T H O D S
    #-----------------------------------------------------------------------

    # read events until we find one with a channel above specified ADC cut

    def readEvent(self, which=R_ONESHOT):
        self.debug("begin:readEvent")

        try:
            reader = self.reader
        except:
            dialog = Dialog(gClient.GetRoot(), self.main)
            dialog.SetText('Oops!', 'First open a root file',
                           230, 24)
            self.debug("end:readEvent")
            return


        # cache previous event number
        self.eventNumberPrev = self.eventNumber

        # loop over events and apply ADC cut

        if   which == R_ONESHOT:
            self.statusBar.SetText('event: %d / %d' % (self.eventNumber, self.nevents-1),
                                   0)
            
            self.reader.read(self.eventNumber)

            self.statusBar.SetText('table(x, y) = (%d, %d) ' % (self.util.tableX, self.util.tableY), 1)

        elif which == R_FORWARD:
            while self.eventNumber < self.nevents-1:
                self.eventNumber += 1
                self.statusBar.SetText('event: %d / %d' % (self.eventNumber, self.nevents-1),
                                       )
                self.reader.read(self.eventNumber)

                ADCmax =  self.reader.maxPadeADC()
                #self.statusBar.SetText('max(ADC): %d' % ADCmax, 1)
                self.statusBar.SetText('table(x, y) = (%d, %d) ' % (self.util.tableX, self.util.tableY), 1)
                if ADCmax < self.ADCcut: continue
                break
        else:
            while self.eventNumber > 0:
                self.eventNumber -= 1
                self.statusBar.SetText('event: %d / %d' % (self.eventNumber, self.nevents-1),
                                       0)
                self.reader.read(self.eventNumber)

                ADCmax =  self.reader.maxPadeADC()
                #self.statusBar.SetText('max(ADC): %d' % ADCmax, 1)
                self.statusBar.SetText('table(x, y) = (%d, %d) ' % (self.util.tableX, self.util.tableY), 1)
                if ADCmax < self.ADCcut: continue				
                break			


        if self.eventNumber <= 0 or self.eventNumber >= self.nevents-1:
            self.stopPlayer()

        # Force a re-drawing of pages of notebook when a page is made
        # visible
        keys = self.noteBook.pages.keys()
        for key in keys:
            self.noteBook.pages[key].redraw = True
        self.util.eventNumber = self.eventNumber

        self.debug("end:readEvent")
    #-----------------------------------------------------------------------
    def displayEvent(self):
        self.progressBar.Reset()
        self.progressBar.SetPosition(self.eventNumber)
        pageNumber = self.noteBook.pageNumber
        page = self.noteBook.pages[pageNumber]
        self.debug("begin:displayEvent - %s" % page.name)
        if self.shutterOpen:
            if not os.path.exists("pdfs_"+self.filename[self.filename.rfind('/')+1:-5]):
                os.system('mkdir pdfs_'+self.filename[self.filename.rfind('/')+1:-5])

            if '3D' in page.name:
                pdfname = 'pdfs_'+self.filename[self.filename.rfind('/')+1:-5]+\
                  '/'+self.pageName+str(self.eventNumber)+'.png'
                gEve.GetDefaultGLViewer().SavePictureScale(pdfname,1.0)
            else:
                pdfname = 'pdfs_'+self.filename[self.filename.rfind('/')+1:-5]+\
                  '/'+self.pageName+str(self.eventNumber)+'.pdf'
                page.canvas.Print(pdfname)
            
            page.redraw = False
            return
        if not page.redraw and not self.shutterOpen and not self.redraw:
            self.debug("end:displayEvent - DO NOTHING")		
            return
        self.refreshFile()

        if '3D' not in page.name:
            self.display[page.name].Draw(self.reader.event(), self.util)
        else:
            self.display['ADC heatmap'].Draw(self.reader.event(), self.util)
            self.util.stealthmode = True
            self.display['Wire chambers'].Draw(self.reader.event(), self.util)
            self.util.stealthmode = False
            self.display[page.name].Draw(self.reader.event(), self.util)
            #if 'TDC' in page.name:
            #self.tdcfile.Close()
        print "displaying with self.util.eventNumber = "+str(self.util.eventNumber)
        self.redraw = False
        page.redraw = False
        self.debug("end:displayEvent")		
