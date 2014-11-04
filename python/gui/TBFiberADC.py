import os, sys, glob, random
from ROOT import *
from time import sleep
from TBUtils import *
import gui.utils

Nsmp = PadeChannel.N_PADE_SAMPLES / 2
Nch = 128

class FiberADC:
	def __init__(self, canvas):
		self.canvas = canvas
	def __del__(self):
		pass

	def Draw(self, event, spill, rechits, tracks, util):
		MakePlots(self, self.canvas, event, rechits, util)

        
		
def MakePlots(object, c1, event, rechits, util):
    try:
        o = object.hMapADCvsFiber
        p = object.hMapCHI2vsFiber

    except:

        c1.Divide(1,2)
        gStyle.SetPalette(1)
        gStyle.SetOptStat(0)
        c1.GetPad(1).SetLogz()
        c1.GetPad(1).SetLogy()

        object.hMapADCvsFiber = TH2F('hMapADCfiber','ADC Samples Vs. Channel', Nch,0,Nch,3000,0,3000)
        
        object.hMapCHI2vsFiber = TH2F('hMapADCvsFiber','Chi2 Vs. Channel', Nch,0,Nch,25,0,1000)
        HistoSamStyle(object.hMapADCvsFiber)
        HistoSamStyle(object.hMapCHI2vsFiber)

    hMapADCvsFiber = object.hMapADCvsFiber
    hMapCHI2vsFiber = object.hMapCHI2vsFiber

    if not util.accumulate:
        hMapADCvsFiber.Reset()
        hMapCHI2vsFiber.Reset()


    print "size of rechits is ", rechits.size()
    if util.FADC_showRecHits:
        for rh in range(rechits.size()):
            pid = rechits[rh].ChannelIndex()
            for ch in range(0,Nch):
                pade = event.GetPadeChan(ch)
                pid2 = pade.GetChannelIndex()
                if pid2==pid:
                    chyes = ch
                    break
            for iwf in range(0,60):
                number = pade.GetWform()[iwf] - pade.GetPedestal()
                hMapADCvsFiber.Fill(pid,number)
            chi2= min(999,rechits[rh].Chi2())
            ndof = rechits[rh].Ndof()
            hMapCHI2vsFiber.Fill(pid,chi2/ndof)

    if util.FADC_showAllHits:
        hit = TBRecHit()
        nSigmaCut=1
        for ch in range(0,Nch):
            pade = event.GetPadeChan(ch)
            pid = pade.GetChannelIndex()
            hit.Init(pade, nSigmaCut)
            if not (hit.Status() and TBRecHit.kZSP):
                for iwf in range(0,60):
                    number = pade.GetWform()[iwf] - pade.GetPedestal()
                    hMapADCvsFiber.Fill(pid,number)
                chi2= min(999,hit.Chi2())
                ndof = hit.Ndof()
                if (ndof!=0):
                    hMapCHI2vsFiber.Fill(pid,1.0*chi2/ndof)
        

    if not util.stealthmode:
        gStyle.SetTitleSize(.08,"t"); 
        c1.cd(1)
        c1.GetPad(1).SetLogz()
        c1.GetPad(1).SetLogy()
        hMapADCvsFiber.Draw('COLZ')
        c1.cd(2)
        hMapCHI2vsFiber.Draw('COLZ')
        c1.Update()

    
    
            
        
                
        
def HistoSamStyle(histo):
    LabelSize = 0.055
    lwidth = 2
    histo.GetXaxis().SetTitle('channel    \t\t  ')
    histo.GetXaxis().SetTitleSize(1.4*LabelSize)
    histo.GetXaxis().SetTitleOffset(0.45)
    histo.GetXaxis().SetLabelSize(1.2*LabelSize)
    histo.GetYaxis().SetLabelSize(1.06*LabelSize)
    histo.GetYaxis().SetTitleOffset(1)
    histo.GetZaxis().SetLabelSize(0.7*LabelSize)

def getWFmax(pade):
	max = 0
	wform = pade.GetWform()
	for iwf in range(0,PadeChannel().N_PADE_SAMPLES): 
		number = pade.GetWform()[iwf]
		if number > max:
			max = number
	return max

