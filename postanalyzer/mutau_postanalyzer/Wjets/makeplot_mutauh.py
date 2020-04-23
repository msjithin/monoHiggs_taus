#!/usr/bin/env python
import ROOT
import re
from array import array
from sys import argv
import ROOT as Root
import csv
from math import sqrt
from math import pi
import datetime

now = datetime.datetime.now()

firstarg=argv[1]
secondarg=argv[2]
thirdarg = argv[3]
fourtharg = ""
fiftharg = ""
if len(argv) > 4 :
  fourtharg=argv[4]

print "length = ",len(argv) 
thirdarg = thirdarg + " " + fourtharg
if firstarg != 'Events_level_' :
  histoname = firstarg+secondarg
else :
  histoname = firstarg
print "Item / histogram name = ", histoname
print "thirdarg = ", thirdarg

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        if argv[0][0] == '-':  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts


myargs = getopts(argv)
if '-logYaxis' in myargs:  # Example usage.
  yaxisLog = 1
else :
  yaxisLog = 0

if '-noLegend' in myargs:  # Example usage.
  noLegend = _mc
else :
  noLegend = 0

if (firstarg == "h_dPhi_" or firstarg =="Tau_phi_" or firstarg =="Muon_phi_" or firstarg =="dphi_muMet_" or firstarg =="dphi_tauMet"):
  piAxis = 1
else :
  piAxis = 0



def add_lumi():
    lowX=0.50
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.10, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.04)
    lumi.SetTextFont (   42 )
    lumi.AddText("#mu#tau_{h}   2017, 41.521 fb^{-1} (13 TeV)")
    return lumi

def add_CMS():
    lowX=0.21
    lowY=0.70
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi

def add_Preliminary():
    lowX=0.21
    lowY=0.63
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(52)
    lumi.SetTextSize(0.06)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("Preliminary")
    return lumi

def make_legend():
        output = ROOT.TLegend(0.70, 0.70, 0.92, 0.84, "", "brNDC")
        #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
        output.SetLineWidth(1)
        output.SetLineStyle(1)
        output.SetFillStyle(0)
        output.SetBorderSize(1)
        output.SetTextFont(62)
        return output

ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetLineWidth(1)
ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("canvas","",0,0,600,600)
c.cd()

OutFile=ROOT.TFile("f_mutau_initial.root","r")
#OutFile_fbkg=ROOT.TFile("ff_method/f_mutau_fakebkg.root","r")

adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

categories=[firstarg,firstarg] 
dirName=[firstarg+"/","qcd/","bkg/"+firstarg+"/"] 
ncat=1


if firstarg=="MET" :  # Example usage.
  dirName[0]="MET_0_/"
if firstarg=="Cutflow": 
  dirName[0]="Cutflow_/"

for i in range (0,ncat):
    WJets=OutFile.Get(dirName[0]+"WJetsToLNu_inc_"+histoname)
    

#    WJets=OutFile.Get(dirName[0]+"WJetsToLNu_"+histoname)
#    WJets_100To200=OutFile.Get(dirName[0]+"WJetsToLNu_HT100To200_"+histoname)
#    WJets_200To400=OutFile.Get(dirName[0]+"WJetsToLNu_HT200To400_"+histoname)
#    WJets_400To600=OutFile.Get(dirName[0]+"WJetsToLNu_HT400To600_"+histoname)
#    WJets_600To800=OutFile.Get(dirName[0]+"WJetsToLNu_HT600To800_"+histoname)
#    WJets_800To1200=OutFile.Get(dirName[0]+"WJetsToLNu_HT800To1200_"+histoname)
#    WJets_1200To2500=OutFile.Get(dirName[0]+"WJetsToLNu_HT1200To2500_"+histoname)
#    WJets_2500ToInf=OutFile.Get(dirName[0]+"WJetsToLNu_HT2500ToInf_"+histoname)
#    nbins_ =  WJets.GetNbinsX()
#    print "WJets inclusive events   ", WJets.Integral(0,nbins_+1)
#    print "WJets 100-200 events   ", WJets_100To200.Integral(0,nbins_+1)
#    print "WJets 200-400 events   ", WJets_200To400.Integral(0,nbins_+1)
#    print "WJets 400-600 events   ", WJets_400To600.Integral(0,nbins_+1)
#    print "WJets 600-800 events   ", WJets_600To800.Integral(0,nbins_+1)
#    print "WJets 800-1200 events   ", WJets_800To1200.Integral(0,nbins_+1)
#    print "WJets 1200-2500 events   ", WJets_1200To2500.Integral(0,nbins_+1)
#    print "WJets 2500-Inf events   ", WJets_2500ToInf.Integral(0,nbins_+1)

#    WJets.Add(WJets_100To200)
#    WJets.Add(WJets_200To400)
#    WJets.Add(WJets_400To600)
 #   WJets.Add(WJets_600To800)
 #   WJets.Add(WJets_800To1200)
 #   WJets.Add(WJets_1200To2500)
 #   WJets.Add(WJets_2500ToInf)
 #   WJets.Add(EWKW)



    WJets.SetFillColor(ROOT.TColor.GetColor("#ffc866"))


    #EWKZ.SetLineColor(1)
    WJets.SetLineColor(1)
    
    stack=ROOT.THStack("stack","stack")
   # stack.Add(ZJetsToNuNu)
  #  stack.Add(WWW)
    stack.Add(WJets)
   # stack.Add(VV)
#    stack.Add(TTJ)
   # stack.Add(F_bkg)
   # stack.Add(TTJ)
   # stack.Add(ZTT)
   
    


    WJets.GetXaxis().SetTitle("")
    WJets.GetXaxis().SetTitleSize(0)


    if piAxis == 1:
      WJets.GetXaxis().SetNdivisions(-405)
    elif histoname=='Events_level_':
      WJets.GetXaxis().SetNdivisions(-115)
    else :
      WJets.GetXaxis().SetNdivisions(-405)

    WJets.GetXaxis().SetLabelSize(0.05)
    WJets.GetYaxis().SetLabelFont(42)
    WJets.GetYaxis().SetLabelOffset(0.01)
    WJets.GetYaxis().SetLabelSize(0.05)
    WJets.GetYaxis().SetTitleSize(0.05)
    WJets.GetYaxis().SetTitleOffset(1.1)
    WJets.GetXaxis().SetTitleSize(0.05)
    #WJets.GetXaxis().SetTitleOffset(1.04)
    
    WJets.SetTitle("")
    WJets.GetYaxis().SetTitle("Events")
    WJets.GetXaxis().SetTitle(thirdarg)
        
    errorBand=WJets.Clone()


    errorBand.SetMarkerSize(0)
    errorBand.SetFillColor(1)
    errorBand.SetFillStyle(3001)
    errorBand.SetLineWidth(1)

pad1 = ROOT.TPad("pad1","pad1",0,0.1,1,1)
pad1.Draw()
pad1.cd()
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetBorderSize(1)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad1.SetGridx()
pad1.SetLeftMargin(0.18)
pad1.SetRightMargin(0.05)
pad1.SetTopMargin(0.122)
pad1.SetBottomMargin(0.122)
pad1.SetFrameFillStyle(0)
pad1.SetFrameLineStyle(0)
pad1.SetFrameLineWidth(1)
pad1.SetFrameBorderMode(0)
pad1.SetFrameBorderSize(1)
if yaxisLog == 1 :
  pad1.SetLogy()  

#for k in range(1,Data.GetSize()-1):
    #s=ZH125.GetBinContent(k)
    ##b=VV.GetBinContent(k)+Fake.GetBinContent(k)
    #b=VV.GetBinContent(k)+ggH125.GetBinContent(k)+W.GetBinContent(k)
    #if (b<0):
	#b=0.000001
    #if (s/(0.00001+0.05*s+b)**0.5 > 0.8):
	#Data.SetBinContent(k,-1)
	#Data.SetBinError(k,-1)
#WJets.SetMarkerStyle(20)

#WJets.SetMarkerSize(0.75)
WJets.GetXaxis().SetTitle(thirdarg)
WJets.GetYaxis().SetTitle("Events")
if yaxisLog == 1 :
  WJets.SetMaximum(1000*max(WJets.GetMaximum(),stack.GetMaximum()))
  WJets.SetMinimum(0.1)
else :
  WJets.SetMaximum(1.35*max(WJets.GetMaximum(),stack.GetMaximum()))
  WJets.SetMinimum(0.0)
#stack.SetLineColor(9)
WJets.Draw()
stack.Draw("histsame")

#errorBand.Draw("e2same")
#ZH125_mc.Draw("histsame")

#WJets.Draw("e1same")
#errorBand.Draw("e2same")

legende=make_legend()
legende.AddEntry(WJets,"WJets","f")
#legende.AddEntry(WJets,"Data runB","elp")
#legende.AddEntry(WJets,"Data runC","elp")
#legende.AddEntry(WJets,"Data runD","elp")
#legende.AddEntry(WJets,"Data runE","elp")
#legende.AddEntry(WJets,"Data runF","elp")

#legende.AddEntry(WJets,"fake background","f")
#legende.AddEntry(W_mc,"W","f")

if noLegend == 0:
  legende.Draw()

l1=add_lumi()
l1.Draw("same")
l2=add_CMS()
l2.Draw("same")
l3=add_Preliminary()
l3.Draw("same")

pad1.RedrawAxis()

categ  = ROOT.TPaveText(0.21, 0.5+0.013, 0.43, 0.70+0.155, "NDC")
categ.SetBorderSize(   0 )
categ.SetFillStyle(    0 )
categ.SetTextAlign(   12 )
categ.SetTextSize ( 0.06 )
categ.SetTextColor(    1 )
categ.SetTextFont (   42 )
if i+1==1:       
  categ.AddText("OS")
elif i+1==2:       
  categ.AddText("SS")
#categ.Draw()

c.cd()
c.cd()
pad1.Draw()

ROOT.gPad.RedrawAxis()

c.Modified()
#c.SaveAs("plots/"+firstarg+"_fr.pdf")
c.SaveAs("plots/"+firstarg+secondarg+".png")


