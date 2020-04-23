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
histoname = firstarg+secondarg
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

if (firstarg == "h_dPhi_" or firstarg =="Tau_phi_" or firstarg =="Electron_phi_" or firstarg =="dphi_muMet_" or firstarg =="dphi_tauMet"):
  piAxis = 1
else :
  piAxis = 0



def add_lumi():
    lowX=0.58
    lowY=0.835
    lumi  = ROOT.TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.06)
    lumi.SetTextFont (   42 )
    lumi.AddText("2016, 35.9 fb^{-1} (13 TeV)")
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
        output = ROOT.TLegend(0.70, 0.40, 0.92, 0.84, "", "brNDC")
        #output = ROOT.TLegend(0.2, 0.1, 0.47, 0.65, "", "brNDC")
        output.SetLineWidth(0)
        output.SetLineStyle(0)
        output.SetFillStyle(0)
        output.SetBorderSize(0)
        output.SetTextFont(62)
        return output

ROOT.gStyle.SetFrameLineWidth(3)
ROOT.gStyle.SetLineWidth(3)
ROOT.gStyle.SetOptStat(0)

c=ROOT.TCanvas("canvas","",0,0,600,600)
c.cd()

OutFile=ROOT.TFile("f_etau_initial.root","r")

adapt=ROOT.gROOT.GetColor(12)
new_idx=ROOT.gROOT.GetListOfColors().GetSize() + 1
trans=ROOT.TColor(new_idx, adapt.GetRed(), adapt.GetGreen(),adapt.GetBlue(), "",0.5)

categories=[firstarg,firstarg] 
dirName=[firstarg+"/","qcd/","bkg/"+firstarg+"/"] 
ncat=1


for i in range (0,ncat):
    HistoData=OutFile.Get(dirName[0]+"data_obs_"+histoname)
    ZTT=OutFile.Get(dirName[0]+"ZTT_"+histoname)
    VV=OutFile.Get(dirName[0]+"VV_"+histoname)
    ZTT1=OutFile.Get(dirName[0]+"ZTT1_"+histoname)
    ZTT2=OutFile.Get(dirName[0]+"ZTT2_"+histoname)
    ZTT3=OutFile.Get(dirName[0]+"ZTT3_"+histoname)
    ZTT4=OutFile.Get(dirName[0]+"ZTT4_"+histoname)
    ZTT.Add(ZTT1)
    ZTT.Add(ZTT2)
    ZTT.Add(ZTT3)
    ZTT.Add(ZTT4)
  
    EWKZ=OutFile.Get(dirName[0]+"EWKZ_"+histoname)
    EWKW=OutFile.Get(dirName[0]+"EWKW_"+histoname)
    EWKZ.Add(EWKW)

    TTJ=OutFile.Get(dirName[0]+"TTJ_"+histoname)

    ggH125=OutFile.Get(dirName[0]+"GluGluHToTauTau_"+histoname)
    ggWW5=OutFile.Get(dirName[0]+"GluGluWWTo2L2Nu_"+histoname)
    qqH125=OutFile.Get(dirName[0]+"VBFHToTauTau_"+histoname)
    ZHToTauTau=OutFile.Get(dirName[0]+"ZHToTauTau_"+histoname)
    ggH125.Add(ZHToTauTau)
    ggH125.Add(ggWW5)
    ggH125.Add(qqH125)

    WJets=OutFile.Get(dirName[0]+"WJetsToLNu_2J_"+histoname)
    WJets_200To400=OutFile.Get(dirName[0]+"WJetsToLNu_200To400_"+histoname)
    WJets_400To600=OutFile.Get(dirName[0]+"WJetsToLNu_400To600_"+histoname)
    WJets_600To800=OutFile.Get(dirName[0]+"WJetsToLNu_600To800_"+histoname)
    WJets_800To1200=OutFile.Get(dirName[0]+"WJetsToLNu_800To1200_"+histoname)
    WJets_1200To2500=OutFile.Get(dirName[0]+"WJetsToLNu_1200To2500_"+histoname)
    WH125=OutFile.Get(dirName[0]+"WHToTauTau_"+histoname)
#    WPlusH125=OutFile.Get(dirName[0]+"WPlusHToTauTau_"+histoname)
    WJets.Add(WJets_200To400)
    WJets.Add(WJets_400To600)
    WJets.Add(WJets_600To800)
    WJets.Add(WJets_800To1200)
    WJets.Add(WJets_1200To2500)
    WJets.Add(WH125)
 #   WJets.Add(WPlusH125)


    WWW=OutFile.Get(dirName[0]+"WWW_"+histoname)
    WWZ=OutFile.Get(dirName[0]+"WWZ_"+histoname)
    WZZ=OutFile.Get(dirName[0]+"WZZ_"+histoname)
    ZZZ=OutFile.Get(dirName[0]+"ZZZ_"+histoname)
    WWW.Add(WWZ)
    WWW.Add(WZZ)
    WWW.Add(ZZZ)

    ZJetsToNuNu=OutFile.Get(dirName[0]+"ZJetsToNuNu"+histoname)



    ZTT.SetFillColor(ROOT.TColor.GetColor("#EA865F"))
    EWKZ.SetFillColor(ROOT.TColor.GetColor("#27AE60"))
    #QCD_mc.SetFillColor(ROOT.TColor.GetColor("#3422f7"))
    VV.SetFillColor(ROOT.TColor.GetColor("#1F618D"))
    TTJ.SetFillColor(ROOT.TColor.GetColor("#9B59B6"))
    ggH125.SetFillColor(ROOT.TColor.GetColor("#7B241C"))
    WJets.SetFillColor(ROOT.TColor.GetColor("#48C9B0"))
    WWW.SetFillColor(ROOT.TColor.GetColor("#F7DC6F"))
#    ZJetsToNuNu.SetFillColor(ROOT.TColor.GetColor("#76448A"))


    EWKZ.SetLineColor(1)
    ZTT.SetLineColor(1)
    VV.SetLineColor(1)
    TTJ.SetLineColor(1)
    ggH125.SetLineColor(1)
    WJets.SetLineColor(1)
    WWW.SetLineColor(1)
 #   ZJetsToNuNu.SetLineColor(1)
    #QCD_mc.SetLineColor(1)
   # ZH125_mc.SetLineColor(4)
  #  ZH125_mc.SetLineWidth(5)

    stack=ROOT.THStack("stack","stack")
#    stack.Add(ZJetsToNuNu)
    stack.Add(WWW)
    stack.Add(EWKZ)
    stack.Add(WJets)
    stack.Add(ggH125)
    stack.Add(TTJ)
    stack.Add(VV)
    stack.Add(ZTT)
 

    HistoData.GetXaxis().SetTitle("")
    HistoData.GetXaxis().SetTitleSize(0)


    if piAxis == 1:
      HistoData.GetXaxis().SetNdivisions(-405)
    else :
      HistoData.GetXaxis().SetNdivisions(-405)

    HistoData.GetXaxis().SetLabelSize(0)
    HistoData.GetYaxis().SetLabelFont(42)
    HistoData.GetYaxis().SetLabelOffset(0.01)
    HistoData.GetYaxis().SetLabelSize(0.05)
    HistoData.GetYaxis().SetTitleSize(0.075)
    HistoData.GetYaxis().SetTitleOffset(1.04)
    HistoData.SetTitle("")
    HistoData.GetYaxis().SetTitle("Events")


    errorBand=OutFile.Get(dirName[0]+"ZTT_"+histoname).Clone()
    errorBand.Add(VV)
    errorBand.Add(TTJ)
    errorBand.Add(ggH125)
    errorBand.Add(WJets)
    errorBand.Add(EWKZ)
    errorBand.Add(WWW)
 #   errorBand.Add(ZJetsToNuNu)


    errorBand.SetMarkerSize(0)
    errorBand.SetFillColor(1)
    errorBand.SetFillStyle(3001)
    errorBand.SetLineWidth(1)

pad1 = ROOT.TPad("pad1","pad1",0,0.35,1,1)
pad1.Draw()
pad1.cd()
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetBorderSize(10)
pad1.SetTickx(1)
pad1.SetTicky(1)
pad1.SetGridx()
pad1.SetLeftMargin(0.18)
pad1.SetRightMargin(0.05)
pad1.SetTopMargin(0.122)
pad1.SetBottomMargin(0.026)
pad1.SetFrameFillStyle(0)
pad1.SetFrameLineStyle(0)
pad1.SetFrameLineWidth(3)
pad1.SetFrameBorderMode(0)
pad1.SetFrameBorderSize(10)
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
HistoData.SetMarkerStyle(20)
HistoData.SetMarkerSize(1)
HistoData.GetXaxis().SetTitle(thirdarg)
HistoData.GetYaxis().SetTitle("Events")
if yaxisLog == 1 :
  HistoData.SetMaximum(1000*max(HistoData.GetMaximum(),stack.GetMaximum()))
  HistoData.SetMinimum(0.1)
else :
  HistoData.SetMaximum(1.35*max(HistoData.GetMaximum(),stack.GetMaximum()))
#stack.SetLineColor(9)
HistoData.Draw("e1")
stack.Draw("histsame")
errorBand.Draw("e2same")
#ZH125_mc.Draw("histsame")
HistoData.Draw("e1same")
errorBand.Draw("e2same")

legende=make_legend()
legende.AddEntry(HistoData,"Data","elp")
#legende.AddEntry(HistoData,"fake background","f")
#legende.AddEntry(W_mc,"W","f")
legende.AddEntry(ZTT,"Drell-Yan","f")
legende.AddEntry(TTJ,"TT jets","f")
legende.AddEntry(VV,"VV","f")
#legende.AddEntry(QCD_mc,"QCD","f")
legende.AddEntry(ggH125,"ggH125","f")
legende.AddEntry(WJets,"Wjets","f")
legende.AddEntry(EWKZ,"EWK","f")
legende.AddEntry(WWW,"WWW","f")
#legende.AddEntry(ZJetsToNuNu,"ZJetsToNuNu","f")

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
pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35);
pad2.SetTopMargin(0.05);
pad2.SetBottomMargin(0.35);
pad2.SetLeftMargin(0.18);
pad2.SetRightMargin(0.05);
pad2.SetTickx(1)
pad2.SetTicky(1)
pad2.SetFrameLineWidth(3)
pad2.SetGridx()
pad2.SetGridy()
pad2.Draw()
pad2.cd()
h1=HistoData.Clone()
h1.SetMaximum(2.0)#FIXME(1.5)
h1.SetMinimum(0.0)#FIXME(0.5)
h1.SetMarkerStyle(20)
h3=errorBand.Clone()
hwoE=errorBand.Clone()
for iii in range (1,hwoE.GetSize()-2):
  hwoE.SetBinError(iii,0)
h3.Sumw2()
h1.Sumw2()
h1.SetStats(0)
h1.Divide(hwoE)
h3.Divide(hwoE)
h1.GetXaxis().SetTitle(thirdarg)#(#vec{p_{T}}(#tau_{1})+#vec{p_{T}}(#tau_{2}))/(p_{T}(#tau_{1})+p_{T}(#tau_{2}))")#("m_{vis} (GeV)")#(#vec{p_{T}(#mu)}+#vec{p_{T}(#tau)})/(p_{T}(#mu)+p_{T}(#tau))")
#if (i+1==1 or i+1==2 or i+1==7 or i+1==8):
#	h1.GetXaxis().SetTitle("Electron p_{T} (GeV)")
#if (i+1==4 or i+1==10):
#     h1.GetXaxis().SetTitle("Muon p_{T} (GeV)")
#if (i+1==6 or i+1==12 or i+1==3 or i+1==5 or i+1==9 or i+1==11):
#     h1.GetXaxis().SetTitle("Tau p_{T} (GeV)")
h1.GetXaxis().SetLabelSize(0.1)
h1.GetYaxis().SetTitle("Obs./Exp.")
h1.GetYaxis().SetLabelSize(0.11)
h1.SetMaximum(2.0)#FIXME(1.5)
h1.SetMinimum(0.0)#FIXME(0.5)
if piAxis ==1 :
  if firstarg == "Tau_phi_" or firstarg == "Electron_phi_" :
    h1.GetXaxis().SetBinLabel(64,"#pi");
    h1.GetXaxis().SetBinLabel(48,"#frac{#pi}{2}");
    h1.GetXaxis().SetBinLabel(32,"0");
    h1.GetXaxis().SetBinLabel(16,"#frac{-#pi}{2}");
    h1.GetXaxis().SetBinLabel(1,"-#pi");
  else :
    h1.GetXaxis().SetBinLabel(40,"#pi");
    h1.GetXaxis().SetBinLabel(30,"#frac{3#pi}{4}");
    h1.GetXaxis().SetBinLabel(20,"#frac{#pi}{2}");
    h1.GetXaxis().SetBinLabel(10,"#frac{#pi}{4}");
    h1.GetXaxis().SetBinLabel(1,"0");

if piAxis == 1:
  h1.GetXaxis().SetLabelOffset(0.01)
  h1.GetXaxis().SetLabelSize(0.15)
  h1.GetXaxis().SetNdivisions(-405)
  
else :
  h1.GetXaxis().SetNdivisions(-405)

  
h1.GetYaxis().SetNdivisions(5)
h1.GetXaxis().SetTitleFont(42)
h1.GetYaxis().SetTitleFont(42)


h1.GetXaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleSize(0.15)
h1.GetYaxis().SetTitleOffset(0.56)
h1.GetXaxis().SetTitleOffset(1.1)






h1.Draw("e0p")
h3.Draw("e2same")

c.cd()
#  pad1.Draw()

ROOT.gPad.RedrawAxis()

c.Modified()
#c.SaveAs("plots/"+firstarg+"_fr.pdf")
c.SaveAs("plots/"+firstarg+secondarg+".png")


