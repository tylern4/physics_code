"""
Usage: script [old_file.root] [new_file.root]

This function compares an old root file to a new root file, it checks the contents 
of the waveform data and the event header and outputs an error with information of
the entry number.  

The event header check uses the operator== of the class, which checks to ensure 
they are exactly the same.

The waveform data check uses EXOWaveformData::ContainsAllWaveformsIn() which checks
to ensure that the old waveform data contains all the waveforms in the new waveform data.
Alternatively one could call the operator== of the class in which case they must have 
the equal number of waveforms.
"""
import ROOT
import sys

if len(sys.argv) != 3:
    print "Usage: script [old_file.root] [new_file.root]"
    sys.exit(1)

ROOT.gSystem.Load("libEXOUtilities")
ROOT.gROOT.SetBatch()

# Make sure we can see debug output
ROOT.EXOErrorLogger.GetLogger().SetOutputThreshold(0)

old_file = ROOT.TFile(sys.argv[1])
new_file = ROOT.TFile(sys.argv[2])

for tree, branch in [("glitch", "GlitchBranch"),
                     ("veto", "VetoBranch")]: 
    old_tmp_tree = old_file.Get(tree)
    new_tmp_tree = new_file.Get(tree)
    if old_tmp_tree.GetEntries() != new_tmp_tree.GetEntries():
        print "Trees (%s) do not have equal numbers of entries" % tree
        sys.exit(1)
    for i in range(old_tmp_tree.GetEntries()):
        old_tmp_tree.GetEntry(i)
        new_tmp_tree.GetEntry(i)
        if not ( getattr(old_tmp_tree, branch) == getattr(new_tmp_tree, branch) ):
            print "Branches (%s) do not match at entry (%i)" % (branch, i)
    print "Done: %s:%s" % (tree, branch)



old_tree = old_file.Get("tree")
new_tree = new_file.Get("tree")

if old_tree.GetEntries() != new_tree.GetEntries():
    print "Trees do not have equal numbers of entries"
    sys.exit(1)

total_entries = old_tree.GetEntries()
c1 = ROOT.TCanvas()
name_of_file = "~/Dropbox/BadEvents.ps"
c1.Print(name_of_file + "[")
for i in range(total_entries):

    if i % 1000 == 0:
        print "Done %i of %i" % (i, total_entries)

    old_tree.GetEntry(i)
    new_tree.GetEntry(i)

    ed_old = old_tree.EventBranch
    ed_new = new_tree.EventBranch

    wf_data_old = ed_old.GetWaveformData()
    wf_data_new = ed_new.GetWaveformData()
    # Can also compare if the waveform data shouldn't have changed.
    #if not wf_data_old == wf_data_new :

    # compare waveforms, the old waveform data just needs to contain the new one, 
    # since the new bin package might throw away more waveforms.
    if not wf_data_old.ContainsAllWaveformsIn( wf_data_new ) :
        sys.stderr.write("Waveform data does not match at entry: %i\n" % i)
        wf_data_new.Decompress()
        wf_data_old.Decompress()
        for wf in wf_data_new.GetWaveformArray():
            old_wf = wf_data_old.GetWaveformWithChannel(wf.fChannel)
            new_hist = wf.GimmeHist("new")
            new_hist.SetLineColor(ROOT.kRed)
            new_hist.Draw("L")
            if old_wf:
                old_wf.GimmeHist("old").Draw("same L")
            c1.Update()
            c1.Print(name_of_file)
        break

    # Can output which waveforms don't exist
    channels_dropped = [str(wf.fChannel) for wf in wf_data_old.GetWaveformArray()
                                     if not wf_data_new.GetWaveformWithChannel(wf.fChannel)] 

    if len(channels_dropped) > 0:
        print "Channels dropped (Entry %i): " % (i)  + ' '.join(channels_dropped)
    
    

    # Check the event header 
    if not ed_old.fEventHeader == ed_new.fEventHeader:
        sys.stderr.write("Event Header does not match at entry: %i\n" % i)

c1.Print(name_of_file + "]")