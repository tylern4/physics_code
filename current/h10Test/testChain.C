{
    TChain* c = new TChain("h10");
    TFileCollection fc("fileList","","3.lis");
    c->AddFileInfoList((TCollection*)fc.GetList());
    TProof::Open("workers=3");
    c->SetProof();
    TStopwatch sw;
    sw.Start();
    c->Process("h10.C++"); //,"",1000000);
    sw.Stop();
    sw.Print();
}
