
#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h" // ADICIONADO
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"

void investiga_protons(Int_t nev = 10000)
{
  gSystem->Load("libEG");
  gSystem->Load("libEgPythia8");

  Double_t E_inicial = 14000;   // energia total dos dois feixes (GeV)

  // eixo x = Delta E ; eixo y = multiplicidade / energia cinetica
  TH2F *hist_mult = new TH2F("hist_mult","Multiplicidade vs Delta E;Delta E (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec = new TH2F("hist_ec","Energia cinetica vs Delta E;Delta E (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13000, 14000);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");
  pythia8->Initialize(2212, 2212, 14000);

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  for (Int_t iev = 0; iev < nev; iev++)
  {
    Double_t soma_ec = 0;
    Int_t multiplicidade = 0;
    Double_t E_lead_fwd = 0, E_lead_bwd = 0;   // proton lider para cada sentido, um pz positivo e outro negativo

    pythia8->GenerateEvent();
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle *)particles->At(ip);
      if (part->GetStatusCode() <= 0) continue;

      Double_t E   = part->Energy();
      Double_t m   = part->GetMass();
      Int_t    pdg = part->GetPdgCode();
      Double_t pz  = part->Pz();

      multiplicidade += 1;
      soma_ec += E - m;                    // energia cinetica da particula

      if (pdg == 2212 || pdg == 2112) {   // proton OU neutron guarda o mais energetico de cada lado
        if (pz > 0 && E > E_lead_fwd) E_lead_fwd = E;
        if (pz < 0 && E > E_lead_bwd) E_lead_bwd = E;
      }
    }

    Double_t E_final = E_lead_fwd + E_lead_bwd;   // energia dos dois protons lideres
    Double_t deltaE  = E_inicial - E_final;       // energia que o feixe "perdeu"

    hist_mult->Fill(deltaE, multiplicidade);
    hist_ec->Fill(deltaE, soma_ec);
  }

  TCanvas *c1 = new TCanvas("c1", "Delta E", 1400, 600);
  c1->Divide(2, 1);

  c1->cd(1);
  hist_mult->Draw("COLZ");

  c1->cd(2);
  hist_ec->Draw("COLZ");

  c1->SaveAs("deltaE.png");
}