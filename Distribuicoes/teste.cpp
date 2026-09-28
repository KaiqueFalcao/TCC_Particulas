
#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"

void distribuicoes_de_energia_e_momento(Int_t nev = 10000, Int_t ndeb = 1 /* Listagem */ )
{
  gSystem->Load("libEG");
  gSystem->Load("libEgPythia8");

  Double_t soma_de_energia;
  Double_t soma_de_momento;
  Double_t px, py;
  Double_t energia_part_final;
  Double_t pt, cont_pt = 0, cont_pt2 = 0;

  TH1F *hist_soma_energia_final = new TH1F("hist_soma_energia_final", "Soma das energias das particulas de estado final do evento", 100, 1, 1);
  TH1F *hist_momentos_particulas_finais = new TH1F("hist_pt_part_final", "Distribuicao dos momentos transversais das particulas de estado final", 100, 1, 1);
  TH1F *hist_px_part_final = new TH1F("hist_px_part_final", "Momentos transversais das particulas finais no eixo X", 100, 1, 1);
  TH1F *hist_py_part_final = new TH1F("hist_py_part_final", "Momentos transversais das particulas finais no eixo Y", 100, 1, 1);
  TH1F *hist_soma_momento = new TH1F("hist_soma_momento","Soma dos momentos das particulas finais dos eventos", 100, 1, 1);
  TH2F *dispersao_momentos = new TH2F("dispersao_momentos", "Dispersao dos momentos de estado central e final", 100, 1, 1, 100, 1, 1);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");

  pythia8->Initialize(2212 /* Proton */, 2212 /* Proton */, 14000 /* TeV */); /* 14000 TeV = 14000000 GeV */

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  for (Int_t iev = 0; iev < nev; iev++)
  {
    cont_pt = 0;
    cont_pt2 = 0;
    soma_de_momento = 0;
    pythia8->GenerateEvent();

    if (iev < ndeb) pythia8->EventListing();

    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particles->At(ip); /* Criacao do ponteiro */

      Int_t status = part->GetStatusCode();
      Double_t eta = part->Eta();
      Double_t phi = part->Phi();

      if (status > 0) /* Particula de estado final */
      {
        pt = part->Pt();
        px = part->Px();
        py = part->Py();

        energia_part_final = part->Energy();

        if (TMath::Cos(phi) * TMath::Sin(phi) < 0 )
        {
          pt *= -1;
        }

        soma_de_momento += pt;
        soma_de_energia += energia_part_final;
        cont_pt += 1;
        hist_momentos_particulas_finais->Fill(pt);
        hist_px_part_final->Fill(px);
        hist_py_part_final->Fill(py);

        if ( -0.9 < eta && eta < 0.9 )
        {
          cont_pt2 += 1;
        }

      }
    }
    dispersao_momentos->Fill(cont_pt, cont_pt2);
    hist_soma_momento->Fill(soma_de_momento);
  }

  hist_soma_energia_final->Fill(soma_de_energia);

  TCanvas *c1 = new TCanvas("c1", "Histogramas e distribuicao", 2500, 2500);
  c1->Divide(2, 3);

  c1->cd(1);
  hist_soma_energia_final->SetTitle("Soma das energias das particulas de estado final do evento");
  hist_soma_energia_final->GetXaxis()->SetTitle("Soma");
  hist_soma_energia_final->GetYaxis()->SetTitle("Frequencia");
  hist_soma_energia_final->Draw();

  c1->cd(2);
  hist_momentos_particulas_finais->SetTitle("Distribuicao de momento transversal das particulas de estado final do evento");
  hist_momentos_particulas_finais->GetXaxis()->SetTitle("Momento");
  hist_momentos_particulas_finais->GetYaxis()->SetTitle("Frequencia");
  hist_momentos_particulas_finais->Draw();

  c1->cd(3);
  hist_px_part_final->SetTitle("Distribuicao dos momentos em X das particulas finais");
  hist_px_part_final->GetXaxis()->SetTitle("Px");
  hist_px_part_final->GetYaxis()->SetTitle("Frequencia");
  hist_px_part_final->Draw();

  c1->cd(4);
  hist_py_part_final->SetTitle("Distribuicao dos momentos em Y das particulas finais");
  hist_py_part_final->GetXaxis()->SetTitle("Py");
  hist_py_part_final->GetYaxis()->SetTitle("Frequencia");
  hist_py_part_final->Draw();

  c1->cd(5);
  hist_soma_momento->SetTitle("Soma momentos das particulas finais");
  hist_soma_momento->GetXaxis()->SetTitle("Soma");
  hist_soma_momento->GetYaxis()->SetTitle("Frequencia");
  hist_soma_momento->Draw();

  c1->cd(6);
  dispersao_momentos->SetTitle("Dispersao dos momentos de estado central e final");
  dispersao_momentos->GetXaxis()->SetTitle("Momento 1");
  dispersao_momentos->GetYaxis()->SetTitle("Momento 2");
  dispersao_momentos->Draw();



}