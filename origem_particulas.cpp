// origem_particulas.cpp
// ---------------------------------------------------------------------------
// Pergunta: cada particula final descende do REMANESCENTE do feixe (status 63),
// do ESPALHAMENTO DURO (status 23), dos dois, ou de nenhum?
//
// Se a maioria cair em "so remanescente" ou "so processo duro", a separacao
// entre feixe e produzido e' possivel.
// Se a maioria cair em "os dois", a separacao nao existe.
// ---------------------------------------------------------------------------

#include "TSystem.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "Pythia8/Pythia.h"
#include <vector>

void origem_particulas(Int_t nev = 100)
{
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");
  pythia8->Initialize(2212, 2212, 14000);

  Pythia8::Pythia *py = pythia8->Pythia8();

  TH1F *h = new TH1F("h", "De onde descende cada particula final;;fracao", 4, 0, 4);
  h->GetXaxis()->SetBinLabel(1, "so remanescente");
  h->GetXaxis()->SetBinLabel(2, "so processo duro");
  h->GetXaxis()->SetBinLabel(3, "os dois");
  h->GetXaxis()->SetBinLabel(4, "nenhum");

  Long64_t total = 0;
  Long64_t cat[4] = {0, 0, 0, 0};

  for (Int_t iev = 0; iev < nev; iev++)
  {
    pythia8->GenerateEvent();

    for (Int_t i = 0; i < py->event.size(); i++)
    {
      if (!py->event[i].isFinal()) continue;
      total++;

      bool tem63 = false, tem23 = false;

      // percorre TODOS os ancestrais (seguindo as duas maes)
      std::vector<int> fila;
      fila.push_back(i);
      Int_t visitados = 0;

      while (!fila.empty() && visitados < 300)
      {
        Int_t idx = fila.back();
        fila.pop_back();
        visitados++;

        if (idx <= 0 || idx >= py->event.size()) continue;

        Int_t st = TMath::Abs(py->event[idx].status());
        if (st == 63) tem63 = true;   // remanescente do feixe
        if (st == 23) tem23 = true;   // saida do espalhamento duro

        Int_t m1 = py->event[idx].mother1();
        Int_t m2 = py->event[idx].mother2();

        if (m2 > m1 && m1 > 0) {
          // faixa de maes (tipico da hadronizacao: a string inteira)
          for (Int_t k = m1; k <= m2 && k - m1 < 50; k++) fila.push_back(k);
        } else {
          if (m1 > 0) fila.push_back(m1);
          if (m2 > 0 && m2 != m1) fila.push_back(m2);
        }
      }

      if      ( tem63 && !tem23) cat[0]++;
      else if (!tem63 &&  tem23) cat[1]++;
      else if ( tem63 &&  tem23) cat[2]++;
      else                       cat[3]++;
    }
  }

  const char *nomes[4] = {"so remanescente ", "so processo duro", "os dois         ", "nenhum          "};

  printf("\n===== Origem das particulas finais (%lld particulas, %d eventos) =====\n", total, nev);
  for (Int_t k = 0; k < 4; k++) {
    Double_t frac = (total > 0) ? (Double_t)cat[k] / total : 0;
    printf("  %s : %8lld   (%.1f%%)\n", nomes[k], cat[k], 100 * frac);
    h->SetBinContent(k + 1, frac);
  }
  printf("======================================================================\n");
  printf("Se 'os dois' domina, a separacao feixe/produzido nao e' bem definida.\n");

  TCanvas *c = new TCanvas("c", "Origem das particulas finais", 900, 600);
  h->SetFillColor(38);
  h->SetStats(0);
  h->Draw("BAR");
  c->SaveAs("origem_particulas.png");
}
