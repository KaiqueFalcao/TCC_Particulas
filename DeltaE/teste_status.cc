// teste_status.cc
// Pythia nativo + histogramas do ROOT.
// Objetivo: usar os status codes reais e a genealogia para identificar
// as particulas FINAIS que descendem do remanescente do feixe (status 63).

#include "Pythia8/Pythia.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace Pythia8;

int main() {

  Pythia pythia;
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 43");
  pythia.readString("Beams:eCM = 14000.");
  pythia.init();

  double E_inicial = 14000.;

  TH2F *hist_mult = new TH2F("hist_mult",
      "Multiplicidade vs Delta E (genealogia);Delta E (GeV);N particulas",
      100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec = new TH2F("hist_ec",
      "Energia cinetica vs Delta E (genealogia);Delta E (GeV);Soma Ecin (GeV)",
      100, 0, 14000, 100, 13850, 14000);

  int nev = 10000;

  for (int iev = 0; iev < nev; ++iev) {

    if (!pythia.next()) continue;

    double soma_ec = 0., E_feixe = 0.;
    int multiplicidade = 0;

    // --- diagnostico: no primeiro evento, mostra os status de verdade ---
    if (iev == 0) {
      printf("Status codes reais no primeiro evento (primeiras 30 entradas):\n");
      for (int i = 0; i < 30 && i < pythia.event.size(); ++i)
        printf("  i=%2d  id=%6d  status=%4d\n",
               i, pythia.event[i].id(), pythia.event[i].status());
      printf("\n");
    }

    for (int i = 0; i < pythia.event.size(); ++i) {

      if (!pythia.event[i].isFinal()) continue;   // so estado final

      double E = pythia.event[i].e();
      double m = pythia.event[i].m();

      multiplicidade += 1;
      soma_ec += E - m;

      // --- sobe a arvore genealogica procurando um ancestral com status 63 ---
      bool vem_do_remanescente = false;
      int idx = i;
      int passos = 0;
      while (idx > 0 && passos < 200) {          // limite de passos por seguranca
        if (abs(pythia.event[idx].status()) == 63) { vem_do_remanescente = true; break; }
        idx = pythia.event[idx].mother1();       // sobe para a mae
        passos++;
      }

      if (vem_do_remanescente) E_feixe += E;
    }

    double deltaE = E_inicial - E_feixe;

    hist_mult->Fill(deltaE, multiplicidade);
    hist_ec->Fill(deltaE, soma_ec);

    if (iev < 5)
      printf("Evento %d:  E_feixe=%.1f  deltaE=%.1f  mult=%d  Ecin=%.1f\n",
             iev, E_feixe, deltaE, multiplicidade, soma_ec);
  }

  TCanvas *c1 = new TCanvas("c1", "Delta E genealogia", 1400, 600);
  c1->Divide(2, 1);
  c1->cd(1); hist_mult->Draw("COLZ");
  c1->cd(2); hist_ec->Draw("COLZ");
  c1->SaveAs("deltaE_genealogia.png");

  return 0;
}