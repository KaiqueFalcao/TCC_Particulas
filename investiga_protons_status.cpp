// investiga_protons.cpp
// ---------------------------------------------------------------------------
// Compara duas formas de definir a "energia final dos feixes" para o Delta E:
//
//   METODO A (cinematico)  : barion (p/n) mais energetico de cada sentido.
//   METODO B (status code) : sobe a cadeia de maes de cada particula final
//                            procurando um ancestral com status 63
//                            (outgoing beam remnant, no registro do Pythia).
//
// O metodo B so e' possivel acessando o registro NATIVO do Pythia. A interface
// TParticle do ROOT converte o status para 0/1 (nao-final / final) e descarta
// o codigo original -- isso e' mostrado no diagnostico.
//
// Os dois primeiros paineis sao a EVIDENCIA de por que o metodo B nao serve:
//   - distribuicao dos status reais (mostra que eles existem e sao acessiveis)
//   - fracao das particulas finais classificadas como "vindas do remanescente"
//     (se essa fracao fica perto de 1, o criterio nao discrimina nada)
// ---------------------------------------------------------------------------

#include "TSystem.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "Pythia8/Pythia.h"

void investiga_protons(Int_t nev = 10000)
{
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia8");

  Double_t E_inicial = 14000;   // energia total dos dois feixes (GeV)

  // --------- METODO A: barion lider (cinematico) ---------
  TH2F *hist_mult = new TH2F("hist_mult","A) Multiplicidade vs #DeltaE (barion lider);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec   = new TH2F("hist_ec","A) Energia cinetica vs #DeltaE (barion lider);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // --------- METODO B: status 63 + genealogia ---------
  TH2F *hist_mult_st = new TH2F("hist_mult_st","B) Multiplicidade vs #DeltaE (status 63);#DeltaE (GeV);N particulas",100, 0, 14000, 100, 0, 800);
  TH2F *hist_ec_st   = new TH2F("hist_ec_st","B) Energia cinetica vs #DeltaE (status 63);#DeltaE (GeV);Soma Ecin (GeV)",100, 0, 14000, 100, 13800, 14000);

  // --------- Diagnosticos ---------
  TH1F *hist_status = new TH1F("hist_status","Status codes REAIS do Pythia (primeiros 100 eventos);status;particulas", 201, -100.5, 100.5);
  TH1F *hist_status_root = new TH1F("hist_status_root","Status que a TParticle do ROOT entrega;status;particulas", 21, -10.5, 10.5);
  TH1F *hist_fracao = new TH1F("hist_fracao","Fracao das finais classificadas como 'do remanescente';fracao;eventos", 102, -0.01, 1.01);

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");
  pythia8->ReadString("Random:seed = 43");
  pythia8->Initialize(2212, 2212, 14000);

  // Ponteiro para o objeto Pythia interno: da acesso ao registro nativo,
  // com os status codes de verdade.
  Pythia8::Pythia *py = pythia8->Pythia8();

  TClonesArray *particles = new TClonesArray("TParticle", 10000);

  Double_t soma_fracao = 0;

  for (Int_t iev = 0; iev < nev; iev++)
  {
    Double_t soma_ec = 0;
    Int_t multiplicidade = 0;
    Double_t E_lead_fwd = 0, E_lead_bwd = 0;

    pythia8->GenerateEvent();

    // ===================================================================
    // LOOP A -- via TParticle (ROOT). E' o caminho que voce ja usava.
    // ===================================================================
    pythia8->ImportParticles(particles, "All");
    Int_t np = particles->GetEntriesFast();

    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle *)particles->At(ip);

      // registra o que o ROOT entrega como status (amostra)
      if (iev < 100) hist_status_root->Fill(part->GetStatusCode());

      if (part->GetStatusCode() <= 0) continue;

      Double_t E   = part->Energy();
      Double_t m   = part->GetMass();
      Int_t    pdg = part->GetPdgCode();
      Double_t pz  = part->Pz();

      multiplicidade += 1;
      soma_ec += E - m;

      if (pdg == 2212 || pdg == 2112) {
        if (pz > 0 && E > E_lead_fwd) E_lead_fwd = E;
        if (pz < 0 && E > E_lead_bwd) E_lead_bwd = E;
      }
    }

    Double_t deltaE = E_inicial - (E_lead_fwd + E_lead_bwd);

    // ===================================================================
    // LOOP B -- via registro NATIVO do Pythia (status codes reais).
    // ===================================================================
    Double_t E_remanescente = 0;
    Int_t n_finais = 0, n_do_remanescente = 0;

    for (Int_t i = 0; i < py->event.size(); i++)
    {
      if (iev < 100) hist_status->Fill(py->event[i].status());

      if (!py->event[i].isFinal()) continue;
      n_finais++;

      // sobe a cadeia de maes procurando um ancestral marcado como remanescente
      Int_t idx = i, passos = 0;
      bool vem_do_remanescente = false;
      while (idx > 0 && passos < 200)
      {
        if (TMath::Abs(py->event[idx].status()) == 63) { vem_do_remanescente = true; break; }
        idx = py->event[idx].mother1();
        passos++;
      }

      if (vem_do_remanescente) {
        E_remanescente += py->event[i].e();
        n_do_remanescente++;
      }
    }

    Double_t deltaE_st = E_inicial - E_remanescente;
    Double_t fracao = (n_finais > 0) ? (Double_t)n_do_remanescente / n_finais : 0;

    soma_fracao += fracao;
    hist_fracao->Fill(fracao);

    if (E_lead_fwd > 0 && E_lead_bwd > 0) {
      hist_mult->Fill(deltaE, multiplicidade);
      hist_ec->Fill(deltaE, soma_ec);
    }

    hist_mult_st->Fill(deltaE_st, multiplicidade);
    hist_ec_st->Fill(deltaE_st, soma_ec);

    if (iev < 10)
      printf("Evento %2d:  A) deltaE=%8.1f   B) deltaE=%8.1f   fracao_remanescente=%.3f   mult=%3d\n",
             iev, deltaE, deltaE_st, fracao, multiplicidade);
  }

  printf("\n================= RESUMO =================\n");
  printf("A) barion lider : <deltaE> = %8.1f GeV  (RMS %.1f)\n",
         hist_mult->GetMean(1), hist_mult->GetRMS(1));
  printf("B) status 63    : <deltaE> = %8.1f GeV  (RMS %.1f)\n",
         hist_mult_st->GetMean(1), hist_mult_st->GetRMS(1));
  printf("\nFracao media das particulas finais classificadas\n");
  printf("como descendentes do remanescente: %.3f\n", soma_fracao / nev);
  printf("(se essa fracao esta perto de 1, o criterio nao separa nada)\n");
  printf("==========================================\n");

  TCanvas *c1 = new TCanvas("c1", "Delta E: cinematico vs status", 1600, 1400);
  c1->Divide(2, 3);

  c1->cd(1); gPad->SetLogy(); hist_status->Draw();
  c1->cd(2); gPad->SetLogy(); hist_status_root->Draw();
  c1->cd(3); hist_fracao->Draw();
  c1->cd(4); hist_mult->Draw("COLZ");
  c1->cd(5); hist_mult_st->Draw("COLZ");
  c1->cd(6); hist_ec_st->Draw("COLZ");

  c1->SaveAs("deltaE_status.png");
}