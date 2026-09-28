// Jatos e funcao de fragmentacao do pion neutro.
// Pythia gera, FastJet agrupa em jatos, e medimos z = pT(pi0)/pT(jato).
// O pi0 e' mantido estavel, entao ele proprio e' constituinte do jato.

#include "TSystem.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TPythia8.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include <vector>

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina, entao da' para mexer nos histogramas.
TH1F   *hist_pt_jato = 0;
TH1F   *hist_n_part  = 0;
TH1F   *hist_pt_pi0  = 0;
TH1F   *hist_z       = 0;
TH2F   *hist_eta_phi = 0;
TCanvas *c1 = 0;
TCanvas *c2 = 0;

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

void jatos_pi0_root()
{
  // ---- parametros 
  Int_t    nev          = 5000;
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 20.;
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;
  Int_t    semente      = 43;
  Double_t pt_hat_min   = 20.;   // sem isso quase nao sai jato acima de 20 GeV

  // ---- Pythia 
  char comando[100];

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("111:mayDecay = off");   // pi0 estavel: vira constituinte
  pythia8->ReadString("Random:setSeed = on");

  sprintf(comando, "PhaseSpace:pTHatMin = %f", pt_hat_min);
  pythia8->ReadString(comando);
  sprintf(comando, "Random:seed = %d", semente);
  pythia8->ReadString(comando);

  pythia8->Initialize(2212, 2212, energia);

  TClonesArray *particulas = new TClonesArray("TParticle", 10000);

  // ---- FastJet -------------------------------------------------------------
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  std::vector<fastjet::PseudoJet> entradas;

  // ---- histogramas ---------------------------------------------------------
  hist_pt_jato = new TH1F("hist_pt_jato",
      "pT dos jatos;pT^{jato} [GeV/c];jatos", 100, 0, 150);
  hist_n_part  = new TH1F("hist_n_part",
      "particulas por jato;N_{const};jatos", 60, 0, 60);
  hist_pt_pi0  = new TH1F("hist_pt_pi0",
      "pT dos #pi^{0};pT^{#pi^{0}} [GeV/c];#pi^{0}", 100, 0, 40);
  hist_z       = new TH1F("hist_z",
      "funcao de fragmentacao;z = pT^{#pi^{0}}/pT^{jato};#pi^{0}", 50, 0, 1);
  hist_eta_phi = new TH2F("hist_eta_phi",
      "mapa #eta-#phi dos #pi^{0};#eta;#phi [rad]",
      60, -eta_part_max, eta_part_max, 60, -TMath::Pi(), TMath::Pi());

  Long64_t n_jatos = 0, n_pi0 = 0;

  // ---- laco de eventos -----------------------------------------------------
  for (Int_t iev = 0; iev < nev; iev++)
  {
    if (iev % 500 == 0) printf("  evento %d / %d\n", iev, nev);

    pythia8->GenerateEvent();
    pythia8->ImportParticles(particulas, "All");
    Int_t np = particulas->GetEntriesFast();

    // monta a lista de entrada do FastJet
    entradas.clear();
    for (Int_t ip = 0; ip < np; ip++)
    {
      TParticle *part = (TParticle*) particulas->At(ip);

      if (part->GetStatusCode() <= 0) continue;
      if (eh_neutrino(part->GetPdgCode())) continue;
      if (part->Pt() < 1e-9) continue;
      if (TMath::Abs(part->Eta()) > eta_part_max) continue;

      fastjet::PseudoJet pj(part->Px(), part->Py(), part->Pz(), part->Energy());
      pj.set_user_index(ip);     // ponte de volta para o TClonesArray
      entradas.push_back(pj);
    }
    if (entradas.size() < 2) continue;

    // agrupa: sempre com pT minimo, senao cada particula mole vira um "jato"
    fastjet::ClusterSequence cluster_seq(entradas, jet_def);
    std::vector<fastjet::PseudoJet> jatos =
        fastjet::sorted_by_pt(cluster_seq.inclusive_jets(pt_jato_min));

    // ---- laco de jatos -----------------------------------------------------
    Int_t n_jatos_ev = jatos.size();
    for (Int_t ij = 0; ij < n_jatos_ev; ij++)
    {
      if (TMath::Abs(jatos[ij].rap()) > y_jato_max) continue;

      Double_t pt_jato = jatos[ij].pt();
      n_jatos++;

      std::vector<fastjet::PseudoJet> constituintes = jatos[ij].constituents();
      Int_t n_const = constituintes.size();

      hist_pt_jato->Fill(pt_jato);
      hist_n_part ->Fill(n_const);

      // ---- laco de constituintes: procura os pi0 ---------------------------
      for (Int_t ic = 0; ic < n_const; ic++)
      {
        Int_t indice = constituintes[ic].user_index();
        TParticle *part = (TParticle*) particulas->At(indice);

        if (part->GetPdgCode() != 111) continue;   // so pi0

        Double_t pt_pi0 = part->Pt();
        Double_t z      = pt_pi0 / pt_jato;
        n_pi0++;

        // TParticle::Phi() devolve [0, 2pi); o eixo do mapa vai de -pi a +pi
        Double_t phi = part->Phi();
        if (phi > TMath::Pi()) phi -= 2 * TMath::Pi();

        hist_pt_pi0 ->Fill(pt_pi0);
        hist_z      ->Fill(z);
        hist_eta_phi->Fill(part->Eta(), phi);
      }
    }
  }

  // ---- canvas 1: quatro paineis --------------------------------------------
  c1 = new TCanvas("c1", "Jatos e fragmentacao do pi0", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato->Draw();
  c1->cd(2); hist_n_part ->Draw();
  c1->cd(3); hist_pt_pi0 ->Draw();
  // ADICAO 1: eixo Y em escala log. A funcao de fragmentacao cai varias
  // ordens de grandeza, entao em escala linear so' se enxerga o primeiro bin.
  c1->cd(4); gPad->SetLogy(); hist_z->Draw();
  c1->Update();                      // desenha na tela agora
  c1->SaveAs("jatos_pi0.png");

  // ---- canvas 2: mapa eta-phi ----------------------------------------------
  c2 = new TCanvas("c2", "Mapa eta-phi dos pi0", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("mapa_eta_phi.png");

  // ---- resumo --------------------------------------------------------------
  printf("\n===================== RESUMO =====================\n");
  printf("  eventos gerados           : %d\n", nev);
  printf("  jatos reconstruidos       : %lld\n", n_jatos);
  printf("  pi0 dentro de jatos       : %lld\n", n_pi0);
  printf("  pi0 por jato              : %.3f\n",
         n_jatos ? (Double_t) n_pi0 / n_jatos : 0.);
  printf("  <pT> dos jatos            : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <z> dos pi0               : %.4f\n", hist_z->GetMean());

}
