// Jatos e funcao de fragmentacao do pion neutro -- VERSAO COM O pi0 DECAINDO.
// Pythia gera, FastJet agrupa em jatos, e medimos z = pT(pi0)/pT(jato).
//

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
#include <set>     // MUDANCA 1: necessario para nao contar o mesmo pi0 duas vezes

// Ponteiros em escopo de arquivo: continuam acessiveis no prompt do ROOT
// depois que a macro termina, entao da' para mexer nos histogramas.
TH1F   *hist_pt_jato = 0;
TH1F   *hist_n_part  = 0;
TH1F   *hist_pt_pi0  = 0;
TH1F   *hist_z       = 0;
TH2F   *hist_eta_phi = 0;
TH1F   *hist_assimetria = 0;
TH1F   *hist_z_gama     = 0;
TCanvas *c1 = 0;
TCanvas *c2 = 0;
TCanvas *c3 = 0;

Bool_t eh_neutrino(Int_t pdg) // devolve true se o PDG for de neutrino
{
  Int_t a = TMath::Abs(pdg);
  return (a == 12 || a == 14 || a == 16 || a == 18);
}

// MUDANCA 2: funcao nova. Sobe a cadeia de maes a partir da particula i e
// devolve o indice do primeiro pi0 ancestral (-1 se nao houver).
//
// Convencao dos indices, ja' deslocada pelo ROOT no ImportParticles:
//   m1 <  0     -> nao tem mae, chegou no topo do registro
//   m2 == m1    -> "carbon copy" (mesma particula recopiada); continua subindo
//   m2 >  m1    -> FAIXA de maes = hadronizacao da string de cor. PARA AQUI,

Int_t mae_pi0(TClonesArray *lista, Int_t i)
{
  Int_t n = lista->GetEntriesFast();
  Int_t idx = i;

  for (Int_t passos = 0; passos < 200; passos++)   // trava contra laco infinito
  {
    if (idx < 0 || idx >= n) return -1;

    TParticle *p = (TParticle*) lista->At(idx);
    Int_t m1 = p->GetFirstMother();
    Int_t m2 = p->GetSecondMother();

    if (m1 < 0 || m1 >= n) return -1;   // topo do registro
    if (m2 > m1)           return -1;   // faixa de maes -> string de cor

    if (((TParticle*) lista->At(m1))->GetPdgCode() == 111) return m1;  // achou

    idx = m1;
  }
  return -1;
}

void jatos_pi0_decai()
{
  // ---- parametros 
  Int_t    nev          = 5000;
  Double_t R            = 0.4;
  Double_t pt_jato_min  = 5.;
  Double_t eta_part_max = 0.9;
  Double_t y_jato_max   = 0.9 - R;
  Double_t energia      = 14000.;

  // ---- Pythia 

  TPythia8 *pythia8 = new TPythia8();
  pythia8->ReadString("HardQCD:all = on");
  pythia8->ReadString("Random:setSeed = on");

  pythia8->ReadString("PhaseSpace:pTHatMin = 20."); // sem isso quase nao sai jato acima de 20 GeV
  pythia8->ReadString("Random:seed = 43");

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

  
  hist_assimetria = new TH1F("hist_assimetria",
      "assimetria do decaimento;A = |E_{1}-E_{2}|/(E_{1}+E_{2});#pi^{0}",
      50, 0, 1);
  hist_z_gama = new TH1F("hist_z_gama",
      "fragmentacao do foton lider;z_{#gamma} = pT^{#gamma}/pT^{jato};#gamma",
      50, 0, 1);

  Long64_t n_jatos = 0, n_pi0 = 0;
  Long64_t n_gama = 0;   // ADICAO 3: quantos pares gama gama foram analisados

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
      std::set<Int_t> indices_pi0;

      for (Int_t ic = 0; ic < n_const; ic++)
      {
        Int_t indice = constituintes[ic].user_index();
        Int_t indice_pi0 = mae_pi0(particulas, indice);

        if (indice_pi0 < 0) continue;   // esse constituinte nao veio de pi0

        indices_pi0.insert(indice_pi0);
      }

      // agora percorremos so' os pi0 unicos encontrados acima
      std::set<Int_t>::iterator it;
      for (it = indices_pi0.begin(); it != indices_pi0.end(); it++)
      {
        TParticle *part = (TParticle*) particulas->At(*it);

        Double_t pt_pi0 = part->Pt();
        Double_t z      = pt_pi0 / pt_jato;
        n_pi0++;

        // TParticle::Phi() devolve [0, 2pi); o eixo do mapa vai de -pi a +pi
        Double_t phi = part->Phi();
        if (phi > TMath::Pi()) phi -= 2 * TMath::Pi();

        hist_pt_pi0 ->Fill(pt_pi0);
        hist_z      ->Fill(z);
        hist_eta_phi->Fill(part->Eta(), phi);

        // O pi0 encontrado pelo mae_pi0() e' justamente a copia que decaiu,
        // entao as filhas dele sao os dois fotons.
        Int_t d1 = part->GetFirstDaughter();
        Int_t d2 = part->GetLastDaughter();

        if (d1 < 0 || d2 <= d1 || d2 >= np) continue;

        TParticle *gama1 = (TParticle*) particulas->At(d1);
        TParticle *gama2 = (TParticle*) particulas->At(d2);

        if (gama1->GetPdgCode() != 22) continue;
        if (gama2->GetPdgCode() != 22) continue;

        Double_t e1 = gama1->Energy();
        Double_t e2 = gama2->Energy();

        Double_t assimetria = TMath::Abs(e1 - e2) / (e1 + e2);
        hist_assimetria->Fill(assimetria);

        // funcao de fragmentacao do foton MAIS ENERGETICO dos dois
        TParticle *gama_lider = gama1;
        if (e2 > e1) gama_lider = gama2;

        Double_t z_gama = gama_lider->Pt() / pt_jato;
        hist_z_gama->Fill(z_gama);

        n_gama++;
      }
    }
  }

  // ---- canvas 1: quatro paineis --------------------------------------------
  c1 = new TCanvas("c1", "Jatos e fragmentacao do pi0", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); hist_pt_jato->Draw();
  c1->cd(2); hist_n_part ->Draw();
  c1->cd(3); hist_pt_pi0 ->Draw();
  c1->cd(4); hist_z->Draw();
  c1->Update();                 
  c1->SaveAs("jatos_pi0_decai.png");

  // ---- canvas 2: mapa eta-phi ----------------------------------------------
  c2 = new TCanvas("c2", "Mapa eta-phi dos pi0", 900, 700);
  hist_eta_phi->Draw("COLZ");
  c2->Update();
  c2->SaveAs("mapa_eta_phi_decai.png");

  // ---- canvas 3: fotons do decaimento --------------------------------------
  // ADICAO 6: canvas novo, para nao mexer no layout 2x2 do c1.
  c3 = new TCanvas("c3", "Fotons do decaimento do pi0", 1200, 500);
  c3->Divide(2, 1);
  c3->cd(1); hist_assimetria->Draw();                  // deve sair PLANA
  c3->cd(2); gPad->SetLogy(); hist_z_gama->Draw();     // tambem em log
  c3->Update();
  c3->SaveAs("fotons_decai.png");

  // ---- resumo --------------------------------------------------------------
  printf("\n===================== RESUMO =====================\n");
  printf("  eventos gerados           : %d\n", nev);
  printf("  jatos reconstruidos       : %lld\n", n_jatos);
  printf("  pi0 dentro de jatos       : %lld\n", n_pi0);
  printf("  pi0 por jato              : %.3f\n",
         n_jatos ? (Double_t) n_pi0 / n_jatos : 0.);
  printf("  <pT> dos jatos            : %.2f GeV/c\n", hist_pt_jato->GetMean());
  printf("  <z> dos pi0               : %.4f\n", hist_z->GetMean());
  // ADICAO 7: linhas novas no resumo
  printf("  pares gama gama           : %lld\n", n_gama);
  printf("  <A> assimetria            : %.4f   (esperado ~0,5 se plana)\n",
         hist_assimetria->GetMean());
  printf("  RMS da assimetria         : %.4f   (esperado ~0,289 = 1/sqrt12)\n",
         hist_assimetria->GetRMS());
  printf("  <z> do foton lider        : %.4f\n", hist_z_gama->GetMean());
}
